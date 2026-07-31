/*
  ==============================================================================

    Standalone tests for TapeChannel.

    Nothing under Source/dsp depends on JUCE, so this builds and runs without the
    plugin wrapper:

        clang++ -std=c++20 -O2 -Wall -Wextra TapeTests.cpp ../dsp/TapeChannel.cpp \
                -o tapetests && ./tapetests

    Every requirement of the tape model is measured rather than asserted by
    construction: the level of each defect, where the frequency-dependent
    headroom actually falls, that the saturation cannot fold, and that the noise
    is reproducible, per channel, and free of allocation.

  ==============================================================================
*/

#include "../dsp/TapeChannel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <string>
#include <utility>
#include <vector>

using namespace lime;

//==============================================================================
// Allocation counting, so process() can be shown to be allocation free.
namespace
{
    bool countingAllocations = false;
    long allocationCount = 0;
}

void* operator new (std::size_t size)
{
    if (countingAllocations)
        ++allocationCount;

    void* p = std::malloc (size == 0 ? 1 : size);

    if (p == nullptr)
        throw std::bad_alloc();

    return p;
}

void* operator new[] (std::size_t size)
{
    return operator new (size);
}

void operator delete (void* p) noexcept { std::free (p); }
void operator delete[] (void* p) noexcept { std::free (p); }
void operator delete (void* p, std::size_t) noexcept { std::free (p); }
void operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

//==============================================================================
namespace
{

int failures = 0;
int checks = 0;

constexpr double twoPi = 6.283185307179586476925286766559;
constexpr double sampleRate = 48000.0;

/** Reports a level in dB relative to a reference magnitude, with an exact zero
    printed as a large negative number rather than -inf. */
double relativeDb (double error, double reference)
{
    if (reference <= 0.0)
        return error <= 0.0 ? -999.0 : 0.0;

    if (error <= 0.0)
        return -999.0;

    return 20.0 * std::log10 (error / reference);
}

void expectBelow (const std::string& what, double measuredDb, double limitDb)
{
    ++checks;
    const bool ok = measuredDb <= limitDb;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-56s %9.2f dB  (limit %7.2f dB)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measuredDb, limitDb);
}

void expectWithin (const std::string& what, double measured, double expected, double tolerance)
{
    ++checks;
    const bool ok = std::abs (measured - expected) <= tolerance;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-56s %9.3f    (want %7.3f +/- %.3f)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measured, expected, tolerance);
}

void expectTrue (const std::string& what, bool ok, const std::string& detail = {})
{
    ++checks;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-56s %s\n", ok ? "pass" : "FAIL", what.c_str(), detail.c_str());
}

//==============================================================================
double rms (const std::vector<double>& v, size_t from = 0)
{
    double sum = 0.0;

    for (size_t n = from; n < v.size(); ++n)
        sum += v[n] * v[n];

    const size_t count = v.size() - from;

    return count > 0 ? std::sqrt (sum / (double) count) : 0.0;
}

double peak (const std::vector<double>& v)
{
    double p = 0.0;

    for (auto x : v)
        p = std::max (p, std::abs (x));

    return p;
}

double maxDifference (const std::vector<double>& a, const std::vector<double>& b)
{
    double worst = 0.0;

    for (size_t n = 0; n < std::min (a.size(), b.size()); ++n)
        worst = std::max (worst, std::abs (a[n] - b[n]));

    return worst;
}

/** Amplitude of the component sitting exactly on bin `cycles` of a window of
    `count` samples. Exact for that bin, and orthogonal to every harmonic of it,
    so saturation products cannot contaminate the reading. */
double binAmplitude (const double* data, int count, int cycles)
{
    double re = 0.0, im = 0.0;

    for (int n = 0; n < count; ++n)
    {
        const double phase = twoPi * (double) cycles * (double) n / (double) count;
        re += data[n] * std::cos (phase);
        im -= data[n] * std::sin (phase);
    }

    return 2.0 * std::sqrt (re * re + im * im) / (double) count;
}

/** Welch estimate of the narrowband power at a frequency: mean squared amplitude
    over Hann-windowed segments. Used only for spectral shape, in dB relative to
    another such reading, so the analysis bandwidth cancels. */
double bandPower (const std::vector<double>& v, double freqHz, int segment)
{
    double total = 0.0;
    int segments = 0;

    for (size_t start = 0; start + (size_t) segment <= v.size(); start += (size_t) segment)
    {
        double re = 0.0, im = 0.0, windowSum = 0.0;

        for (int n = 0; n < segment; ++n)
        {
            const double w = 0.5 - 0.5 * std::cos (twoPi * (double) n / (double) segment);
            const double phase = twoPi * freqHz * (double) n / sampleRate;

            re += v[start + (size_t) n] * w * std::cos (phase);
            im -= v[start + (size_t) n] * w * std::sin (phase);
            windowSum += w;
        }

        const double amplitude = 2.0 * std::sqrt (re * re + im * im) / windowSum;

        total += amplitude * amplitude;
        ++segments;
    }

    return segments > 0 ? total / (double) segments : 0.0;
}

//==============================================================================
TapeChannelSettings allOff()
{
    TapeChannelSettings s;

    s.hissEnabled = false;
    s.saturationEnabled = false;
    s.hfLossEnabled = false;
    s.modulationNoiseEnabled = false;

    return s;
}

/** Runs one channel of a prepared TapeChannel over `input`, in blocks. */
std::vector<double> runMono (TapeChannel& tape, const std::vector<double>& input, int blockSize)
{
    std::vector<double> out = input;
    int n = 0;

    while (n < (int) out.size())
    {
        const int count = std::min (blockSize, (int) out.size() - n);
        double* pointers[1] = { out.data() + n };

        tape.process (pointers, 1, count);
        n += count;
    }

    return out;
}

std::vector<double> makeSine (int count, int cycles, double amplitude)
{
    std::vector<double> v ((size_t) count, 0.0);

    for (int n = 0; n < count; ++n)
        v[(size_t) n] = amplitude * std::sin (twoPi * (double) cycles * (double) n
                                                  / (double) count);

    return v;
}

std::vector<double> makeNoise (int count, uint64_t seed)
{
    std::vector<double> v ((size_t) count, 0.0);
    uint64_t s = seed;

    for (int n = 0; n < count; ++n)
    {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        v[(size_t) n] = 0.5 * (2.0 * ((double) (s >> 11) * 0x1.0p-53) - 1.0);
    }

    return v;
}

const double referenceAmplitude = tapeModel::referenceAmplitude;

//==============================================================================
void testTransparency()
{
    std::printf ("\nEverything disabled: the channel is a wire\n");

    TapeChannel tape;
    tape.prepare (sampleRate, 512, 2);
    tape.setSeed (1234);
    tape.setSettings (allOff());

    const auto input = makeNoise (4096, 99);
    std::vector<double> left = input, right = input;
    double* pointers[2] = { left.data(), right.data() };

    tape.process (pointers, 2, (int) left.size());

    expectBelow ("null vs input, channel 0",
                 relativeDb (maxDifference (left, input), peak (input)), -280.0);
    expectBelow ("null vs input, channel 1",
                 relativeDb (maxDifference (right, input), peak (input)), -280.0);
}

//==============================================================================
void testHissLevel()
{
    std::printf ("\nHiss level and spectral shape\n");

    const int count = 1 << 19;

    for (double asked : { -50.0, -65.0, -80.0 })
    {
        for (double speed : { 15.0, 30.0 })
        {
            TapeChannel tape;
            tape.prepare (sampleRate, 512, 1);
            tape.setSeed (0xBEEF);

            auto s = allOff();
            s.hissEnabled = true;
            s.hissLevelDb = asked;
            s.tapeSpeedIps = speed;
            tape.setSettings (s);

            const auto out = runMono (tape, std::vector<double> ((size_t) count, 0.0), 512);
            const double measured = 20.0 * std::log10 (rms (out) / referenceAmplitude);

            expectWithin ("hiss RMS re reference, " + std::to_string ((int) asked)
                              + " dB at " + std::to_string ((int) speed) + " ips",
                          measured, asked, 1.0);
        }
    }

    // Shape: the papers stress that hiss is not white, so check the replay lift at
    // both ends and the bandwidth limit at the top.
    TapeChannel tape;
    tape.prepare (sampleRate, 512, 1);
    tape.setSeed (7);

    auto s = allOff();
    s.hissEnabled = true;
    s.hissLevelDb = -60.0;
    tape.setSettings (s);

    const auto out = runMono (tape, std::vector<double> ((size_t) count, 0.0), 512);
    const double atOneK = bandPower (out, 1000.0, 4096);

    for (double hz : { 30.0, 100.0, 1000.0, 8000.0, 15000.0, 21000.0 })
        std::printf ("        hiss at %6.0f Hz  %+7.2f dB re 1 kHz\n",
                     hz, 10.0 * std::log10 (bandPower (out, hz, 4096) / atOneK));

    expectTrue ("hiss rises below 100 Hz",
                bandPower (out, 30.0, 4096) > 4.0 * atOneK, "(> 6 dB above 1 kHz)");
    expectTrue ("hiss rises above 4.5 kHz",
                bandPower (out, 8000.0, 4096) > 1.5 * atOneK, "(> 1.8 dB above 1 kHz)");
    expectTrue ("hiss rolls off past the chain bandwidth",
                bandPower (out, 21000.0, 4096) < 0.1 * bandPower (out, 8000.0, 4096),
                "(> 10 dB below 8 kHz)");
}

//==============================================================================
void testSaturationTransparency()
{
    std::printf ("\nSaturation leaves small signals alone\n");

    const int count = 1 << 15;

    for (double levelDb : { -40.0, -80.0 })
    {
        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);

        auto s = allOff();
        s.saturationEnabled = true;
        tape.setSettings (s);

        const double amplitude = referenceAmplitude * std::pow (10.0, levelDb / 20.0);
        const auto input = makeSine (count, 683, amplitude);
        const auto out = runMono (tape, input, 512);

        // Ignore the first eighth: the equalisation pair is exactly inverse, but
        // both halves still have to settle from zero state.
        double worst = 0.0;

        for (size_t n = count / 8; n < (size_t) count; ++n)
            worst = std::max (worst, std::abs (out[n] - input[n]));

        expectBelow ("null at " + std::to_string ((int) levelDb) + " dB re reference",
                     relativeDb (worst, amplitude), levelDb <= -80.0 ? -180.0 : -100.0);
    }
}

//==============================================================================
void testSaturationFrequencyDependence()
{
    std::printf ("\nGradual saturation: compression of the fundamental, dB\n");
    std::printf ("        input level        40 Hz     1 kHz    15 kHz\n");

    const int window = 1 << 15;          // analysis window
    const int count = 2 * window;        // process twice that, analyse the second half

    struct Tone { const char* name; int cycles; };
    const Tone tones[] = { { "40 Hz", 27 }, { "1 kHz", 683 }, { "15 kHz", 10240 } };

    double compressionAtSixDb[3] = { 0.0, 0.0, 0.0 };

    for (double levelDb : { 0.0, 6.0, 12.0 })
    {
        double compression[3] = { 0.0, 0.0, 0.0 };

        for (int t = 0; t < 3; ++t)
        {
            TapeChannel tape;
            tape.prepare (sampleRate, 512, 1);

            auto s = allOff();
            s.saturationEnabled = true;
            tape.setSettings (s);

            const double amplitude = referenceAmplitude * std::pow (10.0, levelDb / 20.0);
            std::vector<double> input ((size_t) count, 0.0);

            for (int n = 0; n < count; ++n)
                input[(size_t) n] = amplitude * std::sin (twoPi * (double) tones[t].cycles
                                                               * (double) n / (double) window);

            const auto out = runMono (tape, input, 512);

            const double in = binAmplitude (input.data() + window, window, tones[t].cycles);
            const double got = binAmplitude (out.data() + window, window, tones[t].cycles);

            compression[t] = 20.0 * std::log10 (got / in);
        }

        std::printf ("        %+5.0f dB re ref  %8.2f  %8.2f  %8.2f\n",
                     levelDb, compression[0], compression[1], compression[2]);

        if (levelDb == 6.0)
            for (int t = 0; t < 3; ++t)
                compressionAtSixDb[t] = compression[t];
    }

    expectTrue ("40 Hz compresses at least 2 dB more than 1 kHz",
                compressionAtSixDb[0] <= compressionAtSixDb[1] - 2.0,
                "(" + std::to_string (compressionAtSixDb[1] - compressionAtSixDb[0]) + " dB more)");
    expectTrue ("15 kHz compresses at least 2 dB more than 1 kHz",
                compressionAtSixDb[2] <= compressionAtSixDb[1] - 2.0,
                "(" + std::to_string (compressionAtSixDb[1] - compressionAtSixDb[2]) + " dB more)");
}

//==============================================================================
void testSaturationMonotonic()
{
    std::printf ("\nSaturation is monotonic and cannot fold\n");

    // The nonlinearity itself, over a range no signal will ever reach. Its slope
    // falls as exp(-2x), so past about |x| = 15 one grid step of the sweep no
    // longer moves the result by a whole double-precision ULP and the curve reads
    // as flat. Strict increase is therefore asked for only where it is resolvable;
    // over the whole range the requirement is that it never turns back down, which
    // is what "no folding" actually means.
    bool increasing = true, nonDecreasing = true, bounded = true;
    double previous = tapeSaturationCurve (-60.0);

    for (int i = 1; i <= 240000; ++i)
    {
        const double x = -60.0 + 120.0 * (double) i / 240000.0;
        const double y = tapeSaturationCurve (x);

        if (std::abs (x) <= 12.0)
            increasing = increasing && (y > previous);

        nonDecreasing = nonDecreasing && (y >= previous);
        bounded = bounded && (std::abs (y) <= 1.0);
        previous = y;
    }

    expectTrue ("curve strictly increasing over -12 to +12", increasing);
    expectTrue ("curve never turns back down, -60 to +60", nonDecreasing);
    expectTrue ("curve bounded by unity, so no folding", bounded);

    // And the whole channel: the fundamental must never turn back down.
    const int window = 1 << 14;
    const int count = 2 * window;
    double lastAmplitude = 0.0;
    bool channelMonotonic = true;
    double worstStep = 1.0e30;

    for (int step = 0; step <= 60; ++step)
    {
        const double levelDb = -80.0 + 2.0 * (double) step;

        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);

        auto s = allOff();
        s.saturationEnabled = true;
        tape.setSettings (s);

        const double amplitude = referenceAmplitude * std::pow (10.0, levelDb / 20.0);
        std::vector<double> input ((size_t) count, 0.0);

        for (int n = 0; n < count; ++n)
            input[(size_t) n] = amplitude * std::sin (twoPi * 341.0 * (double) n
                                                           / (double) window);

        const auto out = runMono (tape, input, 512);
        const double got = binAmplitude (out.data() + window, window, 341);

        if (step > 0)
        {
            channelMonotonic = channelMonotonic && (got > lastAmplitude);
            worstStep = std::min (worstStep, 20.0 * std::log10 (got / lastAmplitude));
        }

        lastAmplitude = got;
    }

    expectTrue ("channel output rises for every 2 dB input step, -80 to +40 dB",
                channelMonotonic,
                "(smallest step " + std::to_string (worstStep) + " dB)");
}

//==============================================================================
void testHighFrequencyLoss()
{
    std::printf ("\nHead and gap loss: corner lands where asked\n");

    const int window = 1 << 15;
    const int count = 2 * window;

    struct Case { double askedHz, speedIps, effectiveHz; int cornerCycles; };
    const Case cases[] =
    {
        { 18000.0, 15.0, 18000.0, 12288 },   // 18 kHz at 48 kHz is bin 12288 exactly
        {  9000.0, 30.0, 18000.0, 12288 },   // 30 ips doubles it
        {  9000.0, 15.0,  9000.0,  6144 }
    };

    for (const auto& c : cases)
    {
        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);

        auto s = allOff();
        s.hfLossEnabled = true;
        s.hfLossHz = c.askedHz;
        s.tapeSpeedIps = c.speedIps;
        tape.setSettings (s);

        expectWithin ("effective corner, " + std::to_string ((int) c.askedHz) + " Hz at "
                          + std::to_string ((int) c.speedIps) + " ips",
                      tape.effectiveHfLossHz(), c.effectiveHz, 1.0);

        const double amplitude = referenceAmplitude;
        double gainDb[3] = { 0.0, 0.0, 0.0 };
        const int cycleSet[3] = { c.cornerCycles / 4, c.cornerCycles, c.cornerCycles * 2 };

        for (int i = 0; i < 3; ++i)
        {
            if (cycleSet[i] >= window / 2)
            {
                gainDb[i] = 0.0;
                continue;
            }

            TapeChannel fresh;
            fresh.prepare (sampleRate, 512, 1);
            fresh.setSettings (s);

            std::vector<double> input ((size_t) count, 0.0);

            for (int n = 0; n < count; ++n)
                input[(size_t) n] = amplitude * std::sin (twoPi * (double) cycleSet[i]
                                                               * (double) n / (double) window);

            const auto out = runMono (fresh, input, 512);

            gainDb[i] = 20.0 * std::log10 (binAmplitude (out.data() + window, window, cycleSet[i])
                                               / binAmplitude (input.data() + window, window,
                                                               cycleSet[i]));
        }

        std::printf ("        f/4 %+6.2f dB   corner %+6.2f dB   2f %+6.2f dB\n",
                     gainDb[0], gainDb[1], gainDb[2]);

        expectWithin ("gain at the corner is -3 dB", gainDb[1], -3.01, 0.15);
        expectWithin ("two octaves below the corner is flat", gainDb[0], 0.0, 0.10);
    }
}

//==============================================================================
void testModulationNoise()
{
    std::printf ("\nModulation noise: absent in silence, scaled by signal\n");

    const int count = 1 << 19;

    // Silence in, silence out: this defect exists only when a signal is present.
    {
        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);
        tape.setSeed (5);

        auto s = allOff();
        s.modulationNoiseEnabled = true;
        s.modulationNoiseDb = -40.0;
        tape.setSettings (s);

        const auto out = runMono (tape, std::vector<double> ((size_t) count, 0.0), 512);

        expectTrue ("output is exactly silent with no input", peak (out) == 0.0,
                    "(peak " + std::to_string (peak (out)) + ")");
    }

    // With signal: the added noise must sit the asked number of dB below it.
    for (double asked : { -55.0, -40.0 })
    {
        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);
        tape.setSeed (5);

        auto s = allOff();
        s.modulationNoiseEnabled = true;
        s.modulationNoiseDb = asked;
        tape.setSettings (s);

        const auto input = makeSine (count, 10923, referenceAmplitude);   // ~1 kHz
        const auto out = runMono (tape, input, 512);

        std::vector<double> residual ((size_t) count, 0.0);

        for (int n = 0; n < count; ++n)
            residual[(size_t) n] = out[(size_t) n] - input[(size_t) n];

        const double measured = 20.0 * std::log10 (rms (residual, (size_t) count / 8)
                                                       / rms (input, (size_t) count / 8));

        expectWithin ("modulation noise re signal, " + std::to_string ((int) asked) + " dB",
                      measured, asked, 0.5);

        if (asked == -40.0)
        {
            // The skirt should hug the signal rather than fill the band.
            const double atCarrier = bandPower (residual, 1200.0, 4096);

            for (double hz : { 200.0, 700.0, 1200.0, 2000.0, 4000.0, 10000.0 })
                std::printf ("        sideband at %6.0f Hz  %+7.2f dB re 1.2 kHz\n",
                             hz, 10.0 * std::log10 (bandPower (residual, hz, 4096) / atCarrier));

            expectTrue ("sidebands are local to the signal",
                        bandPower (residual, 10000.0, 4096) < 0.01 * atCarrier,
                        "(> 20 dB down at 10 kHz)");
        }
    }
}

//==============================================================================
void testDeterminismAndChannels()
{
    std::printf ("\nDeterminism, reset and channel independence\n");

    const int count = 1 << 14;
    const auto input = makeNoise (count, 4242);

    TapeChannelSettings everything;   // all four defects on
    everything.hissLevelDb = -60.0;
    everything.modulationNoiseDb = -45.0;

    // For the size of the difference between two streams, hiss alone: the other
    // three defects are deterministic and cancel, but modulation noise would
    // otherwise dominate the residual and confuse the arithmetic.
    TapeChannelSettings hissOnly = everything;
    hissOnly.saturationEnabled = false;
    hissOnly.hfLossEnabled = false;
    hissOnly.modulationNoiseEnabled = false;

    auto runStereo = [&] (uint64_t seed, bool doubleReset, const TapeChannelSettings& s)
    {
        auto tape = std::make_unique<TapeChannel>();
        tape->prepare (sampleRate, 512, 2);
        tape->setSeed (seed);
        tape->setSettings (s);

        std::vector<double> left = input, right = input;

        if (doubleReset)
        {
            std::vector<double> throwLeft = input, throwRight = input;
            double* pointers[2] = { throwLeft.data(), throwRight.data() };
            tape->process (pointers, 2, count);
            tape->reset();
        }

        double* pointers[2] = { left.data(), right.data() };
        tape->process (pointers, 2, count);

        return std::pair { left, right };
    };

    const auto first = runStereo (0xA5A5, false, everything);
    const auto again = runStereo (0xA5A5, false, everything);
    const auto afterReset = runStereo (0xA5A5, true, everything);

    expectTrue ("same seed gives bit-identical output",
                maxDifference (first.first, again.first) == 0.0
                    && maxDifference (first.second, again.second) == 0.0);
    expectTrue ("reset restarts the noise streams exactly",
                maxDifference (first.first, afterReset.first) == 0.0
                    && maxDifference (first.second, afterReset.second) == 0.0);
    expectTrue ("the two channels are not the same signal",
                maxDifference (first.first, first.second) > 0.0);

    const auto hissLeft = runStereo (0xA5A5, false, hissOnly);
    const auto hissOther = runStereo (0x1234, false, hissOnly);
    const double hissRms = referenceAmplitude * std::pow (10.0, hissOnly.hissLevelDb / 20.0);

    std::vector<double> interChannel ((size_t) count, 0.0), interSeed ((size_t) count, 0.0);

    for (int n = 0; n < count; ++n)
    {
        interChannel[(size_t) n] = hissLeft.first[(size_t) n] - hissLeft.second[(size_t) n];
        interSeed[(size_t) n] = hissLeft.first[(size_t) n] - hissOther.first[(size_t) n];
    }

    // Two independent noise streams differ by sqrt(2) times their own RMS.
    expectWithin ("channel 0 vs channel 1 difference, dB re hiss RMS",
                  20.0 * std::log10 (rms (interChannel) / hissRms), 3.01, 1.5);
    expectWithin ("seed A vs seed B difference, dB re hiss RMS",
                  20.0 * std::log10 (rms (interSeed) / hissRms), 3.01, 1.5);
}

//==============================================================================
void testBlockSizeIndependence()
{
    std::printf ("\nBlock size independence\n");

    const int count = 4096;
    const auto input = makeNoise (count, 31337);

    TapeChannelSettings s;
    s.hissLevelDb = -60.0;

    auto run = [&] (int blockSize)
    {
        TapeChannel tape;
        tape.prepare (sampleRate, 512, 1);
        tape.setSeed (0xC0FFEE);
        tape.setSettings (s);

        return runMono (tape, input, blockSize);
    };

    const auto reference = run (count);

    for (int blockSize : { 1, 7, 512, 1000 })
        expectBelow ("blocks of " + std::to_string (blockSize) + " match one big block",
                     relativeDb (maxDifference (reference, run (blockSize)), peak (reference)),
                     -300.0);
}

//==============================================================================
void testNoAllocation()
{
    std::printf ("\nReal-time safety\n");

    TapeChannel tape;
    tape.prepare (sampleRate, 512, 2);
    tape.setSeed (11);

    TapeChannelSettings s;   // everything on
    tape.setSettings (s);

    std::vector<double> left ((size_t) 512, 0.0), right ((size_t) 512, 0.0);
    double* pointers[2] = { left.data(), right.data() };

    tape.process (pointers, 2, 512);   // warm up any lazy one-time initialisation

    allocationCount = 0;
    countingAllocations = true;

    for (int block = 0; block < 32; ++block)
    {
        for (int n = 0; n < 512; ++n)
            left[(size_t) n] = right[(size_t) n] = 0.2 * std::sin (0.05 * (double) n);

        tape.process (pointers, 2, 512);
    }

    const long duringProcess = allocationCount;

    s.hissLevelDb = -70.0;
    s.tapeSpeedIps = 30.0;
    tape.setSettings (s);

    const long duringSetSettings = allocationCount - duringProcess;
    countingAllocations = false;

    expectTrue ("process() allocates nothing", duringProcess == 0,
                "(" + std::to_string (duringProcess) + " allocations)");
    expectTrue ("setSettings() allocates nothing once prepared", duringSetSettings == 0,
                "(" + std::to_string (duringSetSettings) + " allocations)");
}

//==============================================================================
void testSampleRates()
{
    std::printf ("\nEvery rate the plugin can be asked for\n");

    for (double rate : { 11025.0, 22050.0, 44100.0, 48000.0, 96000.0, 192000.0 })
    {
        const int count = 1 << 17;

        TapeChannel tape;
        tape.prepare (rate, 512, 1);
        tape.setSeed (77);

        auto s = allOff();
        s.hissEnabled = true;
        s.hissLevelDb = -65.0;
        tape.setSettings (s);

        const auto hiss = runMono (tape, std::vector<double> ((size_t) count, 0.0), 512);
        const double measured = 20.0 * std::log10 (rms (hiss) / referenceAmplitude);

        // And once more with all four defects on a signal, to be sure nothing goes
        // non-finite where a corner has had to be clamped below Nyquist.
        TapeChannel full;
        full.prepare (rate, 512, 1);
        full.setSeed (77);
        full.setSettings (TapeChannelSettings {});

        std::vector<double> input ((size_t) count, 0.0);

        for (int n = 0; n < count; ++n)
            input[(size_t) n] = referenceAmplitude * std::sin (twoPi * 1000.0 * (double) n / rate);

        const auto out = runMono (full, input, 512);
        bool finite = true;

        for (auto x : out)
            finite = finite && std::isfinite (x);

        expectWithin ("hiss level at " + std::to_string ((int) rate) + " Hz",
                      measured, -65.0, 1.0);
        expectTrue ("all four defects stay finite at " + std::to_string ((int) rate) + " Hz",
                    finite && peak (out) < 1.0,
                    "(peak " + std::to_string (peak (out)) + ")");
    }
}

//==============================================================================
void testCombined()
{
    std::printf ("\nEverything at once, for a sanity check on levels\n");

    const int count = 1 << 18;

    TapeChannel tape;
    tape.prepare (sampleRate, 512, 1);
    tape.setSeed (2024);

    TapeChannelSettings s;
    tape.setSettings (s);

    const auto input = makeSine (count, 5461, referenceAmplitude);   // ~1 kHz at reference
    const auto out = runMono (tape, input, 512);

    std::printf ("        input RMS  %+7.2f dB re ref\n",
                 20.0 * std::log10 (rms (input) / referenceAmplitude));
    std::printf ("        output RMS %+7.2f dB re ref\n",
                 20.0 * std::log10 (rms (out) / referenceAmplitude));

    expectWithin ("output stays within a dB of the input level",
                  20.0 * std::log10 (rms (out) / rms (input)), 0.0, 1.0);

    // Noise floor with the same settings but no signal, for comparison.
    TapeChannel quiet;
    quiet.prepare (sampleRate, 512, 1);
    quiet.setSeed (2024);
    quiet.setSettings (s);

    const auto silence = runMono (quiet, std::vector<double> ((size_t) count, 0.0), 512);

    std::printf ("        no-signal noise floor %+7.2f dB re ref\n",
                 20.0 * std::log10 (rms (silence) / referenceAmplitude));
}

} // namespace

//==============================================================================
int main()
{
    std::printf ("TapeChannel tests, %.0f Hz\n", sampleRate);

    testTransparency();
    testHissLevel();
    testSaturationTransparency();
    testSaturationFrequencyDependence();
    testSaturationMonotonic();
    testHighFrequencyLoss();
    testModulationNoise();
    testDeterminismAndChannels();
    testBlockSizeIndependence();
    testNoAllocation();
    testSampleRates();
    testCombined();

    std::printf ("\n%d checks, %d failures\n", checks, failures);

    return failures == 0 ? 0 : 1;
}
