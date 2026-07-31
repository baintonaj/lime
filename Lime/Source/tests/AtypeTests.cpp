/*
  ==============================================================================

    Standalone checks for AtypeEngine.

    Builds without JUCE and without the rest of the DSP core:

        clang++ -std=c++20 -O2 -Wall -Wextra AtypeTests.cpp ../dsp/AtypeEngine.cpp \
            -o /tmp/atypetests && /tmp/atypetests

    Every check prints the number it measured, not just a verdict, because the
    acceptance criteria are the paper's published figures and the point is to be
    able to read off how close the engine gets to them.

  ==============================================================================
*/

#include "../dsp/AtypeEngine.h"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <random>
#include <vector>

//==============================================================================
// Allocation counting. Replacing the global operators is the only way to prove
// that process() does not allocate; the counter is armed only around the calls
// under test, so start-up and printing are not counted.

namespace
{
    std::atomic<long> allocationCount { 0 };
    std::atomic<bool> countingAllocations { false };
}

void* operator new (std::size_t size)
{
    if (countingAllocations.load())
        ++allocationCount;

    void* p = std::malloc (size > 0 ? size : 1);

    if (p == nullptr)
        throw std::bad_alloc();

    return p;
}

void* operator new[] (std::size_t size)
{
    return operator new (size);
}

void operator delete (void* p) noexcept                    { std::free (p); }
void operator delete[] (void* p) noexcept                  { std::free (p); }
void operator delete (void* p, std::size_t) noexcept       { std::free (p); }
void operator delete[] (void* p, std::size_t) noexcept     { std::free (p); }

//==============================================================================
namespace
{

int failures = 0;
int checksRun = 0;

void check (bool condition, const char* what)
{
    ++checksRun;

    if (! condition)
    {
        ++failures;
        std::printf ("  FAIL  %s\n", what);
    }
    else
    {
        std::printf ("  pass  %s\n", what);
    }
}

constexpr double twoPi = 6.28318530717958647692;

double toDb (double amplitude)
{
    return 20.0 * std::log10 (std::max (1.0e-300, amplitude));
}

double rms (const std::vector<double>& v, size_t from, size_t to)
{
    double sum = 0.0;

    for (size_t i = from; i < to; ++i)
        sum += v[i] * v[i];

    return (to > from) ? std::sqrt (sum / (double) (to - from)) : 0.0;
}

double peakAbs (const std::vector<double>& v)
{
    double p = 0.0;

    for (double x : v)
        p = std::max (p, std::fabs (x));

    return p;
}

std::vector<double> sine (double hz, double amplitude, double sampleRate, int numSamples)
{
    std::vector<double> out ((size_t) numSamples);

    for (int n = 0; n < numSamples; ++n)
        out[(size_t) n] = amplitude * std::sin (twoPi * hz * (double) n / sampleRate);

    return out;
}

/** Pink-ish noise: white through a one-pole cascade, normalised to a target peak.
    Deterministic seed so every run measures the same signal. */
std::vector<double> pinkNoise (double peak, int numSamples, unsigned long long seed)
{
    std::mt19937_64 rng (seed);
    std::uniform_real_distribution<double> dist (-1.0, 1.0);

    std::vector<double> out ((size_t) numSamples);
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;

    for (int n = 0; n < numSamples; ++n)
    {
        const double w = dist (rng);

        s1 = 0.99765 * s1 + w * 0.0990460;
        s2 = 0.96300 * s2 + w * 0.2965164;
        s3 = 0.57000 * s3 + w * 1.0526913;

        out[(size_t) n] = s1 + s2 + s3 + w * 0.1848;
    }

    const double p = peakAbs (out);

    if (p > 0.0)
        for (auto& x : out)
            x *= peak / p;

    return out;
}

/** Transient material: steps, clicks and gated full-scale noise bursts, which is
    the hardest case for a companding decoder to track. */
std::vector<double> transients (int numSamples, unsigned long long seed)
{
    std::mt19937_64 rng (seed);
    std::uniform_real_distribution<double> dist (-1.0, 1.0);

    std::vector<double> out ((size_t) numSamples, 0.0);

    for (int n = 0; n < numSamples; ++n)
    {
        const int phase = n % 4800;
        double x = 0.0;

        if (phase == 0)
            x = 1.0;                                   // isolated full-scale click
        else if (phase < 240)
            x = 0.9 * dist (rng);                      // full-scale noise burst
        else if (phase < 1200)
            x = 1.0;                                   // sustained step
        else if (phase < 2400)
            x = 0.001 * dist (rng);                    // deep low-level passage
        else if (phase < 3600)
            x = 0.5 * std::sin (twoPi * 100.0 * (double) n / 48000.0);

        out[(size_t) n] = x;
    }

    return out;
}

void runInBlocks (lime::AtypeEngine& engine, std::vector<double>& data, int blockSize)
{
    int n = 0;

    while (n < (int) data.size())
    {
        const int count = std::min (blockSize, (int) data.size() - n);
        double* channel = data.data() + n;
        double* const pointers[1] = { channel };

        engine.process (pointers, 1, count);
        n += count;
    }
}

/** Encode gain of a steady sine, measured from the settled tail. */
double measureEncodeGainDb (lime::AtypeEngine& engine, double hz, double amplitude,
                            double sampleRate, double seconds)
{
    const int numSamples = (int) (seconds * sampleRate);

    engine.reset();
    engine.setMode (lime::AtypeMode::encode);

    std::vector<double> data = sine (hz, amplitude, sampleRate, numSamples);
    const double inputRms = rms (data, (size_t) (numSamples / 2), (size_t) numSamples);

    runInBlocks (engine, data, 512);

    const double outputRms = rms (data, (size_t) (numSamples / 2), (size_t) numSamples);

    return toDb (outputRms) - toDb (inputRms);
}

/** Residual of encode followed by decode, in dB relative to the input. */
double measureNullDb (const std::vector<double>& input, double sampleRate)
{
    lime::AtypeEngine encoder, decoder;

    encoder.prepare (sampleRate, 512, 1);
    decoder.prepare (sampleRate, 512, 1);
    encoder.setMode (lime::AtypeMode::encode);
    decoder.setMode (lime::AtypeMode::decode);

    std::vector<double> work = input;

    runInBlocks (encoder, work, 512);
    runInBlocks (decoder, work, 512);

    std::vector<double> residual (input.size());

    for (size_t i = 0; i < input.size(); ++i)
        residual[i] = work[i] - input[i];

    return toDb (rms (residual, 0, residual.size())) - toDb (rms (input, 0, input.size()));
}

//==============================================================================
void testLowLevelGainCurve()
{
    std::printf ("\nLow-level gain curve, -60 dB sine, encode\n");
    std::printf ("  (the headline criterion: 10 dB flat to 5 kHz, rising to 15 dB at 15 kHz)\n\n");

    constexpr double sampleRate = 48000.0;
    const double probeHz[] = { 50.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0, 8000.0, 10000.0, 15000.0 };

    lime::AtypeEngine engine;
    engine.prepare (sampleRate, 512, 1);

    std::printf ("      Hz    measured   designed     target      error\n");

    double worstFlatError = 0.0;
    double gainAt15k = 0.0;

    for (double hz : probeHz)
    {
        const double measured = measureEncodeGainDb (engine, hz, 1.0e-3, sampleRate, 2.0);
        const double designed = engine.lowLevelGainDb (hz);
        const double target = lime::atypeBands::targetGainDb (hz);

        std::printf ("  %7.0f   %8.3f   %8.3f   %8.3f   %8.3f\n",
                     hz, measured, designed, target, measured - target);

        if (hz <= lime::atypeBands::flatUpToHz)
            worstFlatError = std::max (worstFlatError, std::fabs (measured - target));

        if (hz >= 15000.0)
            gainAt15k = measured;

        check (std::fabs (measured - designed) < 0.02, "measured gain matches the designed response");
    }

    std::printf ("\n  worst deviation from 10 dB at or below 5 kHz: %.3f dB\n", worstFlatError);
    std::printf ("  gain at 15 kHz: %.3f dB (target 15)\n", gainAt15k);

    check (worstFlatError < 0.6, "10 dB region uniform to better than 0.6 dB");
    check (gainAt15k > 14.0 && gainAt15k < 15.6, "15 kHz gain within 1 dB of 15 dB");

    std::printf ("\n  Fitted band gains (side chain, linear):");

    for (int b = 0; b < lime::AtypeEngine::numBands; ++b)
        std::printf (" %.6f", engine.getBandGain (b));

    std::printf ("\n  Threshold %.6f (%.1f dBFS), clipper %.6f (%.2f dBFS, %.2f dB above threshold)\n",
                 engine.getThresholdLevel(), toDb (engine.getThresholdLevel()),
                 engine.getClipLevel(), toDb (engine.getClipLevel()),
                 toDb (engine.getClipLevel()) - toDb (engine.getThresholdLevel()));
}

void testLowLevelCurveAcrossRates()
{
    std::printf ("\nLow-level curve at other sample rates (designed response)\n\n");

    for (double sampleRate : { 44100.0, 88200.0, 96000.0 })
    {
        lime::AtypeEngine engine;
        engine.prepare (sampleRate, 512, 1);

        double worst = 0.0;

        for (int i = 0; i < 200; ++i)
        {
            const double hz = 20.0 * std::pow (5000.0 / 20.0, (double) i / 199.0);
            worst = std::max (worst, std::fabs (engine.lowLevelGainDb (hz) - 10.0));
        }

        const double measured15k = measureEncodeGainDb (engine, 15000.0, 1.0e-3, sampleRate, 2.0);

        std::printf ("  %7.0f Hz:  worst error <= 5 kHz %.3f dB, measured 15 kHz gain %.3f dB\n",
                     sampleRate, worst, measured15k);

        check (worst < 0.7, "flat region holds at this sample rate");
        check (measured15k > 14.0 && measured15k < 15.6, "15 kHz gain holds at this sample rate");
    }

    // At 22.05 kHz the top of the published curve does not exist: 15 kHz is above
    // Nyquist, so the fit has nothing to aim at up there. The flat region is what
    // must survive.
    {
        lime::AtypeEngine engine;
        engine.prepare (22050.0, 512, 1);

        double worst = 0.0;

        for (int i = 0; i < 200; ++i)
        {
            const double hz = 20.0 * std::pow (5000.0 / 20.0, (double) i / 199.0);
            worst = std::max (worst, std::fabs (engine.lowLevelGainDb (hz) - 10.0));
        }

        std::printf ("    22050 Hz:  worst error <= 5 kHz %.3f dB (15 kHz is above Nyquist)\n", worst);
        check (worst < 1.0, "flat region degrades gracefully where the curve is truncated");
    }
}

void testEncodeDecodeNull()
{
    std::printf ("\nEncode then decode, residual against the input\n\n");

    constexpr double sampleRate = 48000.0;
    const int numSamples = (int) (4.0 * sampleRate);

    struct Case { const char* name; std::vector<double> signal; };

    std::vector<Case> cases;

    cases.push_back ({ "1 kHz sine at -20 dB", sine (1000.0, 0.1, sampleRate, numSamples) });
    cases.push_back ({ "1 kHz sine at   0 dB", sine (1000.0, 1.0, sampleRate, numSamples) });
    cases.push_back ({ "60 Hz sine at  -50 dB", sine (60.0, 0.00316, sampleRate, numSamples) });
    cases.push_back ({ "pink noise, peak -6 dB", pinkNoise (0.5, numSamples, 20250730u) });
    cases.push_back ({ "pink noise, peak -45 dB", pinkNoise (0.0056, numSamples, 991u) });
    cases.push_back ({ "transients, clicks and steps to full scale", transients (numSamples, 4242u) });

    for (const auto& c : cases)
    {
        const double residual = measureNullDb (c.signal, sampleRate);

        std::printf ("  %-42s  %9.2f dB\n", c.name, residual);
        check (residual < -200.0, "null better than -200 dB");
    }

    // The fit differs at other rates, so the decoder has to track a different
    // encoder there too.
    for (double rate : { 44100.0, 96000.0 })
    {
        const int count = (int) (2.0 * rate);
        const double residual = measureNullDb (transients (count, 8u), rate);

        std::printf ("  transients at %.0f Hz                          %9.2f dB\n", rate, residual);
        check (residual < -200.0, "null better than -200 dB at this sample rate");
    }
}

void testHighLevelTransparency()
{
    std::printf ("\nStatic encode characteristic, 1 kHz (side chain must vanish at high level)\n\n");

    constexpr double sampleRate = 48000.0;

    lime::AtypeEngine engine;
    engine.prepare (sampleRate, 512, 1);

    std::printf ("   input dBFS    encode gain dB\n");

    double gainAtZero = 0.0;

    for (double levelDb : { -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, -10.0, 0.0 })
    {
        const double amplitude = std::pow (10.0, levelDb / 20.0);
        const double gainDb = measureEncodeGainDb (engine, 1000.0, amplitude, sampleRate, 4.0);

        std::printf ("  %9.1f      %10.3f\n", levelDb, gainDb);

        if (levelDb == 0.0)
            gainAtZero = gainDb;
    }

    check (std::fabs (gainAtZero) < 0.5, "encode gain at peak operating level is within 0.5 dB of unity");
}

void testOvershoot()
{
    std::printf ("\nOvershoot, peak-amplitude step from silence\n\n");

    constexpr double sampleRate = 48000.0;
    const int numSamples = (int) (0.5 * sampleRate);

    lime::AtypeEngine engine;
    engine.prepare (sampleRate, 512, 1);
    engine.setMode (lime::AtypeMode::encode);

    std::vector<double> data ((size_t) numSamples, 1.0);
    runInBlocks (engine, data, 512);

    // The step amplitude is unity, so both figures are already relative to it.
    const double overshootDb = toDb (peakAbs (data));
    const double settled = toDb (rms (data, (size_t) (numSamples * 3 / 4), (size_t) numSamples));

    std::printf ("  peak output %.4f  =>  overshoot %.3f dB (criterion: about 2 dB, not more)\n",
                 peakAbs (data), overshootDb);
    std::printf ("  settled level after 0.4 s: %.3f dB relative to the step\n", settled);

    check (overshootDb <= lime::atypeBands::overshootLimitDb + 0.05, "overshoot no worse than 2 dB");
    check (overshootDb > 0.5, "overshoot is real, so the clipper is not simply muting the side chain");

    // A step at low level must overshoot far less: the compressors are inside
    // their linear region there, so nothing reaches the clipper.
    engine.reset();
    std::vector<double> quiet ((size_t) numSamples, 1.0e-3);
    runInBlocks (engine, quiet, 512);

    const double quietDb = toDb (peakAbs (quiet)) - toDb (1.0e-3);
    std::printf ("  low-level step (-60 dB) settles at %.3f dB, peak %.3f dB\n",
                 toDb (rms (quiet, (size_t) (numSamples / 2), (size_t) numSamples)) - toDb (1.0e-3), quietDb);
}

void testBandFilters()
{
    std::printf ("\nBand-defining filters\n\n");

    constexpr double sampleRate = 96000.0;   // headroom so the slopes are measurable

    lime::AtypeEngine engine;
    engine.prepare (sampleRate, 512, 1);

    const auto slopeDbPerOctave = [&] (int band, double hz)
    {
        const double a = toDb (std::abs (engine.bandResponseAt (band, hz)));
        const double b = toDb (std::abs (engine.bandResponseAt (band, 2.0 * hz)));

        return b - a;
    };

    const double band1Slope = slopeDbPerOctave (0, 800.0);     // an octave up, well into the stop band
    const double band3Slope = -slopeDbPerOctave (2, 187.5);    // an octave down
    const double band4Slope = -slopeDbPerOctave (3, 562.5);

    std::printf ("  band 1 (80 Hz LP)  slope %+7.3f dB/octave at 800 Hz\n", band1Slope);
    std::printf ("  band 3 (3 kHz HP)  slope %+7.3f dB/octave at 187 Hz\n", -band3Slope);
    std::printf ("  band 4 (9 kHz HP)  slope %+7.3f dB/octave at 562 Hz\n", -band4Slope);

    check (std::fabs (band1Slope + 12.0) < 0.3, "band 1 rolls off at 12 dB/octave");
    check (std::fabs (band3Slope + 12.0) < 0.3, "band 3 rolls off at 12 dB/octave");
    check (std::fabs (band4Slope + 12.0) < 0.3, "band 4 rolls off at 12 dB/octave");

    // Corner levels: a Butterworth corner is -3.01 dB.
    std::printf ("  corner levels: band 1 %.4f dB, band 3 %.4f dB, band 4 %.4f dB\n",
                 toDb (std::abs (engine.bandResponseAt (0, 80.0))),
                 toDb (std::abs (engine.bandResponseAt (2, 3000.0))),
                 toDb (std::abs (engine.bandResponseAt (3, 9000.0))));

    check (std::fabs (toDb (std::abs (engine.bandResponseAt (0, 80.0))) + 3.0103) < 0.01,
           "band 1 corner is at 80 Hz");
    check (std::fabs (toDb (std::abs (engine.bandResponseAt (2, 3000.0))) + 3.0103) < 0.01,
           "band 3 corner is at 3 kHz");
    check (std::fabs (toDb (std::abs (engine.bandResponseAt (3, 9000.0))) + 3.0103) < 0.01,
           "band 4 corner is at 9 kHz");

    // Bands 1 + 2 + 3 must sum to unity, magnitude and phase, since band 2 is
    // defined as the complement of the other two.
    double worstSumError = 0.0;

    for (int i = 0; i < 400; ++i)
    {
        const double hz = 10.0 * std::pow (0.45 * sampleRate / 10.0, (double) i / 399.0);
        const std::complex<double> sum = engine.bandResponseAt (0, hz)
                                        + engine.bandResponseAt (1, hz)
                                        + engine.bandResponseAt (2, hz);

        worstSumError = std::max (worstSumError, std::abs (sum - std::complex<double> { 1.0, 0.0 }));
    }

    std::printf ("  worst |band1 + band2 + band3 - 1| over 10 Hz .. 43 kHz: %.3e (%.1f dB)\n",
                 worstSumError, toDb (worstSumError));

    check (worstSumError < 1.0e-14, "bands 1 + 2 + 3 sum to unity");

    // Bands 3 and 4 should contribute comparably above 5 kHz.
    for (double hz : { 6000.0, 9000.0, 12000.0, 15000.0 })
    {
        const double c3 = engine.getBandGain (2) * std::abs (engine.bandResponseAt (2, hz));
        const double c4 = engine.getBandGain (3) * std::abs (engine.bandResponseAt (3, hz));

        std::printf ("  %6.0f Hz: band 3 contributes %.3f, band 4 %.3f (ratio %.2f)\n",
                     hz, c3, c4, c4 / c3);
    }
}

void testBlockSizeIndependence()
{
    std::printf ("\nBlock-size independence\n\n");

    constexpr double sampleRate = 48000.0;
    const int numSamples = 6000;
    const std::vector<double> input = transients (numSamples, 7u);

    for (auto mode : { lime::AtypeMode::encode, lime::AtypeMode::decode })
    {
        std::vector<double> reference = input;
        {
            lime::AtypeEngine engine;
            engine.prepare (sampleRate, 4096, 1);
            engine.setMode (mode);
            runInBlocks (engine, reference, 4096);
        }

        for (int blockSize : { 1, 7, 512 })
        {
            std::vector<double> work = input;

            lime::AtypeEngine engine;
            engine.prepare (sampleRate, blockSize, 1);
            engine.setMode (mode);
            runInBlocks (engine, work, blockSize);

            double worst = 0.0;

            for (size_t i = 0; i < work.size(); ++i)
                worst = std::max (worst, std::fabs (work[i] - reference[i]));

            std::printf ("  %s, blocks of %4d: worst difference %.3e\n",
                         mode == lime::AtypeMode::encode ? "encode" : "decode", blockSize, worst);

            check (worst == 0.0, "block size makes no difference at all");
        }
    }
}

void testChannelIndependence()
{
    std::printf ("\nChannel independence\n\n");

    constexpr double sampleRate = 48000.0;
    const int numSamples = 8000;

    const std::vector<double> left = transients (numSamples, 11u);
    const std::vector<double> right = pinkNoise (0.8, numSamples, 13u);

    std::vector<double> stereoLeft = left, stereoRight = right;
    {
        lime::AtypeEngine engine;
        engine.prepare (sampleRate, 512, 2);
        engine.setMode (lime::AtypeMode::encode);

        int n = 0;

        while (n < numSamples)
        {
            const int count = std::min (512, numSamples - n);
            double* const pointers[2] = { stereoLeft.data() + n, stereoRight.data() + n };

            engine.process (pointers, 2, count);
            n += count;
        }
    }

    std::vector<double> monoLeft = left, monoRight = right;
    {
        lime::AtypeEngine engine;
        engine.prepare (sampleRate, 512, 1);
        engine.setMode (lime::AtypeMode::encode);
        runInBlocks (engine, monoLeft, 512);
    }
    {
        lime::AtypeEngine engine;
        engine.prepare (sampleRate, 512, 1);
        engine.setMode (lime::AtypeMode::encode);
        runInBlocks (engine, monoRight, 512);
    }

    double worst = 0.0;

    for (int n = 0; n < numSamples; ++n)
    {
        worst = std::max (worst, std::fabs (stereoLeft[(size_t) n] - monoLeft[(size_t) n]));
        worst = std::max (worst, std::fabs (stereoRight[(size_t) n] - monoRight[(size_t) n]));
    }

    std::printf ("  stereo versus two mono runs: worst difference %.3e\n", worst);
    check (worst == 0.0, "channels do not interact");
}

void testNoAllocationInProcess()
{
    std::printf ("\nReal-time safety\n\n");

    constexpr double sampleRate = 48000.0;

    lime::AtypeEngine engine;
    engine.prepare (sampleRate, 512, 2);

    std::vector<double> left = transients (48000, 3u);
    std::vector<double> right = pinkNoise (0.7, 48000, 5u);

    allocationCount.store (0);
    countingAllocations.store (true);

    for (auto mode : { lime::AtypeMode::encode, lime::AtypeMode::decode })
    {
        engine.setMode (mode);
        engine.reset();

        int n = 0;

        while (n < 48000)
        {
            const int count = std::min (313, 48000 - n);
            double* const pointers[2] = { left.data() + n, right.data() + n };

            engine.process (pointers, 2, count);
            n += count;
        }
    }

    countingAllocations.store (false);
    const long allocations = allocationCount.load();

    std::printf ("  allocations during 2 x 1 s of stereo processing: %ld\n", allocations);
    check (allocations == 0, "process() and reset() do not allocate");

    std::printf ("  reported latency: %d samples\n", engine.getLatencySamples());
    check (engine.getLatencySamples() == 0, "latency is zero");
}

} // namespace

//==============================================================================
int main()
{
    std::printf ("==============================================================\n");
    std::printf (" Lime AtypeEngine checks — the four-band process, four bands\n");
    std::printf ("==============================================================\n");

    testLowLevelGainCurve();
    testLowLevelCurveAcrossRates();
    testHighLevelTransparency();
    testEncodeDecodeNull();
    testOvershoot();
    testBandFilters();
    testBlockSizeIndependence();
    testChannelIndependence();
    testNoAllocationInProcess();

    std::printf ("\n==============================================================\n");
    std::printf (" %d checks, %d failures\n", checksRun, failures);
    std::printf ("==============================================================\n");

    return failures == 0 ? 0 : 1;
}
