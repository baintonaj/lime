/*
  ==============================================================================

    Standalone tests for Lime's the calibration signals.

    Calibration.{h,cpp} depends on nothing but the standard library and
    Biquad64.h, so this builds and runs without JUCE and without the rest of the
    DSP:

        clang++ -std=c++20 -O2 -Wall -Wextra CalibrationTests.cpp \
                ../dsp/Calibration.cpp -o caltests && ./caltests

    Every figure asserted here comes from section 5 of the reference hardware manual:
    pink noise at -3 dB/octave, 20 ms nicks every 2 s at 15 dB below reference level
    level, an 850 Hz tone rising 10 percent for 30 ms every 750 ms, and a 4 s
    Auto Compare alternation.

  ==============================================================================
*/

#include "../dsp/Calibration.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

using namespace lime;

namespace
{

int failures = 0;
int checks = 0;

constexpr double twoPi = 6.283185307179586476925286766559;

void expectTrue (const std::string& what, bool ok)
{
    ++checks;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %s\n", ok ? "pass" : "FAIL", what.c_str());
}

void expectNear (const std::string& what, double measured, double expected,
                 double tolerance, const char* unit)
{
    ++checks;
    const bool ok = std::abs (measured - expected) <= tolerance;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-46s %12.5f %-4s (spec %.5f, tol %.5f)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measured, unit, expected, tolerance);
}

void expectAtMost (const std::string& what, double measured, double limit, const char* unit)
{
    ++checks;
    const bool ok = measured <= limit;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-46s %12.5f %-4s (limit %.5f)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measured, unit, limit);
}

//==============================================================================
/** In-place iterative radix-2 complex FFT. Only used to measure the spectrum of
    the generated noise, so a plain textbook implementation is enough. */
void fft (std::vector<std::complex<double>>& data)
{
    const int n = (int) data.size();

    for (int i = 1, j = 0; i < n; ++i)
    {
        int bit = n >> 1;

        for (; j & bit; bit >>= 1)
            j ^= bit;

        j ^= bit;

        if (i < j)
            std::swap (data[(size_t) i], data[(size_t) j]);
    }

    for (int len = 2; len <= n; len <<= 1)
    {
        const double theta = -twoPi / (double) len;

        for (int i = 0; i < n; i += len)
            for (int k = 0; k < len / 2; ++k)
            {
                const std::complex<double> w { std::cos (theta * (double) k), std::sin (theta * (double) k) };
                const auto u = data[(size_t) (i + k)];
                const auto v = data[(size_t) (i + k + len / 2)] * w;

                data[(size_t) (i + k)] = u + v;
                data[(size_t) (i + k + len / 2)] = u - v;
            }
    }
}

double rmsOf (const std::vector<double>& v, bool skipZeros)
{
    double sum = 0.0;
    long long count = 0;

    for (auto x : v)
    {
        if (skipZeros && x == 0.0)
            continue;

        sum += x * x;
        ++count;
    }

    return count > 0 ? std::sqrt (sum / (double) count) : 0.0;
}

double dbOf (double amplitude)
{
    return amplitude > 0.0 ? 20.0 * std::log10 (amplitude) : -1000.0;
}

//==============================================================================
void testPinkFilterDesign()
{
    std::printf ("\nPinking filter: analytic response against -3 dB/octave\n");

    for (double rate : { 44100.0, 48000.0, 96000.0 })
    {
        PinkFilter64 filter;
        filter.prepare (rate);

        // The filter's own report, measured on a 4001-point grid at design time.
        expectAtMost (std::to_string ((int) rate) + " Hz, " + std::to_string (filter.getNumSections())
                          + " sections, deviation p-p",
                      filter.getWorstDeviationDb(), 0.05, "dB");

        // Independent check: slope between octave-spaced points.
        double worstSlopeError = 0.0;

        for (double hz = 31.25; hz * 2.0 <= filter.getFitHighHz(); hz *= 2.0)
        {
            const double slope = filter.responseDbAt (hz * 2.0) - filter.responseDbAt (hz);
            worstSlopeError = std::max (worstSlopeError, std::abs (slope + 3.0102999566398120));
        }

        expectAtMost (std::to_string ((int) rate) + " Hz, worst octave-to-octave slope error",
                      worstSlopeError, 0.05, "dB");
    }

    std::printf ("\nPinking filter octave table at 48 kHz (dB relative to 1 kHz)\n");
    {
        PinkFilter64 filter;
        filter.prepare (48000.0);
        const double atOneK = filter.responseDbAt (1000.0);

        for (double hz = 31.25; hz <= 16000.0; hz *= 2.0)
            std::printf ("      %8.2f Hz  %+9.4f dB   (ideal %+9.4f)\n",
                         hz, filter.responseDbAt (hz) - atOneK,
                         -10.0 * std::log10 (hz / 1000.0));
    }

    std::printf ("\nPinking filter noise gain: quadrature against sum of squared impulse response\n");
    {
        PinkFilter64 filter;
        filter.prepare (48000.0);

        // sum h[n]^2 is exactly the mean of |H|^2 over [0, pi], which is what
        // getNoiseGain() squares to. Independent route to the same number, so it
        // validates the level calibration of every generator below.
        filter.reset();
        double energy = 0.0;

        for (int n = 0; n < (1 << 21); ++n)
            energy += std::pow (filter.processSample (n == 0 ? 1.0 : 0.0), 2.0);

        const double fromImpulse = std::sqrt (energy);

        std::printf ("      quadrature %.9f, impulse response %.9f\n",
                     filter.getNoiseGain(), fromImpulse);

        expectAtMost ("noise gain agreement",
                      std::abs (dbOf (fromImpulse) - dbOf (filter.getNoiseGain())), 0.01, "dB");
    }
}

//==============================================================================
void testPinkNoiseSpectrum()
{
    std::printf ("\nGenerated pink noise: measured octave-band spectrum\n");

    constexpr double rate = 48000.0;
    constexpr int fftSize = 8192;
    constexpr double seconds = 90.0;

    PinkNoiseGenerator pink;
    pink.prepare (rate);
    pink.setSeed (0xC0FFEEull);
    pink.reset();

    std::vector<double> window ((size_t) fftSize);
    double windowPower = 0.0;

    for (int n = 0; n < fftSize; ++n)
    {
        window[(size_t) n] = 0.5 - 0.5 * std::cos (twoPi * (double) n / (double) fftSize);
        windowPower += window[(size_t) n] * window[(size_t) n];
    }

    const int hop = fftSize / 2;
    const int totalSamples = (int) (seconds * rate);
    const int numSegments = (totalSamples - fftSize) / hop;

    std::vector<double> ring ((size_t) fftSize);
    std::vector<double> power ((size_t) (fftSize / 2 + 1), 0.0);
    std::vector<std::complex<double>> spectrum ((size_t) fftSize);

    pink.process (ring.data(), fftSize);

    for (int segment = 0; segment < numSegments; ++segment)
    {
        for (int n = 0; n < fftSize; ++n)
            spectrum[(size_t) n] = { ring[(size_t) n] * window[(size_t) n], 0.0 };

        fft (spectrum);

        for (int k = 0; k <= fftSize / 2; ++k)
            power[(size_t) k] += std::norm (spectrum[(size_t) k]);

        // Slide by one hop.
        for (int n = 0; n < fftSize - hop; ++n)
            ring[(size_t) n] = ring[(size_t) (n + hop)];

        pink.process (ring.data() + fftSize - hop, hop);
    }

    // Power spectral density, per Hz, averaged over segments.
    const double binWidth = rate / (double) fftSize;
    const double norm = 1.0 / ((double) numSegments * windowPower * binWidth);

    for (auto& p : power)
        p *= norm;

    // Mean PSD in each octave band, then a least-squares slope in dB per octave.
    struct Band { double lowHz, db; };
    std::vector<Band> bands;

    for (double lowHz = 62.5; lowHz * 2.0 <= 8000.0 * 1.001; lowHz *= 2.0)
    {
        const int firstBin = (int) std::ceil (lowHz / binWidth);
        const int lastBin = (int) std::floor (2.0 * lowHz / binWidth);

        double sum = 0.0;
        int count = 0;

        for (int k = firstBin; k <= lastBin; ++k)
        {
            sum += power[(size_t) k];
            ++count;
        }

        bands.push_back ({ lowHz, 10.0 * std::log10 (sum / (double) count) });
    }

    double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;
    const double invLog2 = 1.0 / std::log (2.0);

    for (const auto& b : bands)
    {
        const double x = std::log (b.lowHz) * invLog2;
        sumX += x;
        sumY += b.db;
        sumXX += x * x;
        sumXY += x * b.db;
    }

    const double n = (double) bands.size();
    const double slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
    const double intercept = (sumY - slope * sumX) / n;

    for (const auto& b : bands)
    {
        const double x = std::log (b.lowHz) * invLog2;
        std::printf ("      %7.1f - %7.1f Hz   PSD %+9.4f dB   residual %+7.4f dB\n",
                     b.lowHz, 2.0 * b.lowHz, b.db, b.db - (intercept + slope * x));
    }

    std::printf ("      measured slope %.4f dB/octave over %.0f s of noise\n", slope, seconds);

    expectNear ("pink noise PSD slope", slope, -3.0102999566398120, 0.1, "dB/oct");

    double worstResidual = 0.0;

    for (const auto& b : bands)
    {
        const double x = std::log (b.lowHz) * invLog2;
        worstResidual = std::max (worstResidual, std::abs (b.db - (intercept + slope * x)));
    }

    expectAtMost ("worst octave-band residual", worstResidual, 0.25, "dB");
}

//==============================================================================
void testCalibrationNoise()
{
    std::printf ("\nCalibration noise: nick timing and level\n");

    constexpr double rate = 48000.0;

    CalibrationNoiseGenerator gen;
    gen.prepare (rate);

    expectNear ("nick length", (double) gen.getNickLengthSamples() / rate * 1000.0,
                20.0, 0.02, "ms");
    expectNear ("nick period", (double) gen.getNickPeriodSamples() / rate,
                2.0, 0.001, "s");

    const int totalSamples = (int) (20.0 * rate);
    std::vector<double> signal ((size_t) totalSamples);
    gen.process (signal.data(), totalSamples);

    // Nicks are an exact mute, so runs of exact zeros locate them.
    std::vector<int> starts, lengths;
    int run = 0;

    for (int i = 0; i < totalSamples; ++i)
    {
        if (signal[(size_t) i] == 0.0)
        {
            if (run == 0)
                starts.push_back (i);

            ++run;
        }
        else if (run > 0)
        {
            lengths.push_back (run);
            run = 0;
        }
    }

    if (run > 0)
        lengths.push_back (run);

    expectTrue ("found " + std::to_string (lengths.size()) + " nicks in 20 s (expect 10)",
                lengths.size() == 10);

    int worstLength = 0;

    for (auto len : lengths)
        worstLength = std::max (worstLength, std::abs (len - 960));

    expectNear ("worst nick length", (double) (960 + worstLength) / rate * 1000.0, 20.0, 0.05, "ms");

    double worstIntervalError = 0.0;

    for (size_t i = 1; i < starts.size(); ++i)
        worstIntervalError = std::max (worstIntervalError,
                                       std::abs ((double) (starts[i] - starts[i - 1]) - 2.0 * rate));

    expectAtMost ("worst nick interval error", worstIntervalError / rate * 1000.0, 0.05, "ms");
    std::printf ("      first nick begins at %.4f s\n", (double) starts.front() / rate);

    // Level, over long enough that the estimate settles. The relative standard
    // error of an RMS estimate on a 1/f spectrum is sqrt(2 / (T * fLowest)) /
    // ln(band ratio), about 0.03 dB at two minutes, so 0.15 dB is generous.
    std::printf ("\nCalibration noise: level relative to reference level\n");

    gen.reset();
    const int longSamples = (int) (120.0 * rate);
    std::vector<double> longSignal ((size_t) longSamples);
    gen.process (longSignal.data(), longSamples);

    const double betweenNicks = rmsOf (longSignal, true);
    const double overall = rmsOf (longSignal, false);

    std::printf ("      designed RMS %.9f, measured RMS between nicks %.9f\n",
                 calibration::noiseRms, betweenNicks);

    expectNear ("measured level between nicks", dbOf (betweenNicks / calibration::referenceRms),
                -15.0, 0.15, "dB");
    expectNear ("designed level between nicks",
                dbOf (calibration::noiseRms / calibration::referenceRms),
                -15.0, 1.0e-9, "dB");
    std::printf ("      including the nicks: %+.4f dB rel. reference level, %+.4f dBFS\n",
                 dbOf (overall / calibration::referenceRms), dbOf (overall));

    // Nothing may clip: the whole point of the uniform white source.
    double peak = 0.0;

    for (auto x : longSignal)
        peak = std::max (peak, std::abs (x));

    std::printf ("      peak %.6f (%.2f dBFS), crest factor %.2f dB\n",
                 peak, dbOf (peak), dbOf (peak / betweenNicks));

    expectAtMost ("peak stays below full scale", peak, 1.0, "");
}

//==============================================================================
void testCalibrationTone()
{
    std::printf ("\nCalibration tone: frequencies, timing and envelope\n");

    constexpr double rate = 48000.0;

    CalibrationToneGenerator tone;
    tone.prepare (rate);

    expectNear ("burst length", (double) tone.getBurstSamples() / rate * 1000.0, 30.0, 0.02, "ms");
    expectNear ("burst period", (double) tone.getPeriodSamples() / rate * 1000.0, 750.0, 0.02, "ms");

    const int totalSamples = (int) (6.0 * rate);
    std::vector<double> signal ((size_t) totalSamples);
    tone.process (signal.data(), totalSamples);

    // Per-sample frequency from the three-term sinusoid recurrence
    // x[n+1] + x[n-1] = 2 cos(w) x[n], evaluated only where the divisor is
    // comfortably away from zero. Exact for a constant-frequency sine, so the
    // only samples it misreads are the one or two either side of a transition.
    std::vector<double> freq ((size_t) totalSamples, -1.0);

    for (int n = 1; n + 1 < totalSamples; ++n)
    {
        if (std::abs (signal[(size_t) n]) < 0.5)
            continue;

        const double cosw = (signal[(size_t) (n + 1)] + signal[(size_t) (n - 1)])
                            / (2.0 * signal[(size_t) n]);

        if (cosw >= -1.0 && cosw <= 1.0)
            freq[(size_t) n] = rate * std::acos (cosw) / twoPi;
    }

    // Steady and modulated frequency, sampled well inside each region.
    double baseSum = 0.0, modSum = 0.0;
    int baseCount = 0, modCount = 0;

    const int period = tone.getPeriodSamples();
    const int burst = tone.getBurstSamples();

    for (int n = 1; n + 1 < totalSamples; ++n)
    {
        if (freq[(size_t) n] < 0.0)
            continue;

        const int phase = n % period;
        const int burstFrom = period - burst;

        if (phase > 200 && phase < burstFrom - 200)
        {
            baseSum += freq[(size_t) n];
            ++baseCount;
        }
        else if (phase > burstFrom + 200 && phase < period - 200)
        {
            modSum += freq[(size_t) n];
            ++modCount;
        }
    }

    expectNear ("unmodulated frequency", baseSum / (double) baseCount, 850.0, 0.01, "Hz");
    expectNear ("modulated frequency", modSum / (double) modCount, 935.0, 0.01, "Hz");
    expectNear ("modulation depth",
                100.0 * ((modSum / (double) modCount) / (baseSum / (double) baseCount) - 1.0),
                10.0, 0.01, "%");

    // Burst edges, from a forward-filled classification. Resolution is limited
    // by the stretch of samples near a zero crossing that carry no estimate,
    // about 0.4 ms at 850 Hz.
    std::vector<int> burstStarts, burstEnds;
    bool modulated = false;

    for (int n = 1; n + 1 < totalSamples; ++n)
    {
        if (freq[(size_t) n] < 0.0)
            continue;

        const bool nowModulated = freq[(size_t) n] > 892.5;

        if (nowModulated && ! modulated)
            burstStarts.push_back (n);
        else if (! nowModulated && modulated)
            burstEnds.push_back (n);

        modulated = nowModulated;
    }

    expectTrue ("found " + std::to_string (burstStarts.size()) + " bursts in 6 s (expect 8)",
                burstStarts.size() == 8);

    double worstDuration = 0.0, worstPeriod = 0.0;

    for (size_t i = 0; i < burstEnds.size() && i < burstStarts.size(); ++i)
        worstDuration = std::max (worstDuration,
                                  std::abs ((double) (burstEnds[i] - burstStarts[i]) / rate * 1000.0 - 30.0));

    for (size_t i = 1; i < burstStarts.size(); ++i)
        worstPeriod = std::max (worstPeriod,
                                std::abs ((double) (burstStarts[i] - burstStarts[i - 1]) / rate * 1000.0 - 750.0));

    expectAtMost ("worst measured burst duration error", worstDuration, 1.5, "ms");
    expectAtMost ("worst measured burst period error", worstPeriod, 1.0, "ms");
    std::printf ("      first burst begins at %.4f s\n", (double) burstStarts.front() / rate);

    // Envelope, from a Hann-weighted 20 ms sliding RMS. A frequency change does
    // not change the energy in the window, so this reads flat unless there is
    // genuine amplitude modulation. This is the check the manual's reasoning for
    // choosing FM over AM demands.
    const int windowLength = (int) std::lround (0.020 * rate);
    std::vector<double> envelopeWindow ((size_t) windowLength);
    double windowSum = 0.0;

    for (int n = 0; n < windowLength; ++n)
    {
        envelopeWindow[(size_t) n] = 0.5 - 0.5 * std::cos (twoPi * (double) n / (double) windowLength);
        windowSum += envelopeWindow[(size_t) n];
    }

    double worstEnvelopeDb = 0.0;
    double worstAtSeconds = 0.0;

    for (int start = 0; start + windowLength < totalSamples; start += 4)
    {
        double sum = 0.0;

        for (int n = 0; n < windowLength; ++n)
            sum += envelopeWindow[(size_t) n]
                   * signal[(size_t) (start + n)] * signal[(size_t) (start + n)];

        const double envelope = std::sqrt (2.0 * sum / windowSum);
        const double deviation = std::abs (dbOf (envelope / tone.getAmplitude()));

        if (deviation > worstEnvelopeDb)
        {
            worstEnvelopeDb = deviation;
            worstAtSeconds = (double) (start + windowLength / 2) / rate;
        }
    }

    std::printf ("      worst envelope deviation %.6f dB at %.4f s\n", worstEnvelopeDb, worstAtSeconds);
    expectAtMost ("envelope deviation through the transitions", worstEnvelopeDb, 0.01, "dB");
}

//==============================================================================
void testAutoCompare()
{
    std::printf ("\nAuto Compare: 4 s alternation\n");

    constexpr double rate = 48000.0;

    AutoCompare compare;
    compare.prepare (rate);

    expectNear ("segment length", (double) compare.getSegmentSamples() / rate, 4.0, 0.001, "s");

    compare.setRunning (true);
    expectTrue ("starts on tape", compare.getCurrentSource() == AutoCompareSource::tape);

    const int totalSamples = (int) (20.0 * rate);
    std::vector<double> tape ((size_t) totalSamples), reference ((size_t) totalSamples);
    std::vector<double> out ((size_t) totalSamples);

    // Distinguishable, so the selection can be checked sample by sample.
    for (int n = 0; n < totalSamples; ++n)
    {
        tape[(size_t) n] = 1.0 + (double) n;
        reference[(size_t) n] = -1.0 - (double) n;
    }

    std::vector<int> sourceAt ((size_t) totalSamples);

    for (int n = 0; n < totalSamples; ++n)
    {
        sourceAt[(size_t) n] = compare.getCurrentSource() == AutoCompareSource::tape ? 0 : 1;
        compare.process (tape.data() + n, reference.data() + n, out.data() + n, 1);
    }

    // Transitions.
    std::vector<double> transitions;

    for (int n = 1; n < totalSamples; ++n)
        if (sourceAt[(size_t) n] != sourceAt[(size_t) (n - 1)])
            transitions.push_back ((double) n / rate);

    std::printf ("      transitions at");

    for (auto t : transitions)
        std::printf (" %.4f", t);

    std::printf (" s\n");

    expectTrue ("four transitions in 20 s", transitions.size() == 4);

    double worstTransition = 0.0;

    for (size_t i = 0; i < transitions.size(); ++i)
        worstTransition = std::max (worstTransition,
                                    std::abs (transitions[i] - 4.0 * (double) (i + 1)));

    expectAtMost ("worst transition time error", worstTransition * 1000.0, 0.05, "ms");

    bool selectionCorrect = true;
    bool sourceReportCorrect = true;

    for (int n = 0; n < totalSamples; ++n)
    {
        const int segment = (n / compare.getSegmentSamples()) % 2;
        const double expected = segment == 0 ? tape[(size_t) n] : reference[(size_t) n];

        if (out[(size_t) n] != expected)
            selectionCorrect = false;

        if (sourceAt[(size_t) n] != segment)
            sourceReportCorrect = false;
    }

    expectTrue ("output takes the selected source exactly", selectionCorrect);
    expectTrue ("getCurrentSource() agrees with the output", sourceReportCorrect);

    // Parked state.
    compare.setRunning (false);
    compare.process (tape.data(), reference.data(), out.data(), 1000);

    bool parkedOnTape = compare.getCurrentSource() == AutoCompareSource::tape;

    for (int n = 0; n < 1000; ++n)
        if (out[(size_t) n] != tape[(size_t) n])
            parkedOnTape = false;

    expectTrue ("parks on tape when not running", parkedOnTape);
}

//==============================================================================
/** The multi-channel overload: one segment clock for however many channels.

    The clock advances once per *sample*, not once per channel-sample, so a stereo pair
    switches tape to reference on exactly the same sample index and the segments keep
    their configured length whatever the channel count. Advancing per channel-sample
    would halve the segment in stereo and desynchronise the two sides.
*/
void testAutoCompareStereo()
{
    std::printf ("\nAuto Compare: stereo shares one sample clock\n");

    constexpr double rate = 48000.0;
    const int totalSamples = (int) (10.0 * rate);

    // Every source distinguishable at every sample, so the selection can be decoded
    // from the output alone.
    std::vector<double> tapeLeft ((size_t) totalSamples), tapeRight ((size_t) totalSamples);
    std::vector<double> reference ((size_t) totalSamples);

    for (int n = 0; n < totalSamples; ++n)
    {
        tapeLeft[(size_t) n] = 1.0 + (double) n;
        tapeRight[(size_t) n] = -1.0 - (double) n;
        reference[(size_t) n] = 1.0e7 + (double) n;
    }

    std::vector<double> outLeft ((size_t) totalSamples), outRight ((size_t) totalSamples);

    AutoCompare stereo;
    stereo.prepare (rate);
    stereo.setRunning (true);

    // Irregular blocks, so nothing depends on the block boundary landing anywhere.
    const int blockSizes[] = { 480, 33, 1024, 7, 512 };
    int at = 0, which = 0;

    while (at < totalSamples)
    {
        const int len = std::min (blockSizes[which % 5], totalSamples - at);

        const double* tape[2] = { tapeLeft.data() + at, tapeRight.data() + at };
        double* out[2] = { outLeft.data() + at, outRight.data() + at };

        stereo.process (tape, reference.data() + at, out, 2, len);
        at += len;
        ++which;
    }

    // Decode each channel's selection from its output values.
    const auto switchIndices = [totalSamples] (const std::vector<double>& out,
                                               const std::vector<double>& tape)
    {
        std::vector<int> switches;

        for (int n = 1; n < totalSamples; ++n)
        {
            const bool wasTape = out[(size_t) (n - 1)] == tape[(size_t) (n - 1)];
            const bool isTape = out[(size_t) n] == tape[(size_t) n];

            if (wasTape != isTape)
                switches.push_back (n);
        }

        return switches;
    };

    const auto leftSwitches = switchIndices (outLeft, tapeLeft);
    const auto rightSwitches = switchIndices (outRight, tapeRight);

    // (a) Both channels switch on exactly the same sample index, every time.
    expectTrue ("both channels switch on the same sample",
                leftSwitches == rightSwitches);

    // (b) The segment length in samples is segmentSamples regardless of channel count:
    // every switch lands on a whole multiple of the segment, so nothing about running
    // two channels shortened or stretched the clock.
    const int segment = stereo.getSegmentSamples();

    expectTrue ("two switches in 10 s of stereo", leftSwitches.size() == 2);

    bool onSegmentGrid = ! leftSwitches.empty();

    for (size_t i = 0; i < leftSwitches.size(); ++i)
        if (leftSwitches[i] != (int) (i + 1) * segment)
            onSegmentGrid = false;

    expectTrue ("stereo switches land on the segment grid", onSegmentGrid);

    // The selection itself, on both channels: tape for the first segment, reference for
    // the second, repeating — and the reference segment carries the *mono* reference on
    // every channel.
    bool selectionCorrect = true;

    for (int n = 0; n < totalSamples; ++n)
    {
        const bool onTape = (n / segment) % 2 == 0;

        if (outLeft[(size_t) n] != (onTape ? tapeLeft[(size_t) n] : reference[(size_t) n]))
            selectionCorrect = false;

        if (outRight[(size_t) n] != (onTape ? tapeRight[(size_t) n] : reference[(size_t) n]))
            selectionCorrect = false;
    }

    expectTrue ("both channels take the selected source exactly", selectionCorrect);

    // (c) The single-channel overload still behaves as before: run the same programme
    // mono, in the same blocks, and the switch indices must be identical to the stereo
    // run's — the segment is a length in samples, not in channel-samples.
    std::vector<double> outMono ((size_t) totalSamples);

    AutoCompare mono;
    mono.prepare (rate);
    mono.setRunning (true);

    at = 0;
    which = 0;

    while (at < totalSamples)
    {
        const int len = std::min (blockSizes[which % 5], totalSamples - at);
        mono.process (tapeLeft.data() + at, reference.data() + at, outMono.data() + at, len);
        at += len;
        ++which;
    }

    expectTrue ("mono and stereo switch on identical samples",
                switchIndices (outMono, tapeLeft) == leftSwitches);
    expectTrue ("mono output matches the stereo left channel", outMono == outLeft);
}

//==============================================================================
void testDetector()
{
    std::printf ("\nCalibration noise detector\n");

    constexpr double rate = 48000.0;
    const int totalSamples = (int) (14.0 * rate);

    auto runDetector = [&] (const std::vector<double>& signal, int blockSize)
    {
        CalibrationNoiseDetector detector;
        detector.prepare (rate);

        long long firedAt = -1;

        for (int n = 0; n < (int) signal.size(); n += blockSize)
        {
            const int count = std::min (blockSize, (int) signal.size() - n);
            detector.push (signal.data() + n, count);

            if (firedAt < 0 && detector.isCalibrationNoisePresent())
                firedAt = n + count;
        }

        std::printf ("      nicks %d, valid intervals %d, level %.2f dBFS, noise-likeness %+.2f dB\n",
                     detector.getNickCount(), detector.getValidIntervals(),
                     detector.getLevelDbFs(), detector.getNoiseLikenessDb());

        return std::pair<bool, double> { detector.isCalibrationNoisePresent(),
                                         firedAt < 0 ? -1.0 : (double) firedAt / rate };
    };

    // Real calibration noise.
    std::vector<double> calibrationNoise ((size_t) totalSamples);
    {
        CalibrationNoiseGenerator gen;
        gen.prepare (rate);
        gen.setSeed (0x1234ull);
        gen.reset();
        gen.process (calibrationNoise.data(), totalSamples);
    }

    std::printf ("    calibration noise:\n");
    auto calibrationResult = runDetector (calibrationNoise, 512);
    expectTrue ("fires on calibration noise", calibrationResult.first);
    std::printf ("      armed after %.3f s\n", calibrationResult.second);

    // Silence.
    std::printf ("    silence:\n");
    std::vector<double> silence ((size_t) totalSamples, 0.0);
    expectTrue ("does not fire on silence", ! runDetector (silence, 512).first);

    // Steady 850 Hz tone at reference level.
    std::printf ("    steady 850 Hz tone at reference level:\n");
    std::vector<double> sine ((size_t) totalSamples);

    for (int n = 0; n < totalSamples; ++n)
        sine[(size_t) n] = std::sin (twoPi * 850.0 * (double) n / rate);

    expectTrue ("does not fire on a steady sine", ! runDetector (sine, 512).first);

    // The full calibration tone, which does have periodic events in it.
    std::printf ("    calibration tone (850 Hz, FM burst every 750 ms):\n");
    std::vector<double> calibrationTone ((size_t) totalSamples);
    {
        CalibrationToneGenerator gen;
        gen.prepare (rate);
        gen.process (calibrationTone.data(), totalSamples);
    }

    expectTrue ("does not fire on calibration tone", ! runDetector (calibrationTone, 512).first);

    // Uninterrupted pink noise: the Auto Compare internal reference.
    std::printf ("    uninterrupted pink noise:\n");
    std::vector<double> pink ((size_t) totalSamples);
    {
        PinkNoiseGenerator gen;
        gen.prepare (rate);
        gen.setSeed (0x9999ull);
        gen.reset();
        gen.process (pink.data(), totalSamples);
    }

    expectTrue ("does not fire on uninterrupted pink noise", ! runDetector (pink, 512).first);

    // calibration noise off a recorder in poor shape, which must still be recognised:
    // the noise-likeness test is only there to exclude tones.
    std::printf ("    calibration noise through three cascaded 5 kHz one-poles:\n");
    std::vector<double> lossy = calibrationNoise;
    {
        const double coeff = 1.0 - std::exp (-twoPi * 5000.0 / rate);

        for (int pole = 0; pole < 3; ++pole)
        {
            double state = 0.0;

            for (auto& x : lossy)
            {
                state += coeff * (x - state);
                x = state;
            }
        }
    }

    expectTrue ("still fires on calibration noise with 15 dB/octave of HF loss",
                runDetector (lossy, 512).first);

    // Block-size independence.
    std::printf ("    block-size independence:\n");
    const auto atOne = runDetector (calibrationNoise, 1);
    const auto atSeven = runDetector (calibrationNoise, 7);

    expectTrue ("arming time is block-size independent",
                std::abs (atOne.second - calibrationResult.second) < 512.0 / rate + 1.0e-12
                    && std::abs (atSeven.second - calibrationResult.second) < 512.0 / rate + 1.0e-12);
}

//==============================================================================
template <typename Generator>
std::vector<double> generateInBlocks (double rate, int totalSamples, int blockSize,
                                      std::uint64_t seed, bool useSeed)
{
    Generator gen;
    gen.prepare (rate);

    if constexpr (requires { gen.setSeed (seed); })
        if (useSeed)
            gen.setSeed (seed);

    gen.reset();

    std::vector<double> out ((size_t) totalSamples);

    for (int n = 0; n < totalSamples; n += blockSize)
        gen.process (out.data() + n, std::min (blockSize, totalSamples - n));

    return out;
}

void testDeterminismAndBlockSize()
{
    std::printf ("\nDeterminism and block-size independence\n");

    constexpr double rate = 48000.0;
    const int totalSamples = (int) (5.0 * rate);

    for (int blockSize : { 1, 7, 512 })
    {
        const auto reference = generateInBlocks<CalibrationNoiseGenerator> (rate, totalSamples, totalSamples, 42ull, true);
        const auto blocked = generateInBlocks<CalibrationNoiseGenerator> (rate, totalSamples, blockSize, 42ull, true);

        expectTrue ("calibration noise, blocks of " + std::to_string (blockSize) + ", bit identical",
                    reference == blocked);
    }

    for (int blockSize : { 1, 7, 512 })
    {
        const auto reference = generateInBlocks<PinkNoiseGenerator> (rate, totalSamples, totalSamples, 42ull, true);
        const auto blocked = generateInBlocks<PinkNoiseGenerator> (rate, totalSamples, blockSize, 42ull, true);

        expectTrue ("pink noise, blocks of " + std::to_string (blockSize) + ", bit identical",
                    reference == blocked);
    }

    for (int blockSize : { 1, 7, 512 })
    {
        const auto reference = generateInBlocks<CalibrationToneGenerator> (rate, totalSamples, totalSamples, 0ull, false);
        const auto blocked = generateInBlocks<CalibrationToneGenerator> (rate, totalSamples, blockSize, 0ull, false);

        expectTrue ("calibration tone, blocks of " + std::to_string (blockSize) + ", bit identical",
                    reference == blocked);
    }

    // Same seed twice, and two different seeds.
    const auto seedA = generateInBlocks<CalibrationNoiseGenerator> (rate, 48000, 128, 7ull, true);
    const auto seedAagain = generateInBlocks<CalibrationNoiseGenerator> (rate, 48000, 333, 7ull, true);
    const auto seedB = generateInBlocks<CalibrationNoiseGenerator> (rate, 48000, 128, 8ull, true);

    expectTrue ("same seed reproduces the same noise exactly", seedA == seedAagain);
    expectTrue ("a different seed gives different noise", seedA != seedB);

    // Auto Compare selection is also block-size independent.
    for (int blockSize : { 1, 7, 512 })
    {
        const int samples = (int) (20.0 * rate);
        std::vector<double> tape ((size_t) samples), reference ((size_t) samples), out ((size_t) samples);

        for (int n = 0; n < samples; ++n)
        {
            tape[(size_t) n] = (double) n;
            reference[(size_t) n] = -(double) n;
        }

        AutoCompare compare;
        compare.prepare (rate);
        compare.setRunning (true);

        for (int n = 0; n < samples; n += blockSize)
            compare.process (tape.data() + n, reference.data() + n, out.data() + n,
                             std::min (blockSize, samples - n));

        bool ok = true;
        const int segment = compare.getSegmentSamples();

        for (int n = 0; n < samples; ++n)
        {
            const double expected = ((n / segment) % 2 == 0) ? tape[(size_t) n] : reference[(size_t) n];

            if (out[(size_t) n] != expected)
                ok = false;
        }

        expectTrue ("Auto Compare, blocks of " + std::to_string (blockSize) + ", same selection", ok);
    }

    // Every rate tier the plugin can be asked for.
    std::printf ("\nPrepares and generates at every sample rate\n");

    for (double sampleRate : { 11025.0, 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 192000.0 })
    {
        CalibrationNoiseGenerator noise;
        CalibrationToneGenerator tone;
        noise.prepare (sampleRate);
        tone.prepare (sampleRate);

        const int samples = (int) (3.0 * sampleRate);
        std::vector<double> a ((size_t) samples), b ((size_t) samples);
        noise.process (a.data(), samples);
        tone.process (b.data(), samples);

        bool finite = true;

        for (int n = 0; n < samples; ++n)
            if (! std::isfinite (a[(size_t) n]) || ! std::isfinite (b[(size_t) n]))
                finite = false;

        const double noiseLevel = dbOf (rmsOf (a, true) / calibration::referenceRms);
        const double toneLevel = dbOf (rmsOf (b, false) / calibration::referenceRms);

        std::printf ("      %7.0f Hz: %2d sections, pink fit %.5f dB p-p, noise %+.3f dB, tone %+.3f dB\n",
                     sampleRate, noise.getNoise().getFilter().getNumSections(),
                     noise.getNoise().getFilter().getWorstDeviationDb(), noiseLevel, toneLevel);

        expectTrue (std::to_string ((int) sampleRate) + " Hz output is finite", finite);
        expectNear (std::to_string ((int) sampleRate) + " Hz noise level", noiseLevel, -15.0, 0.6, "dB");
        expectNear (std::to_string ((int) sampleRate) + " Hz tone level", toneLevel, 0.0, 0.01, "dB");
    }
}

} // namespace

//==============================================================================
int main()
{
    std::printf ("Lime calibration tests\n");
    std::printf ("reference level = amplitude %.1f, RMS %.6f\n",
                 calibration::referenceAmplitude, calibration::referenceRms);

    testPinkFilterDesign();
    testPinkNoiseSpectrum();
    testCalibrationNoise();
    testCalibrationTone();
    testAutoCompare();
    testAutoCompareStereo();
    testDetector();
    testDeterminismAndBlockSize();

    std::printf ("\n%d checks, %d failure%s\n", checks, failures, failures == 1 ? "" : "s");

    return failures == 0 ? 0 : 1;
}
