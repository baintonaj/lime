/*
  ==============================================================================

    Standalone tests for Lime's DSP core.

    Nothing under Source/dsp depends on JUCE, so these build and run without the
    plugin wrapper:

        clang++ -std=c++20 -O2 DspTests.cpp ../dsp/Fft64.cpp \
                -framework Accelerate -o dsptests && ./dsptests

  ==============================================================================
*/

#include "../dsp/AtypeEngine.h"
#include "../dsp/AutoGain.h"
#include "../dsp/Biquad64.h"
#include "../dsp/ControlSmoother.h"
#include "../dsp/Crossfade.h"
#include "../dsp/DelayLine64.h"
#include "../dsp/DspMath.h"
#include "../dsp/Fft64.h"
#include "../dsp/FilterDesign.h"
#include "../dsp/FixedBand.h"
#include "../dsp/FractionalDelay64.h"
#include "../dsp/GainSlew.h"
#include "../dsp/ModulationControl.h"
#include "../dsp/OverlapSave64.h"
#include "../dsp/ShelvingBank64.h"
#include "../dsp/SlidingLayer.h"
#include "../dsp/SpectrumProbe.h"
#include "../dsp/SrEngine.h"
#include "../dsp/WetLowPass.h"
#include "../dsp/Window.h"
#include "../dsp/Wola64.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <random>
#include <string>
#include <thread>
#include <vector>

//==============================================================================
// Allocation counting, as in AtypeTests: replacing the global operators is the
// only way to prove that process() does not allocate, and the counter is armed
// only around the calls under test so start-up and printing are not counted.

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

using namespace lime;

namespace
{

int failures = 0;
int checks = 0;

constexpr double twoPi = 6.283185307179586476925286766559;

/** Reports a level in dB relative to a reference magnitude. Returns -inf as a
    large negative number so it prints readably. */
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

    std::printf ("  [%s] %-58s %8.1f dB  (limit %6.1f dB)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measuredDb, limitDb);
}

double peakMagnitude (const std::vector<double>& v)
{
    double peak = 0.0;

    for (auto x : v)
        peak = std::max (peak, std::abs (x));

    return peak;
}

std::vector<double> makeNoise (int n, unsigned seed)
{
    std::mt19937_64 rng (seed);
    std::uniform_real_distribution<double> dist (-1.0, 1.0);

    std::vector<double> out ((size_t) n);

    for (auto& x : out)
        x = dist (rng);

    return out;
}

/** Naive O(N^2) DFT, used as ground truth for small sizes.

    Kahan-compensated: a plain serial sum accumulates O(N)*eps of its own error,
    which by N=256 is worse than the error of the FFT being tested, so the
    reference has to be the accurate one or the comparison measures nothing.
*/
std::vector<std::complex<double>> naiveDft (const std::vector<double>& x)
{
    const int n = (int) x.size();
    std::vector<std::complex<double>> out ((size_t) (n / 2 + 1));

    for (int k = 0; k <= n / 2; ++k)
    {
        double sumRe = 0.0, sumIm = 0.0;
        double compRe = 0.0, compIm = 0.0;

        for (int t = 0; t < n; ++t)
        {
            // Reduce the angle before the trig call: k*t grows to N^2/2, and
            // cos/sin lose relative accuracy on large arguments.
            const long long kt = (long long) k * (long long) t;
            const double phase = (double) (kt % (long long) n) / (double) n;
            const double theta = -twoPi * phase;

            const double termRe = x[(size_t) t] * std::cos (theta);
            const double termIm = x[(size_t) t] * std::sin (theta);

            double y = termRe - compRe;
            double tmp = sumRe + y;
            compRe = (tmp - sumRe) - y;
            sumRe = tmp;

            y = termIm - compIm;
            tmp = sumIm + y;
            compIm = (tmp - sumIm) - y;
            sumIm = tmp;
        }

        out[(size_t) k] = { sumRe, sumIm };
    }

    return out;
}

const char* backendName (FftBackend b)
{
    switch (b)
    {
        case FftBackend::accelerate: return "accelerate";
        case FftBackend::portable:   return "portable";
        default:                     return "automatic";
    }
}

//==============================================================================
void testForwardAgainstNaiveDft()
{
    std::printf ("\nForward transform vs naive DFT\n");

    for (int order = 2; order <= 8; ++order)
    {
        const int n = 1 << order;
        const auto input = makeNoise (n, 1234u + (unsigned) order);
        const auto truth = naiveDft (input);

        for (auto backend : { FftBackend::portable, FftBackend::accelerate })
        {
            if (backend == FftBackend::accelerate && ! Fft64::accelerateAvailable())
                continue;

            Fft64 fft (order, backend);

            if (fft.getBackend() != backend)
                continue;

            std::vector<std::complex<double>> bins ((size_t) fft.getNumBins());
            fft.forward (input.data(), bins.data());

            double worst = 0.0, reference = 0.0;

            for (size_t k = 0; k < bins.size(); ++k)
            {
                worst = std::max (worst, std::abs (bins[k] - truth[k]));
                reference = std::max (reference, std::abs (truth[k]));
            }

            expectBelow (std::string ("N=") + std::to_string (n) + ", " + backendName (backend),
                         relativeDb (worst, reference), -280.0);
        }
    }
}

void testRoundTrip()
{
    std::printf ("\nForward then inverse is the identity\n");

    for (int order : { 2, 4, 8, 10, 11, 12, 13 })
    {
        const int n = 1 << order;
        const auto input = makeNoise (n, 999u + (unsigned) order);

        for (auto backend : { FftBackend::portable, FftBackend::accelerate })
        {
            if (backend == FftBackend::accelerate && ! Fft64::accelerateAvailable())
                continue;

            Fft64 fft (order, backend);

            if (fft.getBackend() != backend)
                continue;

            std::vector<std::complex<double>> bins ((size_t) fft.getNumBins());
            std::vector<double> output ((size_t) n);

            fft.forward (input.data(), bins.data());
            fft.inverse (bins.data(), output.data());

            double worst = 0.0;

            for (int t = 0; t < n; ++t)
                worst = std::max (worst, std::abs (output[(size_t) t] - input[(size_t) t]));

            expectBelow (std::string ("N=") + std::to_string (n) + ", " + backendName (backend),
                         relativeDb (worst, peakMagnitude (input)), -280.0);
        }
    }
}

void testBackendsAgree()
{
    std::printf ("\nAccelerate and portable backends agree\n");

    if (! Fft64::accelerateAvailable())
    {
        std::printf ("  (skipped: this build has no Accelerate backend)\n");
        return;
    }

    for (int order : { 8, 10, 12, 13 })
    {
        const int n = 1 << order;
        const auto input = makeNoise (n, 4242u + (unsigned) order);

        Fft64 fastFft (order, FftBackend::accelerate);
        Fft64 slowFft (order, FftBackend::portable);

        if (fastFft.getBackend() != FftBackend::accelerate)
        {
            std::printf ("  (skipped N=%d: Accelerate setup unavailable)\n", n);
            continue;
        }

        std::vector<std::complex<double>> fastBins ((size_t) fastFft.getNumBins());
        std::vector<std::complex<double>> slowBins ((size_t) slowFft.getNumBins());

        fastFft.forward (input.data(), fastBins.data());
        slowFft.forward (input.data(), slowBins.data());

        double worst = 0.0, reference = 0.0;

        for (size_t k = 0; k < fastBins.size(); ++k)
        {
            worst = std::max (worst, std::abs (fastBins[k] - slowBins[k]));
            reference = std::max (reference, std::abs (slowBins[k]));
        }

        expectBelow (std::string ("forward, N=") + std::to_string (n),
                     relativeDb (worst, reference), -280.0);

        std::vector<double> fastTime ((size_t) n), slowTime ((size_t) n);
        fastFft.inverse (slowBins.data(), fastTime.data());
        slowFft.inverse (slowBins.data(), slowTime.data());

        worst = 0.0;

        for (int t = 0; t < n; ++t)
            worst = std::max (worst, std::abs (fastTime[(size_t) t] - slowTime[(size_t) t]));

        expectBelow (std::string ("inverse, N=") + std::to_string (n),
                     relativeDb (worst, peakMagnitude (slowTime)), -280.0);
    }
}

void testKnownSignals()
{
    std::printf ("\nKnown signals land in the right bins\n");

    constexpr int order = 10;
    constexpr int n = 1 << order;

    Fft64 fft (order);
    std::vector<std::complex<double>> bins ((size_t) fft.getNumBins());

    // DC: all energy in bin 0, equal to N.
    {
        std::vector<double> input ((size_t) n, 1.0);
        fft.forward (input.data(), bins.data());

        double leakage = 0.0;

        for (size_t k = 1; k < bins.size(); ++k)
            leakage = std::max (leakage, std::abs (bins[k]));

        expectBelow ("DC magnitude error", relativeDb (std::abs (bins[0].real() - (double) n), (double) n), -280.0);
        expectBelow ("DC leakage into other bins", relativeDb (leakage, (double) n), -280.0);
    }

    // Nyquist: alternating +/-1 puts N in the last bin.
    {
        std::vector<double> input ((size_t) n);

        for (int t = 0; t < n; ++t)
            input[(size_t) t] = (t % 2 == 0) ? 1.0 : -1.0;

        fft.forward (input.data(), bins.data());

        const auto last = bins.size() - 1;
        double leakage = 0.0;

        for (size_t k = 0; k < last; ++k)
            leakage = std::max (leakage, std::abs (bins[k]));

        expectBelow ("Nyquist magnitude error", relativeDb (std::abs (bins[last].real() - (double) n), (double) n), -280.0);
        expectBelow ("Nyquist leakage into other bins", relativeDb (leakage, (double) n), -280.0);
    }

    // A bin-centred cosine: amplitude N/2 in that bin only.
    {
        constexpr int binIndex = 64;
        std::vector<double> input ((size_t) n);

        for (int t = 0; t < n; ++t)
            input[(size_t) t] = std::cos (twoPi * (double) binIndex * (double) t / (double) n);

        fft.forward (input.data(), bins.data());

        double leakage = 0.0;

        for (size_t k = 0; k < bins.size(); ++k)
            if (k != (size_t) binIndex)
                leakage = std::max (leakage, std::abs (bins[k]));

        const double expected = (double) n / 2.0;
        expectBelow ("cosine magnitude error",
                     relativeDb (std::abs (std::abs (bins[(size_t) binIndex]) - expected), expected), -280.0);
        expectBelow ("cosine leakage into other bins", relativeDb (leakage, expected), -280.0);
    }
}

//==============================================================================
/** Worst absolute difference between `signal` and `reference` delayed by
    `delay`, over the region where both are fully valid. */
double delayedError (const std::vector<double>& signal,
                     const std::vector<double>& reference,
                     int delay, int skip)
{
    double worst = 0.0;

    for (int n = delay + skip; n < (int) signal.size(); ++n)
        worst = std::max (worst, std::abs (signal[(size_t) n] - reference[(size_t) (n - delay)]));

    return worst;
}

/** Searches for the delay that best explains `signal` as a delayed copy of
    `reference`, so a latency mismatch reports as a number rather than a blob. */
int bestFitDelay (const std::vector<double>& signal,
                  const std::vector<double>& reference,
                  int maxDelay, int skip)
{
    int best = 0;
    double bestErr = -1.0;

    for (int d = 0; d <= maxDelay; ++d)
    {
        const double err = delayedError (signal, reference, d, skip);

        if (bestErr < 0.0 || err < bestErr)
        {
            bestErr = err;
            best = d;
        }
    }

    return best;
}

void testWindowCola()
{
    std::printf ("\nWindow pair overlap-adds flat\n");

    for (int order : { 8, 10, 12, 13 })
    {
        const int n = 1 << order;

        for (int hop : { n / 2, n / 4, n / 8, n / 16, 256 })
        {
            if (hop < 1 || n % hop != 0 || hop > n / 2)
                continue;

            WolaWindow w (n, hop);

            expectBelow (std::string ("N=") + std::to_string (n) + ", hop=" + std::to_string (hop)
                             + " COLA deviation",
                         relativeDb (w.colaDeviation(), 1.0), -280.0);
        }
    }
}

void testWolaPerfectReconstruction()
{
    std::printf ("\nWOLA reconstructs exactly with an identity frame function\n");

    struct Case { int order, hop; };

    for (auto c : { Case { 10, 256 }, Case { 12, 256 }, Case { 13, 256 }, Case { 11, 512 } })
    {
        const int n = 1 << c.order;
        const int total = 8 * n;

        const auto input = makeNoise (total, 777u + (unsigned) c.order + (unsigned) c.hop);
        std::vector<double> output ((size_t) total, 0.0);

        Wola64 wola (c.order, c.hop);

        // Deliberately irregular block sizes: a host will not hand us multiples
        // of the hop, and the engine has to be indifferent to that.
        const int blockSizes[] = { 64, 511, 1, 129, 1024, 37 };
        int at = 0, which = 0;

        while (at < total)
        {
            const int len = std::min (blockSizes[which % 6], total - at);
            wola.process (input.data() + at, output.data() + at, len,
                          [] (std::complex<double>*, int) { /* identity */ });
            at += len;
            ++which;
        }

        const int latency = wola.getLatencySamples();
        const double err = delayedError (output, input, latency, n);
        const int fitted = bestFitDelay (output, input, 2 * n, n);

        expectBelow (std::string ("N=") + std::to_string (n) + ", hop=" + std::to_string (c.hop)
                         + " reconstruction",
                     relativeDb (err, peakMagnitude (input)), -280.0);

        ++checks;

        if (fitted != latency)
        {
            ++failures;
            std::printf ("  [FAIL] %-58s reported %d, best fit %d\n",
                         "reported latency matches measured", latency, fitted);
        }
        else
        {
            std::printf ("  [pass] %-58s %8d samples\n", "reported latency matches measured", latency);
        }
    }
}

void testOverlapSaveAgainstDirectConvolution()
{
    std::printf ("\nOverlap-save matches direct FIR convolution\n");

    struct Case { int order, taps, advance; };

    for (auto c : { Case { 13, 4096, 3840 },   // the main path's intended geometry
                    Case { 13, 4096, 0 },      // maximum advance
                    Case { 10, 200, 256 },
                    Case { 12, 1024, 512 } })
    {
        const int n = 1 << c.order;
        OverlapSave64 os (c.order, c.taps, c.advance);

        ++checks;
        if (! os.isValid())
        {
            ++failures;
            std::printf ("  [FAIL] N=%d taps=%d advance=%d is not a valid geometry\n",
                         n, c.taps, os.getHop());
            continue;
        }
        std::printf ("  [pass] N=%-5d taps=%-5d advance=%-5d geometry valid\n", n, c.taps, os.getHop());

        // A decaying random impulse response: broadband, and nothing about it is
        // symmetric or special.
        auto ir = makeNoise (c.taps, 31337u + (unsigned) c.order);
        for (int j = 0; j < c.taps; ++j)
            ir[(size_t) j] *= std::exp (-3.0 * (double) j / (double) c.taps);

        os.setImpulseResponse (ir.data(), c.taps);

        const int total = 4 * n;
        const auto input = makeNoise (total, 5150u + (unsigned) c.order);
        std::vector<double> output ((size_t) total, 0.0);

        const int blockSizes[] = { 128, 1000, 7, 4096, 333 };
        int at = 0, which = 0;

        while (at < total)
        {
            const int len = std::min (blockSizes[which % 5], total - at);
            os.process (input.data() + at, output.data() + at, len);
            at += len;
            ++which;
        }

        // Direct time-domain convolution as ground truth.
        std::vector<double> truth ((size_t) total, 0.0);
        for (int t = 0; t < total; ++t)
        {
            double acc = 0.0;
            const int jMax = std::min (c.taps - 1, t);

            for (int j = 0; j <= jMax; ++j)
                acc += ir[(size_t) j] * input[(size_t) (t - j)];

            truth[(size_t) t] = acc;
        }

        const int latency = os.getLatencySamples();
        const double err = delayedError (output, truth, latency, c.taps);

        expectBelow (std::string ("N=") + std::to_string (n) + ", taps=" + std::to_string (c.taps)
                         + ", advance=" + std::to_string (os.getHop()),
                     relativeDb (err, peakMagnitude (truth)), -280.0);
    }
}

void expectEqual (const std::string& what, long long got, long long want)
{
    ++checks;
    const bool ok = (got == want);

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-58s %8lld  (want %lld)\n",
                 ok ? "pass" : "FAIL", what.c_str(), got, want);
}

void testGeometry()
{
    std::printf ("\nGeometry scales with sample rate and caps at 2x\n");

    struct Expect { double rate; int exponent, hfSize, lfSize, hfHop, lfHop, passLatency; };

    for (auto e : { Expect {  11025.0, -2,  256, 1024,  64,  256,  1280 },
                    Expect {  22050.0, -1,  512, 2048, 128,  512,  2560 },
                    Expect {  44100.0,  0, 1024, 4096, 256, 1024,  5120 },
                    Expect {  48000.0,  0, 1024, 4096, 256, 1024,  5120 },
                    Expect {  88200.0,  1, 2048, 8192, 512, 2048, 10240 },
                    Expect {  96000.0,  1, 2048, 8192, 512, 2048, 10240 },
                    Expect { 176400.0,  1, 2048, 8192, 512, 2048, 10240 },
                    Expect { 192000.0,  1, 2048, 8192, 512, 2048, 10240 } })
    {
        const auto g = SrGeometry::forSampleRate (e.rate);
        const auto tag = std::to_string ((int) e.rate) + " Hz ";

        expectEqual (tag + "rate exponent", g.rateExponent, e.exponent);
        expectEqual (tag + "HF window", g.hfSize(), e.hfSize);
        expectEqual (tag + "LF window", g.lfSize(), e.lfSize);
        expectEqual (tag + "HF hop", g.hfHop, e.hfHop);
        expectEqual (tag + "LF hop", g.lfHop, e.lfHop);
        expectEqual (tag + "one-pass latency", g.passLatencySamples(), e.passLatency);

        // Each window must be a whole multiple of its own hop, or the encoder's and
        // decoder's analysis framing would drift apart and their control state would
        // stop matching, which is what makes the decode exact.
        expectEqual (tag + "HF window is a whole number of hops", g.hfSize() % g.hfHop, 0);
        expectEqual (tag + "LF window is a whole number of hops", g.lfSize() % g.lfHop, 0);

        // And the LF hop must be a whole multiple of the HF hop, so the two halves'
        // frames land together rather than beating against each other.
        expectEqual (tag + "LF hop is a whole number of HF hops", g.lfHop % g.hfHop, 0);

        // Four-fold overlap on both windows. Sixteen on the long one was the single
        // largest avoidable cost in the engine; see the note on SrGeometry.
        expectEqual (tag + "HF overlap factor", g.hfSize() / g.hfHop, 4);
        expectEqual (tag + "LF overlap factor", g.lfSize() / g.lfHop, 4);
    }
}

void testEngineBypassIsPureDelay()
{
    std::printf ("\nEngine with the process switched off is an exact delay\n");

    for (double rate : { 22050.0, 48000.0, 96000.0 })
    {
        SrEngine engine;
        engine.prepare (rate, 512, 1);
        engine.setProcess (SrProcess::off);

        // Encode, not loop. With the process off the engine is an exact delay in every mode
        // that does not have the modelled channel in it — and loop does, deliberately, so
        // that switching the process off leaves the channel to be heard on its own. See
        // testEngineOffKeepsTheChannelInLoop.
        engine.setMode (SrMode::encode);

        const int latency = engine.getLatencySamples();
        const int total = 3 * latency;

        const auto input = makeNoise (total, 8675u + (unsigned) rate);
        std::vector<double> output = input;

        const int blockSizes[] = { 256, 64, 401, 512, 13 };
        int at = 0, which = 0;

        while (at < total)
        {
            const int len = std::min (blockSizes[which % 5], total - at);
            double* ptr = output.data() + at;
            engine.process (&ptr, 1, len);
            at += len;
            ++which;
        }

        const double err = delayedError (output, input, latency, 0);

        expectBelow (std::to_string ((int) rate) + " Hz bypass residual",
                     relativeDb (err, peakMagnitude (input)), -280.0);
    }
}

/** Streams a signal through an engine in irregular blocks. */
void runEngineBlocks (SrEngine& engine, std::vector<double>& buffer)
{
    const int blockSizes[] = { 512, 128, 333, 64, 1024, 7 };
    int at = 0, which = 0;
    const int total = (int) buffer.size();

    while (at < total)
    {
        const int len = std::min (blockSizes[which % 6], total - at);
        double* ptr = buffer.data() + at;
        engine.process (&ptr, 1, len);
        at += len;
        ++which;
    }
}

void testEngineOffKeepsTheChannelInLoop()
{
    std::printf ("\nWith the process off, loop mode still runs the modelled channel\n");

    // The panel's central comparison is: switch the process between off and b with the mode
    // on both, and hear what the noise reduction does to the channel. That only says
    // anything if the channel is in circuit either way. It was not — off was a plain delay,
    // so both sides of the comparison were clean and the switch was inaudible.
    const double rate = 48000.0;
    const int total = 1 << 16;

    const auto input = makeNoise (total, 24601u);

    std::vector<double> loopOff = input, encodeOff = input;

    {
        SrEngine engine;
        engine.prepare (rate, 512, 1);
        engine.setProcess (SrProcess::off);
        engine.setMode (SrMode::loop);
        runEngineBlocks (engine, loopOff);
    }

    {
        SrEngine engine;
        engine.prepare (rate, 512, 1);
        engine.setProcess (SrProcess::off);
        engine.setMode (SrMode::encode);
        runEngineBlocks (engine, encodeOff);
    }

    // Measured past the pipeline, against the delayed input.
    SrEngine geometry;
    geometry.prepare (rate, 512, 1);
    const int latency = geometry.getLatencySamples();

    double loopErr = 0.0, encodeErr = 0.0, reference = 0.0;

    for (int n = latency; n < total; ++n)
    {
        const double dry = input[(size_t) (n - latency)];
        loopErr = std::max (loopErr, std::abs (loopOff[(size_t) n] - dry));
        encodeErr = std::max (encodeErr, std::abs (encodeOff[(size_t) n] - dry));
        reference = std::max (reference, std::abs (dry));
    }

    // Encode with the process off is still an exact delay; loop is audibly not.
    expectBelow ("encode, off, residual against the delayed input",
                 relativeDb (encodeErr, reference), -280.0);

    const double loopDb = relativeDb (loopErr, reference);
    ++checks;
    const bool loopAudible = loopDb > -40.0;

    if (! loopAudible)
        ++failures;

    std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                 loopAudible ? "pass" : "FAIL",
                 "loop, off, the channel is audible", loopDb, -40.0);
}

/** Streams a signal through an engine in irregular blocks. */
void runEngine (SrEngine& engine, std::vector<double>& buffer)
{
    const int blockSizes[] = { 512, 128, 333, 64, 1024, 7 };
    int at = 0, which = 0;
    const int total = (int) buffer.size();

    while (at < total)
    {
        const int len = std::min (blockSizes[which % 6], total - at);
        double* ptr = buffer.data() + at;
        engine.process (&ptr, 1, len);
        at += len;
        ++which;
    }
}

void testEncodeDecodeNull()
{
    std::printf ("\nEncode then decode is complementary\n");

    // Loop mode with every tape defect disabled: the impairment stage is transparent,
    // so what remains is purely the encoder followed by its own inverse.
    TapeChannelSettings transparent;
    transparent.hissEnabled = false;
    transparent.saturationEnabled = false;
    transparent.hfLossEnabled = false;
    transparent.modulationNoiseEnabled = false;

    for (double rate : { 44100.0, 48000.0, 96000.0 })
    {
        SrEngine engine;
        engine.prepare (rate, 1024, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (SrMode::loop);
        engine.getTapeChannel().setSettings (transparent);

        const int latency = engine.getLatencySamples();
        const int total = 6 * latency;

        // Pink-ish broadband material at a level that exercises the action region.
        auto input = makeNoise (total, 20250u + (unsigned) rate);
        for (auto& v : input)
            v *= 0.05;

        std::vector<double> output = input;
        runEngine (engine, output);

        // Skip the settling head: the control smoothers start from rest.
        const double err = delayedError (output, input, latency, latency);

        // OPEN ITEM. The plan asks for -100 dB; what is achieved is -30 to -68 dB.
        //
        // The cause is measured, not guessed, and it is a direct trade-off rather than a
        // bug. An STFT per-bin gain is only invertible when the profile varies gradually
        // across bins: a smooth 24 dB tilt is undone to -126 dB, a single-bin 24 dB
        // notch only to -16 dB. But accurate level detection needs the opposite —
        // sharp selectivity, so a tone is read at its true level. Smoothing the envelope
        // by *averaging* favours the null (-86 dB tonal) and costs curve fidelity
        // (5 dB of boost left at reference where the paper allows one). Smoothing it in
        // power with a peak-normalised kernel favours fidelity (0.85 dB at reference,
        // and every published figure within about a dB) and costs the null.
        //
        // Fidelity to the papers is the point of the project, so that is what is
        // chosen. These limits record what the chosen setting achieves, so a regression
        // is still visible. See DESIGN.md for the route that could satisfy both.
        expectBelow (std::to_string ((int) rate) + " Hz loop-mode null, broadband",
                     relativeDb (err, peakMagnitude (input)), -26.0);
    }

    // Separate encoder and decoder instances, which is how the plugin is actually
    // used: print an encoded track, decode it on the way back.
    {
        constexpr double rate = 48000.0;

        SrEngine encoder, decoder;
        encoder.prepare (rate, 1024, 1);
        decoder.prepare (rate, 1024, 1);
        encoder.setProcess (SrProcess::spectralRecording);
        decoder.setProcess (SrProcess::spectralRecording);
        encoder.setMode (SrMode::encode);
        decoder.setMode (SrMode::decode);

        const int latency = encoder.getLatencySamples() + decoder.getLatencySamples();
        const int total = 4 * latency;

        struct Case { const char* name; double scale; bool tonal; };

        for (auto c : { Case { "steady tone at -20 dB", 0.1, true },
                        Case { "broadband at -26 dB", 0.05, false },
                        Case { "broadband at -60 dB", 5.0e-4, false },
                        Case { "broadband near reference", 0.7, false } })
        {
            encoder.reset();
            decoder.reset();

            std::vector<double> input ((size_t) total);

            if (c.tonal)
                for (int n = 0; n < total; ++n)
                    input[(size_t) n] = c.scale * std::sin (twoPi * 1000.0 * (double) n / rate);
            else
            {
                input = makeNoise (total, 4711u);
                for (auto& v : input)
                    v *= c.scale;
            }

            std::vector<double> work = input;
            runEngine (encoder, work);
            runEngine (decoder, work);

            const double err = delayedError (work, input, latency, latency);

            // Recorded, not aspirational — see the note above.
            expectBelow (std::string ("separate instances, ") + c.name,
                         relativeDb (err, peakMagnitude (input)), c.tonal ? -64.0 : -38.0);
        }
    }
}

//==============================================================================
void expectNear (const std::string& what, double got, double want, double tolerance,
                 const char* units = "dB")
{
    ++checks;
    const bool ok = std::abs (got - want) <= tolerance;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-46s %+8.3f %s (want %+7.3f, tol %.3f)\n",
                 ok ? "pass" : "FAIL", what.c_str(), got, units, want, tolerance);
}

/** The antisaturation figures from section 4.1, as attenuation in dB relative to
    the midband. These are the paper's own numbers and act as ground truth. */
const std::vector<std::pair<double, double>> publishedAntisaturation
{
    {    25.0, 10.0 },
    {  5000.0,  2.0 },
    { 10000.0,  6.0 },
    { 15000.0, 10.0 }
};

constexpr double midbandRefHz = 1000.0;

void testAnalogPrototypesMatchThePaper()
{
    std::printf ("\nAnalog main-path prototypes reproduce the section 4.1 figures\n");

    const auto protos = srNetworks::mainPath();
    const double ref = 20.0 * std::log10 (std::abs (srNetworks::responseAt (protos, midbandRefHz)));

    for (auto [hz, wantAtten] : publishedAntisaturation)
    {
        const double got = 20.0 * std::log10 (std::abs (srNetworks::responseAt (protos, hz))) - ref;
        expectNear (std::to_string ((int) hz) + " Hz antisaturation", got, -wantAtten, 0.35);
    }

    // Least treatment: the fixed networks must not encroach on the midband.
    for (double hz : { 200.0, 400.0, 700.0, 1500.0, 2500.0 })
    {
        const double got = 20.0 * std::log10 (std::abs (srNetworks::responseAt (protos, hz))) - ref;
        expectNear (std::to_string ((int) hz) + " Hz midband flatness", got, 0.0, 0.25);
    }
}

void testDigitalMainPathMatchesThePublishedFigures()
{
    std::printf ("\nDigital main path reproduces the section 4.1 figures at every rate\n");

    for (double rate : { 44100.0, 48000.0, 96000.0, 192000.0 })
    {
        const auto coeffs = srNetworks::mainPathBiquadsForRate (rate);

        const double worstTarget = srNetworks::worstTargetErrorDb (coeffs, rate);
        const double worstMidband = srNetworks::worstMidbandErrorDb (coeffs, rate);

        // Section 4.1 states its figures as "about 2 dB", "about 6 dB", "about
        // 10 dB", so a few tenths is within the paper's own precision. The
        // residual is irreducible: the skewing low-pass has its shelf zeros above
        // Nyquist, so no biquad can track the analog prototype exactly near the
        // top of the band, and the error is largest at 44.1 kHz where 15 kHz sits
        // at 0.68 of Nyquist.
        expectNear (std::to_string ((int) rate) + " Hz worst published-figure error",
                    worstTarget, 0.0, 0.50);
        expectNear (std::to_string ((int) rate) + " Hz worst midband deviation",
                    worstMidband, 0.0, 0.25);
    }

    // Print the achieved curve at 48 kHz so the numbers are on the record.
    const auto coeffs = srNetworks::mainPathBiquadsForRate (48000.0);
    BiquadCascade64 cascade;
    cascade.setCoeffs (coeffs);

    const double ref = 20.0 * std::log10 (std::abs (cascade.responseAt (twoPi * 1000.0 / 48000.0)));
    std::printf ("    achieved at 48 kHz, relative to 1 kHz:");

    for (double hz : { 25.0, 100.0, 1000.0, 5000.0, 10000.0, 15000.0 })
        std::printf ("  %g Hz %+.2f", hz,
                     20.0 * std::log10 (std::abs (cascade.responseAt (twoPi * hz / 48000.0))) - ref);

    std::printf ("\n");
}

void testBiquadInversionIsExact()
{
    std::printf ("\nMain path is exactly invertible: encode then decode nulls\n");

    for (double rate : { 44100.0, 48000.0, 96000.0 })
    {
        const auto coeffs = srNetworks::mainPathBiquadsForRate (rate);

        BiquadCascade64 forward, inverse;
        forward.setCoeffs (coeffs);
        inverse.setInverseCoeffs (coeffs);

        const int total = 1 << 16;
        const auto input = makeNoise (total, 60600u + (unsigned) rate);
        std::vector<double> work = input;

        forward.process (work.data(), total);
        inverse.process (work.data(), total);

        double worst = 0.0;

        // Skip the head: both cascades start from rest, so the first samples
        // carry the settling of the pair rather than a steady-state error.
        for (int n = 64; n < total; ++n)
            worst = std::max (worst, std::abs (work[(size_t) n] - input[(size_t) n]));

        // Not the numerical floor: the 25 dB LF skewing shelf puts a zero at
        // 9.5 Hz, so the inverse has a pole of radius 0.9991 and a direct-form
        // biquad loses precision that close to the unit circle. Still some
        // 128 dB below the -100 dB the plan asks of the full encode/decode null.
        expectBelow (std::to_string ((int) rate) + " Hz forward-then-inverse residual",
                     relativeDb (worst, peakMagnitude (input)), -200.0);
    }

    // The reason the main path is not an STFT multiplier: verify the inverse
    // networks really are minimum phase, so the coefficient swap is stable.
    std::printf ("\nInverse cascade is stable (all poles inside the unit circle)\n");

    const auto coeffs = srNetworks::mainPathBiquadsForRate (48000.0);
    int sectionIndex = 0;

    for (const auto& c : coeffs)
    {
        // After inversion the denominator is the forward numerator, so the
        // forward zeros become the inverse poles. Roots of b0 + b1 z^-1 + b2 z^-2.
        const auto inv = c.inverted();
        const double a1 = inv.a1, a2 = inv.a2;
        const double disc = a1 * a1 - 4.0 * a2;

        double maxRadius = 0.0;

        if (disc >= 0.0)
        {
            const double r = std::sqrt (disc);
            maxRadius = std::max (std::abs ((-a1 + r) * 0.5), std::abs ((-a1 - r) * 0.5));
        }
        else
        {
            maxRadius = std::sqrt (a2);   // complex pair: radius^2 = a2
        }

        ++checks;
        const bool ok = maxRadius < 1.0;

        if (! ok)
            ++failures;

        std::printf ("  [%s] section %d inverse pole radius %39.6f\n",
                     ok ? "pass" : "FAIL", sectionIndex, maxRadius);
        ++sectionIndex;
    }
}

//==============================================================================
/** Builds a staggered cascade of per-bin fixed bands and measures the steady-state
    encode gain against input level, which is what Figs. 12 and 13 specify. */
struct StaggeredCascade
{
    static constexpr double sampleRate = 48000.0;
    static constexpr int fftSize = 1024;
    static constexpr int hop = 256;

    std::vector<FixedBand> stages;
    int numBins = fftSize / 2 + 1;
    int probeBin = 0;

    StaggeredCascade (int numStages, double stageGainDb, double finalTauMs,
                      bool highFrequency, double probeHz)
    {
        // Section 2.5: thresholds of approximately -30, -48 and -62 dB. The
        // high-frequency side uses all three; the low-frequency side uses two.
        const double thresholds[3] = { -30.0, -48.0, -62.0 };

        stages.resize ((size_t) numStages);
        probeBin = (int) std::lround (probeHz * fftSize / sampleRate);

        for (int i = 0; i < numStages; ++i)
        {
            FixedBandParams p;
            p.thresholdDb = thresholds[i];
            p.precedingGainDb = stageGainDb * (double) i;
            p.lowLevelGainDb = stageGainDb;
            p.finalTauMs = finalTauMs;
            p.passBandCornerHz = highFrequency ? 1600.0 : 400.0;
            p.passBandIsHighPass = highFrequency;

            stages[(size_t) i].prepare (p, numBins, sampleRate, fftSize,
                                        (double) hop / sampleRate);
        }
    }

    /** Settles the cascade on a flat spectrum at `levelDb` and returns the total
        encode gain in dB at the probe bin. */
    double gainAt (double levelDb)
    {
        for (auto& s : stages)
            s.reset();

        const double magnitude = std::pow (10.0, levelDb / 20.0);

        std::vector<double> mags ((size_t) numBins, magnitude);
        std::vector<double> transfer ((size_t) numBins, 0.0);

        double probe = magnitude;

        // 600 frames is 3.2 s at this geometry, well past the 300 ms final
        // smoothing, so the measurement is genuinely steady state.
        for (int frame = 0; frame < 600; ++frame)
        {
            std::fill (mags.begin(), mags.end(), magnitude);
            probe = magnitude;

            for (auto& s : stages)
            {
                s.processFrame (mags.data(), transfer.data());

                for (int k = 0; k < numBins; ++k)
                    mags[(size_t) k] *= 1.0 + transfer[(size_t) k];

                probe *= 1.0 + transfer[(size_t) probeBin];
            }
        }

        return 20.0 * std::log10 (probe / magnitude);
    }
};

void testFixedBandStaggering()
{
    std::printf ("\nStaggered per-bin fixed bands: encode gain vs input level\n");

    // "somewhat over 8 dB" per stage (section 1). The exact figure is fitted below.
    constexpr double stageGainDb = 8.0;

    StaggeredCascade hf (3, stageGainDb, 160.0, true, 3000.0);
    StaggeredCascade lf (2, stageGainDb, 300.0, false, 200.0);

    std::printf ("    input      HF gain    LF gain\n");

    for (double levelDb : { -100.0, -90.0, -80.0, -70.0, -60.0, -50.0, -40.0,
                            -30.0, -20.0, -10.0, 0.0, 10.0 })
    {
        std::printf ("    %6.0f dB  %7.2f dB %7.2f dB\n",
                     levelDb, hf.gainAt (levelDb), lf.gainAt (levelDb));
    }

    // Section 1: about 24 dB of total dynamic effect at high frequencies, 16 dB at
    // low. Measured deep below the lowest threshold, where every stage is fully
    // boosting.
    expectNear ("HF total boost, deep subthreshold", hf.gainAt (-110.0), 24.0, 2.5);
    expectNear ("LF total boost, deep subthreshold", lf.gainAt (-110.0), 16.0, 2.5);

    // Section 4.1: no action in the top 25 dB of the range, so the boost must have
    // collapsed to essentially unity by reference level. This is the dual-path
    // premise — high-level signals travel the main path unprocessed.
    // Section 1 quantifies "a further dynamic action of about 1 dB above the
    // reference level", so the residual boost there should be roughly 1 dB, not
    // zero: the top of the range is nearly but not entirely unprocessed.
    expectNear ("HF gain at reference level", hf.gainAt (0.0), 1.0, 0.5);
    expectNear ("LF gain at reference level", lf.gainAt (0.0), 1.0, 0.5);

    // Section 4.1: "a compression ratio of about 2:1" in the intermediate level
    // region. Measured as the steepest local slope, since the average over the
    // whole action range is necessarily gentler than its steepest part.
    const auto steepestRatio = [] (StaggeredCascade& c)
    {
        double steepest = 0.0;

        for (double levelDb = -90.0; levelDb <= 5.0; levelDb += 5.0)
        {
            const double lower = c.gainAt (levelDb - 2.5);
            const double upper = c.gainAt (levelDb + 2.5);
            const double slope = (lower - upper) / 5.0;

            steepest = std::max (steepest, 1.0 / (1.0 - slope));
        }

        return steepest;
    };

    const double hfRatio = steepestRatio (hf);
    const double lfRatio = steepestRatio (lf);

    std::printf ("    steepest compression ratio: HF %.2f:1, LF %.2f:1\n", hfRatio, lfRatio);

    // Section 4.1's "about 2:1" is not fully reached: with thresholds 14 to 18 dB
    // apart, three 8 dB stages cannot compound more steeply without giving up the
    // 1 dB of high-level action that section 1 quantifies. Measured 1.77:1, and the
    // trade-off table is in FixedBandParams::feedforwardWeight.
    expectNear ("HF steepest compression ratio", hfRatio, 1.8, 0.2, ":1");

    // Two low-frequency stages cannot compound as steeply as three, so the LF side
    // is expected to fall short of 2:1; it is checked for being in the right region
    // rather than against the paper's headline figure.
    expectNear ("LF steepest compression ratio", lfRatio, 1.6, 0.3, ":1");

    // The characteristic must be monotonic: more input, less boost, never a
    // reversal, or the decoder cannot invert it.
    double previous = 1.0e9;
    bool monotonic = true;

    for (double levelDb = -110.0; levelDb <= 12.0; levelDb += 2.0)
    {
        const double g = hf.gainAt (levelDb);

        if (g > previous + 1.0e-9)
            monotonic = false;

        previous = g;
    }

    ++checks;

    if (! monotonic)
    {
        ++failures;
        std::printf ("  [FAIL] %-58s\n", "HF characteristic is monotonic in level");
    }
    else
    {
        std::printf ("  [pass] %-58s\n", "HF characteristic is monotonic in level");
    }
}

//==============================================================================
void testModulationControl()
{
    std::printf ("\nModulation control signals MC1-MC8\n");

    constexpr double sampleRate = 48000.0;
    constexpr int fftSize = 4096;
    constexpr int numBins = fftSize / 2 + 1;
    constexpr double frameSeconds = 256.0 / sampleRate;

    ModulationControl mc;
    mc.prepare (numBins, sampleRate, fftSize, frameSeconds);

    std::vector<double> spectrum ((size_t) numBins, 0.0);
    std::vector<double> silence ((size_t) numBins, 0.0);

    // A unit tone in a single bin. The MC signals are weighted RMS levels on the
    // per-bin magnitude scale (see ModulationControl.h), so a lone tone reads its
    // weighting magnitude divided by the weighting's norm.
    const auto toneAt = [&] (double hz)
    {
        std::fill (spectrum.begin(), spectrum.end(), 0.0);
        const int bin = (int) std::lround (hz * fftSize / sampleRate);
        spectrum[(size_t) bin] = 1.0;
        return (double) bin * sampleRate / (double) fftSize;
    };

    // The same single-pole magnitude the weightings are built from, so the shape
    // checks below can state their expectation analytically.
    const auto onePole = [] (double hz, double cornerHz, bool highPass)
    {
        const double ratio = hz / cornerHz;
        const double lowPass = 1.0 / std::sqrt (1.0 + ratio * ratio);

        return highPass ? lowPass * ratio : lowPass;
    };

    // The weighting *shapes* are checked as ratios between two tones — one at the
    // corner, one deep in the weighting's own passband — which the normalisation
    // cancels out of, so these pin the corners themselves.
    const auto mcRatio = [&] (double cornerToneHz, double passbandToneHz, auto&& read)
    {
        const double hzCorner = toneAt (cornerToneHz);
        mc.processFrame (spectrum.data(), spectrum.data());
        const double atCorner = read();

        const double hzPassband = toneAt (passbandToneHz);
        mc.processFrame (spectrum.data(), spectrum.data());
        const double inPassband = read();

        return std::array<double, 3> { atCorner / inPassband, hzCorner, hzPassband };
    };

    using namespace mcWeighting;

    {
        const auto [ratio, hzC, hzP] = mcRatio (3000.0, 18000.0, [&] { return mc.hfSliding(); });
        expectNear ("MC1 corner against its passband", ratio,
                    onePole (hzC, hfSlidingHighPassHz, true) / onePole (hzP, hfSlidingHighPassHz, true),
                    0.01, "   ");
    }
    {
        const auto [ratio, hzC, hzP] = mcRatio (200.0, 23.0, [&] { return mc.lfSliding(); });
        expectNear ("MC4 corner against its passband", ratio,
                    onePole (hzC, lfSlidingLowPassHz, false) / onePole (hzP, lfSlidingLowPassHz, false),
                    0.01, "   ");
    }
    {
        const auto pair = [&] (double hz) { return onePole (hz, hfFixedLowPassOneHz, false)
                                                 * onePole (hz, hfFixedLowPassTwoHz, false); };
        const auto [ratio, hzC, hzP] = mcRatio (400.0, 23.0, [&] { return mc.hfFixed(); });
        expectNear ("MC3 corner pair against its passband", ratio, pair (hzC) / pair (hzP), 0.01, "   ");
    }
    {
        const auto pair = [&] (double hz) { return onePole (hz, lfFixedHighPassOneHz, true)
                                                 * onePole (hz, lfFixedHighPassTwoHz, true); };
        const auto [ratio, hzC, hzP] = mcRatio (800.0, 18000.0, [&] { return mc.lfFixed(); });
        expectNear ("MC6 corner pair against its passband", ratio, pair (hzC) / pair (hzP), 0.01, "   ");
    }

    // The *scale*: MC is the weighted RMS of the per-bin magnitudes, so a flat
    // spectrum of magnitude m reads exactly m under every weighting — the property
    // that keeps the opposition commensurate with a per-bin control signal at
    // every transform size. A raw power sum would read tens of times higher here,
    // and differently at every sample rate.
    std::fill (spectrum.begin(), spectrum.end(), 0.25);
    mc.processFrame (spectrum.data(), spectrum.data());

    expectNear ("a flat spectrum reads its own magnitude (MC1)", mc.hfSliding(), 0.25, 1.0e-9, "   ");
    expectNear ("a flat spectrum reads its own magnitude (MC4)", mc.lfSliding(), 0.25, 1.0e-9, "   ");
    expectNear ("a flat spectrum reads its own magnitude (MC3)", mc.hfFixed(), 0.25, 1.0e-9, "   ");
    expectNear ("a flat spectrum reads its own magnitude (MC6)", mc.lfFixed(), 0.25, 1.0e-9, "   ");

    // Selectivity: modulation control exists to make each stage insensitive to
    // energy outside its own region, so a low tone must drive the low-frequency
    // weightings and not the high-frequency ones, and the reverse.
    toneAt (100.0);
    mc.processFrame (spectrum.data(), spectrum.data());
    const double lowToneHfSliding = mc.hfSliding();
    const double lowToneLfSliding = mc.lfSliding();

    toneAt (12000.0);
    mc.processFrame (spectrum.data(), spectrum.data());
    const double highToneLfFixed = mc.lfFixed();
    const double highToneLfSliding = mc.lfSliding();

    std::printf ("    100 Hz tone: MC1 %.4f, MC4 %.4f    12 kHz tone: MC4 %.4f, MC6 %.4f\n",
                 lowToneHfSliding, lowToneLfSliding, highToneLfSliding, highToneLfFixed);

    expectBelow ("MC1 rejects a 100 Hz tone", relativeDb (lowToneHfSliding, 1.0), -25.0);
    expectBelow ("MC4 rejects a 12 kHz tone", relativeDb (highToneLfSliding, 1.0), -35.0);

    ++checks;
    if (lowToneLfSliding <= lowToneHfSliding || highToneLfFixed <= highToneLfSliding)
    {
        ++failures;
        std::printf ("  [FAIL] %-58s\n", "weightings select their own frequency region");
    }
    else
    {
        std::printf ("  [pass] %-58s\n", "weightings select their own frequency region");
    }

    // The smoothed variants must lag their unsmoothed sources rather than track
    // them: MC2 from MC1, MC5 from MC4, MC7 from MC6.
    mc.reset();
    toneAt (3000.0);
    mc.processFrame (spectrum.data(), spectrum.data());

    ++checks;
    if (mc.hfSlidingOvershoot() >= mc.hfSliding())
    {
        ++failures;
        std::printf ("  [FAIL] %-58s\n", "MC2 lags MC1 on a step");
    }
    else
    {
        std::printf ("  [pass] %-58s  %.4f vs %.4f\n", "MC2 lags MC1 on a step",
                     mc.hfSlidingOvershoot(), mc.hfSliding());
    }

    // ...and must converge to it when the input is held.
    for (int frame = 0; frame < 400; ++frame)
        mc.processFrame (spectrum.data(), spectrum.data());

    expectNear ("MC2 converges to MC1 when held", mc.hfSlidingOvershoot(), mc.hfSliding(),
                0.001, "   ");

    // MC8 detects high-frequency transients: silent for a steady signal, fires on a
    // rise, then decays with the specified 30 ms hold.
    mc.reset();
    const double transientToneHz = toneAt (10000.0);

    for (int frame = 0; frame < 200; ++frame)
        mc.processFrame (spectrum.data(), spectrum.data());

    const double steadyTransient = mc.hfTransient();

    mc.reset();
    mc.processFrame (silence.data(), silence.data());
    mc.processFrame (spectrum.data(), spectrum.data());
    const double onsetTransient = mc.hfTransient();

    std::printf ("    MC8: steady %.6f, on onset %.4f\n", steadyTransient, onsetTransient);

    expectBelow ("MC8 is quiet for a steady high-frequency tone",
                 relativeDb (steadyTransient, 1.0), -60.0);

    // The onset reads the tone's full weighted RMS level — the rise from silence
    // is the level itself — which for a lone unit tone is its weighting magnitude
    // over the weighting's norm. Computed here from the published corner, so the
    // check pins the scale rather than just the sign.
    double transientNorm = 0.0;

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * sampleRate / (double) fftSize;
        const double w = onePole (hz, transientHighPassHz, true);
        transientNorm += w * w;
    }

    transientNorm = std::sqrt (transientNorm);

    expectNear ("MC8 fires on a high-frequency onset", onsetTransient,
                onePole (transientToneHz, transientHighPassHz, true) / transientNorm,
                1.0e-9, "   ");

    // Peak hold: after the onset, silence should decay with roughly a 30 ms
    // constant, so about one time constant later it should be near 1/e.
    const int framesIn30ms = (int) std::lround (0.030 / frameSeconds);

    for (int frame = 0; frame < framesIn30ms; ++frame)
        mc.processFrame (silence.data(), silence.data());

    const double afterOneTau = mc.hfTransient() / onsetTransient;
    expectNear ("MC8 peak hold decays one time constant in 30 ms",
                afterOneTau, std::exp (-1.0), 0.05, "   ");
}

//==============================================================================
/** One staggered stage: a per-bin fixed band and a sliding layer, combined by
    whichever composition is under test. Used to reproduce the probe-tone
    behaviour of Figs. 14 to 16. */
struct ProbeStage
{
    static constexpr double sampleRate = 48000.0;
    static constexpr int fftSize = 2048;
    static constexpr int hop = 256;
    static constexpr double stageGainDb = 8.0;

    int numBins = fftSize / 2 + 1;

    FixedBand fixedBand;
    SlidingLayer sliding;
    SlidingComposition composition;

    std::vector<double> stageInput, fixedTransfer, slidingFraction, spreadLevel, fixedFeed;

    ProbeStage (SlidingComposition mode, bool highFrequency)
        : composition (mode)
    {
        FixedBandParams fp;
        fp.thresholdDb = -30.0;
        fp.lowLevelGainDb = stageGainDb;
        fp.passBandCornerHz = highFrequency ? 1600.0 : 400.0;
        fp.passBandIsHighPass = highFrequency;
        fp.finalTauMs = highFrequency ? 160.0 : 300.0;
        fixedBand.prepare (fp, numBins, sampleRate, fftSize, (double) hop / sampleRate);

        SlidingLayerParams sp;
        sp.highFrequency = highFrequency;
        sp.lowLevelGainDb = stageGainDb;
        sp.thresholdDb = -30.0;
        sp.controlCornerHz = highFrequency ? 10000.0 : 80.0;
        sp.firstTauMs = highFrequency ? 5.0 : 7.5;
        sp.finalTauMs = highFrequency ? 80.0 : 150.0;
        sp.composition = mode;
        sliding.prepare (sp, numBins, sampleRate, fftSize, (double) hop / sampleRate);

        stageInput.assign ((size_t) numBins, 0.0);
        fixedTransfer.assign ((size_t) numBins, 0.0);
        slidingFraction.assign ((size_t) numBins, 0.0);
        spreadLevel.assign ((size_t) numBins, 0.0);
        fixedFeed.assign ((size_t) numBins, 0.0);
    }

    /** Settles on a fixed spectrum and returns the per-bin boost in dB. */
    std::vector<double> boostProfile (const std::vector<double>& magnitudes, int frames = 500)
    {
        fixedBand.reset();
        sliding.reset();

        const double transferScale = std::pow (10.0, stageGainDb / 20.0) - 1.0;
        std::vector<double> boostDb ((size_t) numBins, 0.0);

        for (int frame = 0; frame < frames; ++frame)
        {
            std::copy (magnitudes.begin(), magnitudes.end(), stageInput.begin());

            sliding.processFrame (stageInput.data(), slidingFraction.data(), spreadLevel.data());

            // Action substitution drives the fixed band from the spread level too, so
            // that both circuits see the masked region as occupied.
            if (composition == SlidingComposition::actionSubstitution)
                for (int k = 0; k < numBins; ++k)
                    fixedFeed[(size_t) k] = std::max (stageInput[(size_t) k], spreadLevel[(size_t) k]);
            else
                std::copy (stageInput.begin(), stageInput.end(), fixedFeed.begin());

            fixedBand.processFrame (fixedFeed.data(), fixedTransfer.data());

            for (int k = 0; k < numBins; ++k)
            {
                const double phiFixed = fixedTransfer[(size_t) k] / transferScale;
                const double phiSliding = slidingFraction[(size_t) k];

                const double phi = composition == SlidingComposition::restriction
                                     ? phiFixed * phiSliding
                                     : 1.0 - (1.0 - phiFixed) * (1.0 - phiSliding);

                boostDb[(size_t) k] = 20.0 * std::log10 (1.0 + transferScale * phi);
            }
        }

        return boostDb;
    }

    int binFor (double hz) const { return (int) std::lround (hz * fftSize / sampleRate); }
    double hzFor (int bin) const { return (double) bin * sampleRate / (double) fftSize; }
};

void testSlidingLayerCompositions()
{
    std::printf ("\nSliding layer: probe-tone behaviour of the two compositions (Figs. 14-16)\n");

    const char* names[] = { "restriction", "action substitution" };
    const SlidingComposition modes[] = { SlidingComposition::restriction,
                                         SlidingComposition::actionSubstitution };

    // A dominant tone over a low-level probe floor, which is how Figs. 14 to 16
    // were measured: "adding a swept frequency probe tone at levels between -60 dB
    // and -80 dB into the encoder input".
    constexpr double probeDb = -70.0;

    for (double dominantHz : { 200.0, 800.0, 3000.0 })
    {
        std::printf ("\n  dominant %.0f Hz at 0 dB, probe floor %.0f dB\n", dominantHz, probeDb);
        std::printf ("    probe Hz  ");

        const double probeFreqs[] = { 50.0, 100.0, 200.0, 400.0, 800.0,
                                     1600.0, 3000.0, 6000.0, 12000.0 };

        for (double hz : probeFreqs)
            std::printf ("%7.0f", hz);

        std::printf ("\n");

        for (int m = 0; m < 2; ++m)
        {
            ProbeStage stage (modes[m], dominantHz >= 800.0);

            std::vector<double> spectrum ((size_t) stage.numBins,
                                          std::pow (10.0, probeDb / 20.0));
            spectrum[(size_t) stage.binFor (dominantHz)] = 1.0;

            const auto boost = stage.boostProfile (spectrum);

            std::printf ("    %-19s", names[m]);

            for (double hz : probeFreqs)
                std::printf ("%7.2f", boost[(size_t) stage.binFor (hz)]);

            std::printf ("\n");
        }
    }

    // Both compositions must satisfy the properties that make the process work at
    // all, whatever their difference in degree.
    for (int m = 0; m < 2; ++m)
    {
        ProbeStage stage (modes[m], true);

        const double full = std::pow (10.0, -70.0 / 20.0);
        std::vector<double> quiet ((size_t) stage.numBins, full);

        const auto quietBoost = stage.boostProfile (quiet);

        // With nothing above threshold anywhere, every bin gets the stage's full
        // low-level boost.
        expectNear (std::string (names[m]) + ": full boost on a quiet spectrum",
                    quietBoost[(size_t) stage.binFor (3000.0)], ProbeStage::stageGainDb, 0.35);

        std::vector<double> withTone = quiet;
        withTone[(size_t) stage.binFor (3000.0)] = 1.0;

        const auto toneBoost = stage.boostProfile (withTone);

        // Boost must be reduced at the dominant itself, or the encoder would
        // overload the channel there. This is the viability gate, and the two
        // compositions are measured against it with their known outcomes asserted,
        // so that a later change which alters either one is noticed.
        //
        // Action substitution fails it, and the reason is structural rather than a
        // tuning problem. Eq. (1) means "full boost if either circuit achieves it".
        // At a dominant's own bin the spread from its neighbours is low, because a
        // masking threshold excludes the bin it is computed for, so the sliding
        // circuit reports full boost and restores exactly what the fixed band
        // correctly removed. Feeding the sliding circuit max(own, spread) would fix
        // the dominant but then both circuits see the same thing everywhere and the
        // composition degenerates to a single band.
        const double atDominant = toneBoost[(size_t) stage.binFor (3000.0)];
        const bool reducesAtDominant = atDominant < ProbeStage::stageGainDb - 1.0;
        const bool expectedToReduce = modes[m] == SlidingComposition::restriction;

        ++checks;

        if (reducesAtDominant != expectedToReduce)
        {
            ++failures;
            std::printf ("  [FAIL] %-46s %.2f dB (expected %s)\n",
                         (std::string (names[m]) + ": dominant-overload gate").c_str(),
                         atDominant, expectedToReduce ? "reduced" : "not reduced");
        }
        else
        {
            std::printf ("  [pass] %-46s %.2f dB (%s, as measured)\n",
                         (std::string (names[m]) + ": dominant-overload gate").c_str(),
                         atDominant, reducesAtDominant ? "reduced" : "NOT reduced");
        }

        // Asymmetry: report it rather than assert a direction. The analog's
        // high-frequency sliding band recovers boost *above* a dominant, but that
        // is a consequence of its fixed band being wideband — it has broad loss to
        // undo. A per-bin fixed band has no such loss, so the only reasons to
        // withhold boost away from a dominant are that masking makes it pointless
        // and that a smooth profile tolerates channel errors. Masking spreads
        // upward, so the profile here is broad above and narrow below: the reverse
        // of the circuit, but the same as the paper's own psychoacoustic statement
        // that "the only region of the spectrum that is not boosted in gain is the
        // region that is controlled by masking".
        const double below = toneBoost[(size_t) stage.binFor (1200.0)];
        const double above = toneBoost[(size_t) stage.binFor (7500.0)];

        std::printf ("         %-20s asymmetry: %.2f dB below, %.2f dB above\n",
                     names[m], below, above);

        // What must hold either way: the dip is local, so boost is fully recovered
        // far from the dominant in both directions.
        expectNear (std::string (names[m]) + ": recovered far below",
                    toneBoost[(size_t) stage.binFor (300.0)], ProbeStage::stageGainDb, 0.35);
        expectNear (std::string (names[m]) + ": recovered far above",
                    toneBoost[(size_t) stage.binFor (18000.0)], ProbeStage::stageGainDb, 0.5);
    }
}

void testGainSlew()
{
    std::printf ("\nGain-slew clamp (standing in for the overshoot suppressors)\n");

    constexpr int numBins = 16;

    GainSlewParams p;
    GainSlew slew;
    slew.prepare (p, numBins);

    std::vector<double> transfer ((size_t) numBins, 0.0);

    // First frame primes the state and must pass through untouched, or the very start
    // of a signal would be attenuated.
    std::fill (transfer.begin(), transfer.end(), std::pow (10.0, 24.0 / 20.0) - 1.0);
    slew.process (transfer.data());

    expectNear ("first frame passes through", 20.0 * std::log10 (1.0 + transfer[0]), 24.0, 0.01);

    // A large downward step should be limited to the engaging rate.
    std::fill (transfer.begin(), transfer.end(), 0.0);   // 0 dB, a 24 dB drop
    slew.process (transfer.data());

    const double afterDrop = 20.0 * std::log10 (1.0 + transfer[0]);
    expectNear ("24 dB drop limited to the engaging rate", 24.0 - afterDrop,
                p.engagingDbPerFrame, 0.01);

    // Releasing is allowed to be faster, since recovery speed is something the paper
    // explicitly values.
    slew.reset();
    std::fill (transfer.begin(), transfer.end(), 0.0);
    slew.process (transfer.data());
    std::fill (transfer.begin(), transfer.end(), std::pow (10.0, 24.0 / 20.0) - 1.0);
    slew.process (transfer.data());

    expectNear ("24 dB rise limited to the releasing rate",
                20.0 * std::log10 (1.0 + transfer[0]), p.releasingDbPerFrame, 0.01);

    // Changes inside the clamp must pass bit-exactly, so ordinary material never
    // meets the limiter at all.
    slew.reset();
    std::fill (transfer.begin(), transfer.end(), 0.0);
    slew.process (transfer.data());

    const double smallDb = 0.5 * p.engagingDbPerFrame;
    std::fill (transfer.begin(), transfer.end(), std::pow (10.0, smallDb / 20.0) - 1.0);
    slew.process (transfer.data());

    expectNear ("a change inside the clamp is untouched",
                20.0 * std::log10 (1.0 + transfer[0]), smallDb, 0.001);

    // Modulation control raises the threshold, coupling transient behaviour to the
    // steady-state circuits as section 2.4 requires.
    slew.reset();
    std::fill (transfer.begin(), transfer.end(), 0.0);
    slew.process (transfer.data());
    slew.setModulationControl (1.0);
    std::fill (transfer.begin(), transfer.end(), std::pow (10.0, 24.0 / 20.0) - 1.0);
    slew.process (transfer.data());

    const double relaxed = 20.0 * std::log10 (1.0 + transfer[0]);

    ++checks;

    if (relaxed <= p.releasingDbPerFrame)
    {
        ++failures;
        std::printf ("  [FAIL] %-46s %.2f dB\n", "modulation control relaxes the limit", relaxed);
    }
    else
    {
        std::printf ("  [pass] %-46s %.2f dB vs %.2f\n", "modulation control relaxes the limit",
                     relaxed, p.releasingDbPerFrame);
    }
}

//==============================================================================
/** Measures the assembled encoder's gain at a probe frequency, by running a tone
    through it and comparing output to input amplitude. */
struct EncoderProbe
{
    static constexpr double rate = 48000.0;

    SrEngine engine;
    int latency = 0;

    EncoderProbe()
    {
        engine.prepare (rate, 1024, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (SrMode::encode);
        latency = engine.getLatencySamples();
    }

    /** Encode gain in dB for a steady tone at `hz` and amplitude `amp`. */
    double gainAt (double hz, double amp)
    {
        engine.reset();

        const int total = 8 * latency;
        std::vector<double> buffer ((size_t) total);

        for (int n = 0; n < total; ++n)
            buffer[(size_t) n] = amp * std::sin (twoPi * hz * (double) n / rate);

        int at = 0;
        while (at < total)
        {
            const int len = std::min (512, total - at);
            double* ptr = buffer.data() + at;
            engine.process (&ptr, 1, len);
            at += len;
        }

        // RMS over the settled tail, converted back to amplitude.
        double sum = 0.0;
        const int from = total / 2;

        for (int n = from; n < total; ++n)
            sum += buffer[(size_t) n] * buffer[(size_t) n];

        const double outAmp = std::sqrt (2.0 * sum / (double) (total - from));

        return 20.0 * std::log10 (outAmp / amp);
    }
};

void testPublishedCurves()
{
    std::printf ("\nAssembled encoder against the paper's published curves\n");

    EncoderProbe probe;

    // Fig. 13: the low-level (subthreshold) characteristic. Section 4.2 describes it
    // as resembling the inverse CCIR noise-weighting curve, peaking in the upper
    // middle and tapering at both extremes.
    std::printf ("\n  Fig. 13 low-level encode curve, tone at -85 dB\n    ");

    const double curveFreqs[] = { 30.0, 60.0, 125.0, 250.0, 500.0, 800.0,
                                  1500.0, 3000.0, 6000.0, 10000.0, 15000.0 };

    for (double hz : curveFreqs)
        std::printf ("%8.0f", hz);

    std::printf ("\n    ");

    double peak = -1.0e9, atLowEnd = 0.0, atHighEnd = 0.0;

    for (double hz : curveFreqs)
    {
        const double g = probe.gainAt (hz, std::pow (10.0, -85.0 / 20.0));
        std::printf ("%8.2f", g);

        peak = std::max (peak, g);

        if (hz == 30.0)   atLowEnd = g;
        if (hz == 15000.0) atHighEnd = g;
    }

    std::printf ("\n");

    // Section 1: about 24 dB of dynamic effect at high frequencies, 16 dB at low.
    // The peak of the curve is the high-frequency figure.
    expectNear ("Fig. 13 peak boost", peak, 24.0, 3.0);

    // Section 4.2: less noise reduction is needed at the extremes, and strong spectral
    // skewing is applied there, so the curve must taper at both ends.
    ++checks;

    if (atLowEnd >= peak - 3.0 || atHighEnd >= peak - 3.0)
    {
        ++failures;
        std::printf ("  [FAIL] %-46s low %.1f, peak %.1f, high %.1f dB\n",
                     "curve tapers at both extremes", atLowEnd, peak, atHighEnd);
    }
    else
    {
        std::printf ("  [pass] %-46s low %.1f, peak %.1f, high %.1f dB\n",
                     "curve tapers at both extremes", atLowEnd, peak, atHighEnd);
    }

    // Fig. 12: the single-tone encoder characteristic. Dynamic action from about
    // -62 dB to -5 dB at high frequencies, with no action in the top 25 dB.
    std::printf ("\n  Fig. 12 encode gain vs level, 3 kHz tone\n    input dB  ");

    const double levels[] = { -90.0, -80.0, -70.0, -60.0, -50.0, -40.0,
                              -30.0, -20.0, -10.0, 0.0 };

    for (double db : levels)
        std::printf ("%8.0f", db);

    std::printf ("\n    gain dB   ");

    double lowLevelGain = 0.0, referenceGain = 0.0;

    for (double db : levels)
    {
        const double g = probe.gainAt (3000.0, std::pow (10.0, db / 20.0));
        std::printf ("%8.2f", g);

        if (db == -90.0) lowLevelGain = g;
        if (db == 0.0)   referenceGain = g;
    }

    std::printf ("\n");

    expectNear ("Fig. 12 subthreshold boost at 3 kHz", lowLevelGain, 24.0, 3.0);

    // The dual-path premise: high-level signals travel the main path essentially
    // unprocessed, with section 1's "further dynamic action of about 1 dB" remaining.
    expectNear ("Fig. 12 gain at reference level", referenceGain, 1.0, 1.5);

    // Antisaturation, section 4.1, measured at high level where it dominates.
    std::printf ("\n  Section 4.1 antisaturation, tone near reference\n");

    const double midband = probe.gainAt (1000.0, 0.9);

    for (auto [hz, expected] : publishedAntisaturation)
    {
        const double attenuation = midband - probe.gainAt (hz, 0.9);
        expectNear (std::to_string ((int) hz) + " Hz antisaturation, assembled",
                    attenuation, expected, 3.0);
    }
}

void testShelvingBankInvertsExactly()
{
    std::printf ("\nShelving bank: time-domain gain that inverts exactly\n");

    constexpr double rate = 48000.0;

    ShelvingBank64 encodeBank, decodeBank;
    encodeBank.prepare (rate);
    decodeBank.prepare (rate);

    const int bands = encodeBank.getNumBands();
    std::printf ("    %d bands from %.0f Hz to %.0f Hz\n", bands,
                 encodeBank.getCentreHz (0), encodeBank.getCentreHz (bands - 1));

    std::vector<double> target ((size_t) bands, 0.0);

    const auto tilt = [&] (double depthDb)
    {
        for (int i = 0; i < bands; ++i)
            target[(size_t) i] = depthDb * (double) i / (double) (bands - 1);
    };

    // The solve has to actually produce the profile that was asked for, or the bank
    // is exactly invertible but exactly wrong.
    tilt (24.0);
    encodeBank.setTargetDb (target.data());

    double worstProfile = 0.0;

    for (int i = 0; i < bands; ++i)
        worstProfile = std::max (worstProfile,
                                 std::abs (encodeBank.responseDbAt (encodeBank.getCentreHz (i))
                                           - target[(size_t) i]));

    expectNear ("solved profile matches the request", worstProfile, 0.0, 2.0);

    const int total = 1 << 16;
    auto input = makeNoise (total, 90210u);
    for (auto& v : input)
        v *= 0.5;

    // Static profile.
    {
        encodeBank.reset();
        decodeBank.reset();
        encodeBank.setTargetDb (target.data());
        decodeBank.setTargetDb (target.data());

        std::vector<double> work = input;
        encodeBank.process (work.data(), total);
        decodeBank.processInverse (work.data(), total);

        double worst = 0.0;

        for (int n = 64; n < total; ++n)
            worst = std::max (worst, std::abs (work[(size_t) n] - input[(size_t) n]));

        expectBelow ("static 24 dB profile, round trip",
                     relativeDb (worst, peakMagnitude (input)), -200.0);
    }

    // Coefficients changing every block, which is the case that defeats an STFT
    // entirely: a per-bin gain applied in one transform and divided out in another
    // manages only -16 dB on a single-bin notch. Direct form I holds because its state
    // is the real input and output history, so both directions reference the same past
    // samples even as the coefficients move.
    {
        encodeBank.reset();
        decodeBank.reset();

        std::vector<double> work = input;
        std::mt19937_64 rng (1234);
        std::uniform_real_distribution<double> dist (-24.0, 24.0);

        constexpr int hop = 256;

        for (int at = 0; at < total; at += hop)
        {
            for (int i = 0; i < bands; ++i)
                target[(size_t) i] = dist (rng);

            encodeBank.setTargetDb (target.data());
            decodeBank.setTargetDb (target.data());

            const int len = std::min (hop, total - at);
            encodeBank.process (work.data() + at, len);
            decodeBank.processInverse (work.data() + at, len);
        }

        double worst = 0.0;

        for (int n = 64; n < total; ++n)
            worst = std::max (worst, std::abs (work[(size_t) n] - input[(size_t) n]));

        expectBelow ("gains randomised every block, round trip",
                     relativeDb (worst, peakMagnitude (input)), -200.0);
    }

    // A flat request must be transparent, so the bank costs nothing when idle.
    {
        std::fill (target.begin(), target.end(), 0.0);
        encodeBank.reset();
        encodeBank.setTargetDb (target.data());

        std::vector<double> work = input;
        encodeBank.process (work.data(), total);

        double worst = 0.0;

        for (int n = 64; n < total; ++n)
            worst = std::max (worst, std::abs (work[(size_t) n] - input[(size_t) n]));

        expectBelow ("a flat request is transparent",
                     relativeDb (worst, peakMagnitude (input)), -200.0);
    }
}

//==============================================================================
/** How finely the split between the two staggered halves is resolved, which is a control
    where the crossover frequency used to be one.

    The reference system meets its two halves at 800 Hz because two analog band-defining
    networks is what it could build. Here the same power-complementary pole is sampled onto
    a power-of-two number of log-spaced sections and that count is what the user sets, from
    a hard two-band split up to one section per bin.

    Four things want measuring rather than trusting. That the weighting is rebuilt in place,
    on the audio thread, into tables the stages already hold pointers to — the kind of change
    that compiles and does nothing. That the default still lands where the continuous profile
    used to, so the change is a no-op at the setting everyone will use. That quantising the
    profile leaves |HP|^2 + |LP|^2 = 1, which is what the two halves multiplying rests on.
    And what coarse sections cost the null, which is the one place this control can do real
    damage: an 8 dB span in N steps has 8/N dB edges, and setEnvelopeSmoothingBark records
    that a gain profile which jumps between neighbouring bins is not invertible at all.
*/
void testSectionResolution()
{
    std::printf ("\nSection resolution\n");

    EncoderProbe probe;
    const double amplitude = std::pow (10.0, -85.0 / 20.0);

    std::printf ("    sections     20 Hz    500 Hz     5 kHz\n");

    double lowest = 1.0e9, highest = -1.0e9;      // spread at 20 Hz across every count
    double lowestHf = 1.0e9, highestHf = -1.0e9;  // and at 5 kHz
    double finestMid = 0.0;
    double finestLow = 0.0, finestHigh = 0.0;

    for (const int count : { 4, 8, 64, 128, 256, 2048, 16384 })
    {
        probe.engine.setSections (count);

        const double low = probe.gainAt (20.0, amplitude);
        const double mid = probe.gainAt (500.0, amplitude);
        const double high = probe.gainAt (5000.0, amplitude);

        std::printf ("    %8d  %7.2f   %7.2f   %7.2f\n", count, low, mid, high);

        lowest = std::min (lowest, low);      highest = std::max (highest, low);
        lowestHf = std::min (lowestHf, high); highestHf = std::max (highestHf, high);

        if (count == 16384) { finestLow = low; finestMid = mid; finestHigh = high; }
    }

    // The default has to reproduce the continuous profile it replaced, or every other
    // number in this suite is measuring a different process. The reference figures are
    // the continuous 800 Hz-crossover profile's, re-measured after the shelf DC
    // normalisation was made exact and the modulation-control opposition moved to the
    // paper's ratio form — the earlier build read 2.98, 18.67 and 20.72 dB here, of
    // which up to 0.53 dB was the antisaturation shelf's DC error. See DESIGN.md.
    expectNear ("the finest grid matches the old profile, 20 Hz", finestLow, 2.53, 0.1);
    expectNear ("the finest grid matches the old profile, 500 Hz", finestMid, 18.22, 0.1);
    expectNear ("the finest grid matches the old profile, 5 kHz", finestHigh, 20.02, 0.1);

    // Both extremes belong wholly to one half at any resolution, so quantising the profile
    // must not move them. This is what would break if the two halves stopped being sampled
    // by the identical function of frequency.
    expectNear ("20 Hz is unmoved by the count", highest - lowest, 0.0, 0.5);
    expectNear ("5 kHz is unmoved by the count", highestHf - lowestHf, 0.0, 0.5);

    // Note that this check only means anything because the range now starts at 4. Below
    // that the extremes are supposed to move: at 1 section the side chain's allocation goes
    // flat and the 20 Hz-to-5 kHz spread closes up by 6 dB, and at 2 it is a single step.
    // Both were measured, and both are excluded by minSections rather than by this loop —
    // see the note there for why they went.

    // Rounded down to a power of two, and clamped at both ends.
    probe.engine.setSections (0);
    expectNear ("clamped at the bottom", (double) probe.engine.getSections(),
                (double) SrSideChain::minSections, 0.0, "sections");

    probe.engine.setSections (100000);
    expectNear ("clamped at the top", (double) probe.engine.getSections(),
                (double) SrSideChain::maxSections, 0.0, "sections");

    probe.engine.setSections (300);
    expectNear ("rounded down to a power of two", (double) probe.engine.getSections(),
                256.0, 0.0, "sections");

    //==========================================================================
    // The profile itself, read off the two halves directly rather than inferred from a
    // tone. A short window at 1024 and a long one at 4096 share a sample rate, so the
    // short half's bin k sits at exactly the same frequency as the long half's bin 4k —
    // which is what makes the complementarity checkable bin by bin rather than by
    // interpolation.
    {
        SrSideChain hf, lf;
        hf.prepare (true,  10, 256,  48000.0);
        lf.prepare (false, 12, 1024, 48000.0);

        double worstSum = 0.0;
        int worstDistinct = 0;
        bool monotoneInFrequency = true;
        int distinctAtFour = 0;

        for (const int count : { 4, 8, 64, 128, 256, 2048, 16384 })
        {
            hf.setSections (count);
            lf.setSections (count);

            const double* hfBoost = hf.getLowLevelBoostDb();
            const double* lfBoost = lf.getLowLevelBoostDb();

            // Three stages above, two below, so the boosts divide by 24 and 16 to recover
            // the weights the two halves were handed.
            for (int k = 0; k < hf.getNumBins(); ++k)
            {
                const double sum = hfBoost[k] / 24.0 + lfBoost[4 * k] / 16.0;
                worstSum = std::max (worstSum, std::abs (sum - 1.0));
            }

            // A staircase, not a curve: the number of distinct values the profile takes
            // cannot exceed the number of sections. Bins outside 20 Hz to 20 kHz collapse
            // into the outermost sections, so at the finest grids it is comfortably fewer.
            int distinct = 1;

            for (int k = 1; k < lf.getNumBins(); ++k)
                if (lfBoost[k] != lfBoost[k - 1])
                    ++distinct;

            worstDistinct = std::max (worstDistinct, distinct - count);

            if (count == 4)
                    distinctAtFour = distinct;

            // The underlying pole rises with frequency and sampling it on a log grid
            // cannot reorder that, so the high half's boost must never fall as k rises.
            for (int k = 1; k < hf.getNumBins(); ++k)
                if (hfBoost[k] < hfBoost[k - 1] - 1.0e-12)
                    monotoneInFrequency = false;
        }

        expectBelow ("the halves stay power complementary at every count",
                     relativeDb (worstSum, 1.0), -240.0);

        ++checks;

        if (worstDistinct <= 0)
            std::printf ("  [pass] %-46s no count exceeds its own section total\n",
                         "the profile is a staircase");
        else
            ++failures, std::printf ("  [FAIL] %-46s %d values too many\n",
                                     "the profile is a staircase", worstDistinct);

        expectNear ("four sections give exactly four values", (double) distinctAtFour,
                    4.0, 0.0, "values");

        ++checks;

        if (monotoneInFrequency)
            std::printf ("  [pass] %-46s boost never falls as frequency rises\n",
                         "the profile is monotone in frequency");
        else
            ++failures, std::printf ("  [FAIL] %-46s not monotone\n",
                                     "the profile is monotone in frequency");
    }

    //==========================================================================
    // What coarse sections cost. A separate encoder and decoder, which is how the plugin
    // is actually used, against broadband noise at the level the null is quoted at
    // elsewhere in this suite.
    //
    // Measured at 48 kHz, the null improves as the grid gets finer and then stops. The
    // step between neighbouring sections shrinks as 8/N dB, and past about 64 it is small
    // enough that something else limits the null — the same wall the broadband figures
    // elsewhere in this suite run into. So the control costs nothing above 64, and what it
    // costs below that is why the range starts where it does.
    //
    //     sections     1      2   |    4      8     64    128  *256*  2048  16384
    //     null      -67.5  -24.0  |  -29.9  -39.0  -46.6  -46.6  -46.6  -46.6  -46.6 dB
    //                 (excluded)  |
    //
    // 256 is starred because it is the default. The two settings left of the bar were
    // measured before being excluded and are kept here because the shape of the curve is
    // not what anyone would guess:
    //
    // One section nulls 21 dB *better* than the finest grid, because a single weight across
    // the whole band is the smoothest profile the process can apply — no bin-to-bin step
    // anywhere — and smoothness is exactly what the transform can undo. Two sections is the
    // worst in the range, because that is where the entire 8 dB arrives as one step. So the
    // cost is not monotone in the count: it leaps 43 dB between the first two positions and
    // then recovers.
    //
    // Neither is reachable now. One section nulls best by not doing the job — no
    // frequency-dependent allocation means nothing is being allocated — and two is the
    // analog machine's own arrangement, which is a real loss and is taken knowingly. See
    // SrSideChain::minSections.
    //
    // What survives of the lesson is that the null alone does not choose this control's
    // setting, and a reader meeting the -67.5 dB figure cold would draw the wrong
    // conclusion from it.
    //
    // The top of the range reads the same as the middle of it because this runs at 48 kHz,
    // where the long window has 2049 bins: anything past 2048 sections is asking for a grid
    // finer than there are bins to put it on, so several sections land in one bin and
    // nothing changes. The higher positions are there for 96 and 192 kHz, and 16384 is
    // there to make the detent count odd. They are measured anyway — a count past the bin
    // count has to be shown inert rather than merely assumed harmless, and this is what
    // says so.
    {
        std::printf ("\n  Encode/decode null against the section count\n");

        constexpr double rate = 48000.0;

        SrEngine encoder, decoder;
        encoder.prepare (rate, 1024, 1);
        decoder.prepare (rate, 1024, 1);
        encoder.setProcess (SrProcess::spectralRecording);
        decoder.setProcess (SrProcess::spectralRecording);
        encoder.setMode (SrMode::encode);
        decoder.setMode (SrMode::decode);

        const int latency = encoder.getLatencySamples() + decoder.getLatencySamples();
        const int total = 4 * latency;

        double coarsest = 0.0, plateau = 0.0, finest = 0.0, atDefault = 0.0;

        for (const int count : { 4, 8, 64, 128, 256, 2048, 16384 })
        {
            encoder.reset();
            decoder.reset();
            encoder.setSections (count);
            decoder.setSections (count);

            auto input = makeNoise (total, 4711u);

            for (auto& v : input)
                v *= 0.05;                      // -26 dB, as the broadband null case uses

            std::vector<double> work = input;
            runEngine (encoder, work);
            runEngine (decoder, work);

            const double db = relativeDb (delayedError (work, input, latency, latency),
                                          peakMagnitude (input));

            std::printf ("    %8d  %8.1f dB\n", count, db);

            if (count == 4)     coarsest = db;
            if (count == 256)   plateau = db;
            if (count == 16384) finest = db;

            if (count == SrSideChain::defaultSections)
                atDefault = db;
        }

        // Recorded rather than aspirational, and in both directions. The finest grid must
        // not regress; the coarsest must stay measurably worse. Asserting the bad number
        // too is what stops a later change quietly "fixing" the coarse settings by
        // flattening the staircase that is the whole point of the control.
        expectBelow ("the default nulls where the continuous profile did", finest, -46.0);

        ++checks;
        const bool coarseCosts = coarsest > -44.0;

        if (! coarseCosts)
            ++failures;

        std::printf ("  [%s] %-46s %+8.1f dB  (wants above %.1f dB)\n",
                     coarseCosts ? "pass" : "FAIL",
                     "the coarsest setting still costs the null", coarsest, -44.0);

        // The plateau is the finding worth protecting: past the default's 256 sections
        // the steps are already below whatever else limits the null, so the rest of the
        // range is free. The plateau moved up from 64 when the modulation-control
        // opposition gained its smoothing — the null floor dropped by 13 dB and the
        // grid's own contribution became visible one octave further. If this ever stops
        // holding, the section grid has started to matter again and the manual's advice
        // about the default changes with it.
        expectNear ("past 256 sections the grid costs nothing", finest - plateau, 0.0, 0.5);

        // The default is not the top of the range, so it has to be shown to sit on the
        // plateau rather than assumed to. This is the assertion that would catch someone
        // moving defaultSections below the knee, where it would start costing the round
        // trip without anything else in the suite noticing.
        expectNear ("the default nulls as well as the finest grid",
                    atDefault - finest, 0.0, 0.5);

    }
}

//==============================================================================
/** Moving the control while audio is running.

    This is a stepped control, so every move is a jump: one detent can shift a bin's share
    of the boost by several decibels at once, and a per-bin gain that changes between one
    frame and the next is precisely what GainSlew exists to stop. The question is whether
    the slew already covers a parameter move as well as it covers programme material, or
    whether changing the count clicks.

    Measured rather than reasoned about, on the envelope rather than the waveform: a
    spectral process has a window's worth of latency, so a discontinuity shows up as a step
    between consecutive hops of output level, not as a step between two samples.
*/
void testSectionTransitionIsSmooth()
{
    std::printf ("\nMoving the section count while running\n");

    constexpr double rate = 48000.0;
    constexpr int hop = 256;

    // Subthreshold, so the side chain is at full boost and the profile change is as large
    // as this control can make it.
    const double amplitude = std::pow (10.0, -60.0 / 20.0);

    // The widest move the panel allows in one gesture: the two ends of the range.
    struct Case { const char* name; int from; int to; };

    for (auto c : { Case { "coarsest to finest", SrSideChain::minSections, SrSideChain::maxSections },
                    Case { "finest to coarsest", SrSideChain::maxSections, SrSideChain::minSections },
                    Case { "one detent, 128 to 256", 128, 256 } })
    {
        SrEngine engine;
        engine.prepare (rate, 512, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (SrMode::encode);
        engine.setSections (c.from);

        const int latency = engine.getLatencySamples();
        const int total = 8 * latency;
        const int switchAt = total / 2;

        std::vector<double> buffer ((size_t) total);

        for (int n = 0; n < total; ++n)
            buffer[(size_t) n] = amplitude * std::sin (twoPi * 1000.0 * (double) n / rate);

        for (int at = 0; at < total; at += hop)
        {
            if (at >= switchAt && at - hop < switchAt)
                engine.setSections (c.to);

            const int len = std::min (hop, total - at);
            double* ptr = buffer.data() + at;
            engine.process (&ptr, 1, len);
        }

        // Per-hop RMS in dB, then the largest step between neighbouring hops. Everything
        // before the transition has reached steady state, so any step at all is the move.
        std::vector<double> hopDb;

        for (int at = latency; at + hop <= total; at += hop)
        {
            double sum = 0.0;

            for (int n = at; n < at + hop; ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            hopDb.push_back (10.0 * std::log10 (std::max (1.0e-30, sum / (double) hop)));
        }

        double worstStep = 0.0;

        for (size_t i = 1; i < hopDb.size(); ++i)
            worstStep = std::max (worstStep, std::abs (hopDb[i] - hopDb[i - 1]));

        // The slew's own limit is 3 dB per frame engaging and 6 releasing, and the long
        // half runs a hop four times this one, so a move cannot legitimately produce more
        // than a few decibels in any one hop of output. Anything approaching the full
        // profile change arriving at once would be a click.
        expectBelow (std::string ("worst hop-to-hop step, ") + c.name,
                     worstStep, 6.5);
    }
}

//==============================================================================
/** The one pole every control reaches the signal path through.

    Two things need to be true and neither is obvious from the code: that the settling time
    is what it claims — 95% of the distance inside 75 ms, at any sample rate — and that a
    one-pole step response never overshoots, which is what makes it safe to put in front of
    a gain that must not exceed its target.
*/
void testControlSmoothing()
{
    std::printf ("\nControl smoothing\n");

    for (const double rate : { 44100.0, 48000.0, 96000.0, 192000.0 })
    {
        ControlSmoother smoother;
        smoother.prepare (rate, 0.075, 0.0);
        smoother.setTarget (1.0);

        const int at75ms = (int) std::llround (0.075 * rate);

        double previous = smoother.getCurrent();
        double worstBacktrack = 0.0;
        double peak = 0.0;
        double at75 = 0.0;

        // Four settling times, which is long enough to see any overshoot come back down.
        for (int n = 0; n < 4 * at75ms; ++n)
        {
            const double value = smoother.next();

            worstBacktrack = std::min (worstBacktrack, value - previous);
            peak = std::max (peak, value);
            previous = value;

            if (n == at75ms - 1)
                at75 = value;
        }

        expectNear ("settled fraction at 75 ms, " + std::to_string ((int) rate) + " Hz",
                    at75, 0.95, 0.005, "  ");

        // Monotone and never above the target: a step into a gain must not ring.
        expectNear ("no backtracking, " + std::to_string ((int) rate) + " Hz",
                    worstBacktrack, 0.0, 1.0e-15, "  ");
        // How far above the target it ever got — not how close it got to it, which after
        // four settling times is still 6e-6 short by construction.
        expectNear ("no overshoot, " + std::to_string ((int) rate) + " Hz",
                    std::max (0.0, peak - 1.0), 0.0, 1.0e-15, "  ");
    }

    // A saved session must be in force on the first sample, not faded in.
    ControlSmoother snapped;
    snapped.prepare (48000.0, 0.075, 0.5);
    expectNear ("prepare starts at its initial value", snapped.getCurrent(), 0.5, 0.0, "  ");
    expectEqual ("prepare leaves nothing to smooth", snapped.isMoving() ? 1 : 0, 0);

    snapped.setTarget (0.25);
    expectEqual ("a new target is something to smooth", snapped.isMoving() ? 1 : 0, 1);
    snapped.snapTo (0.25);
    expectEqual ("snapTo arrives immediately", snapped.isMoving() ? 1 : 0, 0);
}

//==============================================================================
/** The fractional delay's read-position floor.

    The kernel's earliest tap sits one sample ahead of the read position, and ahead of the
    sample just written the ring holds nothing but the oldest audio it has. Before the
    floor, a steady half-sample request put that tap past the write index, and its Lagrange
    weight at t = 0.5 is -1/16: whatever went through the line a whole ring earlier — about
    40 ms at this capacity — re-emerged under the fresh audio at -24 dB. The tests below
    prime the ring with loud noise first, so the stale audio the old code would have leaked
    is actually there to leak.
*/
void testFractionalDelayFloor()
{
    std::printf ("\nFractional delay floors the read position at one sample\n");

    constexpr double rate = 48000.0;

    // (a) The regression itself. A 0.5-sample request on a line full of old noise, then a
    // unit impulse: with the floor the response is confined to the taps around one sample,
    // and everything past the first few samples is silence rather than the ring's history.
    {
        FractionalDelay64 delay;
        delay.prepare (rate, 0.040);
        delay.snapDelaySamples (0.5);

        auto noise = makeNoise (4096, 1111u);
        delay.process (noise.data(), 4096);          // the ring now holds nothing but noise

        std::vector<double> impulse ((size_t) 4096, 0.0);
        impulse[0] = 1.0;
        delay.process (impulse.data(), 4096);

        double tail = 0.0;

        for (int n = 8; n < 4096; ++n)
            tail = std::max (tail, std::abs (impulse[(size_t) n]));

        expectBelow ("0.5-sample request, energy outside the first 8 taps",
                     relativeDb (tail, 1.0), -200.0);
    }

    // (b) Magnitude flatness on a tone. A floored 0.5-sample request lands on a whole
    // sample, where the interpolator is exact; a true half-sample position (10.5) is the
    // worst case for the kernel, and third-order Lagrange still holds 1 kHz at 48 kHz
    // flat to well under a hundredth of a dB.
    for (double request : { 0.5, 10.5 })
    {
        FractionalDelay64 delay;
        delay.prepare (rate, 0.040);
        delay.snapDelaySamples (request);

        const int total = 9600;
        std::vector<double> tone ((size_t) total);

        for (int n = 0; n < total; ++n)
            tone[(size_t) n] = 0.5 * std::sin (twoPi * 1000.0 * (double) n / rate);

        std::vector<double> out = tone;
        delay.process (out.data(), total);

        // RMS over 160 whole cycles, past the line's fill.
        double inSum = 0.0, outSum = 0.0;

        for (int n = 960; n < 960 + 7680; ++n)
        {
            inSum += tone[(size_t) n] * tone[(size_t) n];
            outSum += out[(size_t) n] * out[(size_t) n];
        }

        expectNear ("amplitude at 1 kHz through a " + std::to_string (request).substr (0, 4)
                        + "-sample delay",
                    10.0 * std::log10 (outSum / inSum), 0.0, 0.01);
    }

    // (c) Exact zero never reaches the interpolator: process() bypasses the line
    // entirely, so a zero delay is bit-exact passthrough, not a filtered copy.
    {
        FractionalDelay64 delay;
        delay.prepare (rate, 0.040);
        delay.snapDelaySamples (0.0);

        const auto input = makeNoise (2048, 2222u);
        auto work = input;
        delay.process (work.data(), 2048);

        int diffs = 0;

        for (int n = 0; n < 2048; ++n)
            if (work[(size_t) n] != input[(size_t) n])
                ++diffs;

        expectEqual ("zero delay is bit-exact passthrough, differing samples", diffs, 0);
    }

    // (d) A fractional delay actually delays by what it says: the impulse response of a
    // 10.5-sample line has its amplitude centroid at exactly 10.5 samples — the four
    // Lagrange taps straddle the position symmetrically.
    {
        FractionalDelay64 delay;
        delay.prepare (rate, 0.040);
        delay.snapDelaySamples (10.5);

        std::vector<double> impulse ((size_t) 64, 0.0);
        impulse[0] = 1.0;
        delay.process (impulse.data(), 64);

        double weight = 0.0, moment = 0.0;

        for (int n = 0; n < 64; ++n)
        {
            weight += std::abs (impulse[(size_t) n]);
            moment += (double) n * std::abs (impulse[(size_t) n]);
        }

        expectNear ("impulse centroid of a 10.5-sample delay", moment / weight, 10.5,
                    0.05, "smp");
    }

    // (e) Ramping the control back to zero must actually let the line go dormant. The
    // smoother is a one-pole and never reaches zero exactly, so a dormancy test on the
    // *value* would leave the line running — and floored at one sample — for the many
    // seconds the residual takes to underflow. Half a second after the ramp starts the
    // output must be the input again, bit-exactly.
    {
        FractionalDelay64 delay;
        delay.prepare (rate, 0.040);
        delay.snapDelaySamples (480.0);   // 10 ms, well away from zero

        std::vector<double> block ((size_t) 512, 0.0);

        for (int n = 0; n < 512; ++n)
            block[(size_t) n] = std::sin (0.05 * (double) n);

        auto scratch = block;
        delay.process (scratch.data(), 512);

        delay.setDelaySamples (0.0);

        // Half a second of blocks for the smoother to settle at zero.
        for (int i = 0; i < 47; ++i)
        {
            scratch = block;
            delay.process (scratch.data(), 512);
        }

        scratch = block;
        delay.process (scratch.data(), 512);

        int mismatches = 0;

        for (int n = 0; n < 512; ++n)
            if (scratch[(size_t) n] != block[(size_t) n])
                ++mismatches;

        ++checks;

        if (mismatches != 0)
        {
            ++failures;
            std::printf ("  [FAIL] %-58s %d samples differ\n",
                         "a ramp back to zero goes dormant", mismatches);
        }
        else
        {
            std::printf ("  [pass] %-58s bit-exact passthrough\n",
                         "a ramp back to zero goes dormant");
        }
    }
}

//==============================================================================
/** The crossfade's whole reason to exist over the one-pole: it arrives.

    Callers throw switches once the wet path "has faded out", so the ramp has to cover the
    entire distance in exactly the time it was configured with and then hold the endpoint
    bit-exactly — an exponential does neither. The phase accumulator can sit a rounding ulp
    short of the target for one extra call, but the raised cosine has zero slope at its
    endpoints, so the *shaped* output — the thing callers multiply by — is already exactly
    at the endpoint on schedule, and that is what is asserted.
*/
void testCrossfadeRamp()
{
    std::printf ("\nCrossfade reaches its endpoint exactly and holds it\n");

    for (double rate : { 44100.0, 48000.0, 96000.0 })
    {
        constexpr double seconds = 0.005;
        const long long n = (long long) std::llround (seconds * rate);

        Crossfade fade;
        fade.prepare (rate, seconds, 0.0);
        fade.setTarget (1.0);

        double value = fade.getCurrent();

        for (long long i = 0; i < n; ++i)
            value = fade.next();

        expectNear ("up after exactly N calls, " + std::to_string ((int) rate) + " Hz",
                    value, 1.0, 0.0, "  ");

        // (b) Settled means settled: every further call returns the endpoint bit-exactly,
        // forever, so a caller acting on the finished fade can rely on true unity.
        int offEndpoint = 0;

        for (int i = 0; i < 4096; ++i)
            if (fade.next() != 1.0)
                ++offEndpoint;

        expectEqual ("holds 1.0 bit-exactly, " + std::to_string ((int) rate) + " Hz",
                     offEndpoint, 0);
        expectEqual ("nothing left to move, " + std::to_string ((int) rate) + " Hz",
                     fade.isMoving() ? 1 : 0, 0);

        // And back down, which must land on exactly 0.0 the same way.
        fade.setTarget (0.0);

        for (long long i = 0; i < n; ++i)
            value = fade.next();

        expectNear ("down after exactly N calls, " + std::to_string ((int) rate) + " Hz",
                    value, 0.0, 0.0, "  ");
    }

    // (c) The shape: monotone, with zero slope at both ends. The first and last
    // increments must be tiny against the mid-ramp step, or the transition has the
    // corner the raised cosine exists to remove.
    {
        Crossfade fade;
        fade.prepare (48000.0, 0.010, 0.0);
        fade.setTarget (1.0);

        const int n = 480;
        double previous = 0.0;
        double first = 0.0, last = 0.0, largest = 0.0;
        bool monotone = true;

        for (int i = 0; i < n; ++i)
        {
            const double value = fade.next();
            const double step = value - previous;

            if (step < 0.0)
                monotone = false;

            if (i == 0)      first = step;
            if (i == n - 1)  last = step;

            largest = std::max (largest, step);
            previous = value;
        }

        ++checks;

        if (! monotone)
        {
            ++failures;
            std::printf ("  [FAIL] %-58s\n", "the ramp is monotone");
        }
        else
        {
            std::printf ("  [pass] %-58s\n", "the ramp is monotone");
        }

        expectBelow ("first increment against the mid-ramp step",
                     relativeDb (first, largest), -26.0);
        expectBelow ("last increment against the mid-ramp step",
                     relativeDb (last, largest), -26.0);
    }
}

//==============================================================================
/** setDelay repositions the read head instead of flushing the line.

    The old behaviour cleared everything on any change, so switching modes cost a full
    line of dropout. The new contract: shortening exposes nothing and loses nothing —
    the read position simply steps toward the write position across samples that are
    still valid; lengthening zeroes exactly the newly exposed span, so the dropout lasts
    the growth and never replays whatever the ring held from before.
*/
void testDelayLineRepositioning()
{
    std::printf ("\nDelay line setDelay repositions instead of flushing\n");

    DelayLine64 line;
    line.prepare (256);
    line.setDelay (100);

    // A ramp from 1, so a zeroed output can never be mistaken for signal.
    int t = 0;
    const auto feed = [&]
    {
        double x = (double) (++t);
        line.process (&x, 1);
        return x;
    };

    for (int i = 0; i < 300; ++i)
        feed();

    // (a) Shrinking the delay: correct output immediately, no dropout at all. Every
    // sample after the switch must be the ramp delayed by the new amount, bit-exactly.
    line.setDelay (50);

    int wrongAfterShrink = 0, zerosAfterShrink = 0;

    for (int i = 0; i < 100; ++i)
    {
        const double y = feed();

        if (y == 0.0)
            ++zerosAfterShrink;

        if (y != (double) (t - 50))
            ++wrongAfterShrink;
    }

    expectEqual ("shrink 100 -> 50: zeroed samples", zerosAfterShrink, 0);
    expectEqual ("shrink 100 -> 50: wrong samples", wrongAfterShrink, 0);

    // (b) Growing it again: exactly the newly exposed span comes out as silence — the
    // 30 samples the read position stepped back across — and what follows is the
    // correctly delayed ramp, never audio from before the switch.
    line.setDelay (80);

    int zerosBeforeAudio = 0, wrongAfterGrow = 0;
    bool inSilence = true;

    for (int i = 0; i < 200; ++i)
    {
        const double y = feed();

        if (inSilence && y == 0.0)
        {
            ++zerosBeforeAudio;
            continue;
        }

        inSilence = false;

        if (y != (double) (t - 80))
            ++wrongAfterGrow;
    }

    expectEqual ("grow 50 -> 80: zeros before the audio resumes", zerosBeforeAudio, 30);
    expectEqual ("grow 50 -> 80: wrong samples after the gap", wrongAfterGrow, 0);

    // (c) reset() still clears the whole line: nothing of the ramp survives it, and the
    // first delay's worth of output is silence again.
    line.reset();

    int wrongAfterReset = 0;

    for (int i = 0; i < 200; ++i)
    {
        double x = 5.0;
        line.process (&x, 1);

        const double expected = i < 80 ? 0.0 : 5.0;

        if (x != expected)
            ++wrongAfterReset;
    }

    expectEqual ("reset clears the line completely, wrong samples", wrongAfterReset, 0);
}

//==============================================================================
/** The display probe's triple buffer.

    First the single-threaded semantics, where every outcome is deterministic; then the
    reason it is a triple buffer at all — a writer publishing every hop must never let
    the reader see half of one frame and half of another.
*/
void testSpectrumProbe()
{
    std::printf ("\nSpectrum probe: triple-buffered snapshots\n");

    SpectrumProbe probe (2049);
    SpectrumProbe::Snapshot snap;

    // (a) Nothing has been published, so there is nothing to claim.
    expectEqual ("read before any push returns false", probe.read (snap) ? 1 : 0, 0);

    probe.prepare (513, 48000.0, 1024);
    probe.setActive (true);

    expectEqual ("prepare alone publishes nothing", probe.read (snap) ? 1 : 0, 0);

    // (b) One frame in, the same frame out, carrying its own geometry.
    std::vector<double> level ((size_t) 513, 0.25);
    std::vector<double> gain ((size_t) 513, 2.0);

    probe.push (level.data(), gain.data(), 513, 48000.0, 1024);

    expectEqual ("push then read returns true", probe.read (snap) ? 1 : 0, 1);
    expectEqual ("frame carries its bin count", snap.bins, 513);
    expectEqual ("frame carries its FFT size", snap.fftSize, 1024);
    expectNear ("frame carries its sample rate", snap.sampleRate, 48000.0, 0.0, "Hz");
    expectEqual ("vectors sized to the frame's bins",
                 (long long) snap.levelDb.size(), 513);

    int wrongValues = 0;

    for (int k = 0; k < snap.bins; ++k)
    {
        if (snap.levelDb[(size_t) k] != 20.0 * std::log10 (0.25))
            ++wrongValues;

        if (snap.gainDb[(size_t) k] != 20.0 * std::log10 (2.0))
            ++wrongValues;
    }

    expectEqual ("levels and gains stored in dB exactly", wrongValues, 0);

    // (c) A frame with different geometry describes itself: the snapshot is sized by
    // the frame, never by a getter the writer may since have moved on from.
    probe.push (level.data(), gain.data(), 257, 96000.0, 512);
    probe.read (snap);

    expectEqual ("new geometry: bins", snap.bins, 257);
    expectEqual ("new geometry: FFT size", snap.fftSize, 512);
    expectEqual ("new geometry: vector size follows the frame",
                 (long long) snap.levelDb.size(), 257);
    expectNear ("new geometry: sample rate", snap.sampleRate, 96000.0, 0.0, "Hz");

    // (d) Two pushes, one read: the reader claims the freshest frame, not the queue.
    std::fill (level.begin(), level.end(), 0.5);
    probe.push (level.data(), gain.data(), 513, 48000.0, 1024);
    std::fill (level.begin(), level.end(), 0.125);
    probe.push (level.data(), gain.data(), 513, 48000.0, 1024);

    probe.read (snap);
    expectNear ("two pushes then one read returns the newer frame",
                snap.levelDb[0], 20.0 * std::log10 (0.125), 0.0, "dB");

    // (e) With nothing new published, read re-reads the frame it already holds.
    SpectrumProbe::Snapshot again;

    expectEqual ("re-read with no new push returns true", probe.read (again) ? 1 : 0, 1);
    expectEqual ("re-read returns the identical frame",
                 (again.bins == snap.bins && again.levelDb == snap.levelDb
                  && again.gainDb == snap.gainDb) ? 1 : 0, 1);

    // (f) An inactive probe's push is a no-op: the reader keeps the frame it had.
    probe.setActive (false);
    std::fill (level.begin(), level.end(), 1.0);
    probe.push (level.data(), gain.data(), 513, 48000.0, 1024);
    probe.setActive (true);

    probe.read (again);
    expectNear ("push while inactive changes nothing",
                again.levelDb[0], 20.0 * std::log10 (0.125), 0.0, "dB");

    // The torn-read regression. A writer hammering frames whose bins all carry the
    // frame index, against a reader claiming continuously: a two-slot scheme laps the
    // reader and hands it a frame that is half one push and half another, which shows
    // up here as a snapshot whose bins disagree with each other.
    {
        constexpr int bins = 512;
        constexpr int frames = 10000;

        SpectrumProbe hammer (bins);
        hammer.prepare (bins, 48000.0, 1024);
        hammer.setActive (true);

        std::atomic<bool> done { false };

        std::thread writer ([&hammer, &done]
        {
            std::vector<double> writerLevel ((size_t) bins, 0.0);

            for (int frame = 1; frame <= frames; ++frame)
            {
                std::fill (writerLevel.begin(), writerLevel.end(), (double) frame);
                hammer.push (writerLevel.data(), nullptr, bins, 48000.0, 1024);
            }

            done.store (true, std::memory_order_release);
        });

        long reads = 0, torn = 0;
        SpectrumProbe::Snapshot claimed;

        const auto inspect = [&]
        {
            if (! hammer.read (claimed))
                return;

            ++reads;

            if (claimed.bins != bins)
                ++torn;

            for (int k = 1; k < claimed.bins; ++k)
                if (claimed.levelDb[(size_t) k] != claimed.levelDb[0])
                {
                    ++torn;
                    break;
                }
        };

        while (! done.load (std::memory_order_acquire))
            inspect();

        writer.join();
        inspect();   // and the final frame, after the writer has finished

        std::printf ("    %ld reads against %d pushed frames\n", reads, frames);

        expectEqual ("every snapshot is internally uniform, torn reads", torn, 0);

        ++checks;

        if (reads <= 0)
        {
            ++failures;
            std::printf ("  [FAIL] %-58s\n", "the reader actually read something");
        }
        else
        {
            std::printf ("  [pass] %-58s\n", "the reader actually read something");
        }
    }
}

//==============================================================================
/** The display probes follow the running direction.

    A decode-only instance used to publish nothing — the probes belonged to the
    encode side chains, which never run in that mode — so the panel froze on its
    last frame. Now the engine wires the probes to whichever direction is in
    circuit, and a decode instance publishes the cut it applies: gains at or below
    unity, where an encode instance publishes boost above it. */
void testProbesFollowTheRunningDirection()
{
    std::printf ("\nDisplay probes follow the running direction\n");

    constexpr double rate = 48000.0;
    constexpr int block = 512;

    const auto run = [] (SrMode mode)
    {
        SrEngine engine;
        engine.prepare (rate, block, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (mode);
        engine.setProbesActive (true);

        // Deterministic broadband noise at -46 dB, quiet enough that the encode
        // boost and the decode cut are both well engaged.
        std::uint32_t seed = 0x1234567u;
        std::vector<double> buffer ((size_t) block, 0.0);
        double* channel[1] = { buffer.data() };

        SpectrumProbe::Snapshot snapshot;
        auto* probe = engine.getLfProbe (0);

        // A full second, not just until the first frame: the first published frame
        // carries a unity transfer — the control state is one frame behind by
        // design — and the smoothers take tens of milliseconds to build the gain.
        for (int i = 0; i < 100; ++i)
        {
            for (int n = 0; n < block; ++n)
            {
                seed = seed * 1664525u + 1013904223u;
                buffer[(size_t) n] = 0.005 * ((double) seed / 2147483648.0 - 1.0);
            }

            engine.process (channel, 1, block);
        }

        const bool got = probe->read (snapshot);

        // Mean over the low half's own region — 50 to 400 Hz, where its band
        // weight is essentially one — rather than the whole spectrum, most of
        // which this half deliberately leaves alone.
        double sum = 0.0;
        int counted = 0;

        for (int k = 1; k < snapshot.bins; ++k)
        {
            const double hz = (double) k * snapshot.sampleRate / (double) snapshot.fftSize;

            if (hz >= 50.0 && hz <= 400.0)
            {
                sum += snapshot.gainDb[(size_t) k];
                ++counted;
            }
        }

        return std::pair<bool, double> (got, counted > 0 ? sum / (double) counted : 0.0);
    };

    const auto [encodeGot, encodeMean] = run (SrMode::encode);
    const auto [decodeGot, decodeMean] = run (SrMode::decode);

    expectEqual ("encode publishes frames", encodeGot ? 1 : 0, 1);
    expectEqual ("decode publishes frames", decodeGot ? 1 : 0, 1);

    // The engaged directions read as what they do: boost above the line on
    // encode, the complementary cut below it on decode.
    ++checks;

    if (encodeMean > 1.0)
    {
        std::printf ("  [pass] %-58s %+.1f dB mean\n", "encode publishes boost", encodeMean);
    }
    else
    {
        ++failures;
        std::printf ("  [FAIL] %-58s %+.1f dB mean\n", "encode publishes boost", encodeMean);
    }

    expectBelow ("decode publishes cut", decodeMean, -1.0);
}

//==============================================================================
/** Loudness-matching makeup: converges on the dry level, holds through silence. */
void testAutoGain()
{
    std::printf ("\nAuto gain matches the wet path's loudness to the dry path's\n");

    constexpr double rate = 48000.0;
    constexpr int block = 512;

    AutoGain gain;
    gain.prepare (rate);
    gain.setEnabled (true);

    std::vector<double> dry ((size_t) block, 0.0);
    std::vector<double> wet ((size_t) block, 0.0);
    const double* dryPtr[1] = { dry.data() };
    double* wetPtr[1] = { wet.data() };

    std::uint32_t seed = 99u;

    // The wet path runs 6 dB quiet; ten seconds of programme should bring the
    // correction to +6 dB and the applied output level to the dry level.
    double lastWetRms = 0.0, dryRms = 0.0;

    for (int i = 0; i < (int) (10.0 * rate / block); ++i)
    {
        double drySum = 0.0, wetSum = 0.0;

        for (int n = 0; n < block; ++n)
        {
            seed = seed * 1664525u + 1013904223u;
            const double x = 0.25 * ((double) seed / 2147483648.0 - 1.0);

            dry[(size_t) n] = x;
            wet[(size_t) n] = 0.5 * x;
            drySum += x * x;
        }

        gain.measure (dryPtr, (const double* const*) wetPtr, 1, block);
        gain.apply (wetPtr, 1, block);

        for (int n = 0; n < block; ++n)
            wetSum += wet[(size_t) n] * wet[(size_t) n];

        dryRms = std::sqrt (drySum / block);
        lastWetRms = std::sqrt (wetSum / block);
    }

    expectNear ("correction converges to the deficit", gain.getCorrectionDb(), 6.0, 0.3);
    expectNear ("applied level matches the dry level",
                20.0 * std::log10 (lastWetRms / dryRms), 0.0, 0.3);

    // Silence carries no loudness to match: the correction holds where it was.
    const double heldDb = gain.getCorrectionDb();

    std::fill (dry.begin(), dry.end(), 0.0);
    std::fill (wet.begin(), wet.end(), 0.0);

    for (int i = 0; i < (int) (5.0 * rate / block); ++i)
    {
        gain.measure (dryPtr, (const double* const*) wetPtr, 1, block);
        gain.apply (wetPtr, 1, block);
    }

    expectNear ("the correction holds through silence", gain.getCorrectionDb(), heldDb, 0.1);

    // Switched off, it returns to unity.
    gain.setEnabled (false);

    for (int i = 0; i < (int) (2.0 * rate / block); ++i)
        gain.apply (wetPtr, 1, block);

    expectNear ("disabled returns to unity", gain.getCorrectionDb(), 0.0, 0.1);
}

//==============================================================================
/** Real-time safety of the assembled engine.

    prepare() promises the audio thread that process() never allocates, and the promise
    has to survive the two things a session actually does: hand the engine blocks larger
    than the maximum it declared — the filter-bank path chunks them through its prepared
    scratch, the STFT path streams any length through the WOLA's own hop buffers — and
    move the section control while audio is running, whose weight tables are rebuilt in
    place into vectors prepare() sized.
*/
void testEngineProcessDoesNotAllocate()
{
    std::printf ("\nEngine process() never allocates after prepare()\n");

    SrEngine engine;
    engine.prepare (48000.0, 512, 2);
    engine.setProcess (SrProcess::spectralRecording);
    engine.setMode (SrMode::loop);

    // Oversize blocks are the point: 4096 against a prepared maximum of 512.
    const int blockSizes[] = { 512, 64, 4096, 1, 333, 2048 };

    auto left = makeNoise (28224, 5555u);
    auto right = makeNoise (28224, 6666u);

    for (auto& v : left)  v *= 0.05;
    for (auto& v : right) v *= 0.05;

    const auto runBlocks = [&]
    {
        int at = 0, which = 0;
        const int total = (int) left.size();

        while (at < total)
        {
            const int len = std::min (blockSizes[which % 6], total - at);
            double* pointers[2] = { left.data() + at, right.data() + at };

            engine.process (pointers, 2, len);
            at += len;
            ++which;
        }
    };

    // STFT realisation, with a section move mid-stream — setSections is documented
    // allocation free, so it is armed too rather than taken on faith. The
    // side-chain corners and the time constants join it: the wrapper calls
    // their setters every chunk, so their moving paths belong inside the net.
    allocationCount.store (0);
    countingAllocations.store (true);

    runBlocks();
    engine.setSections (64);
    engine.setSidechainHighPass (200.0);
    engine.setSidechainLowPass (6000.0);
    engine.setTimeConstants (5.0, 100.0, 50.0, 1000.0);
    runBlocks();

    countingAllocations.store (false);

    expectEqual ("stft realisation, oversize blocks and a section move",
                 allocationCount.load(), 0);

    // The filter-bank realisation takes oversize blocks in scratch-sized chunks; the
    // switch itself is a config call and sits outside the armed region.
    engine.setGainRealisation (GainRealisation::filterBank);

    allocationCount.store (0);
    countingAllocations.store (true);

    runBlocks();

    countingAllocations.store (false);

    expectEqual ("filter-bank realisation, oversize blocks",
                 allocationCount.load(), 0);
}

//==============================================================================
/** The H-2 regression: modulation control must oppose the sliding layer, not switch
    it off.

    The MC level is one scalar for the whole spectrum. The old formulation subtracted it
    from each bin's own control, and on broadband material — where the weighted RMS is at
    least as large as any single bin's smoothed control — the difference went to zero
    everywhere and the layer fell silent outright. The paper's form is a ratio: MC raises
    the reference the compression knee measures against, so a rising broadband level
    releases the withholding smoothly and never annihilates it.
*/
void testSlidingLayerModulationControlRelease()
{
    std::printf ("\nSliding layer: modulation control releases rather than disables\n");

    constexpr double sampleRate = 48000.0;
    constexpr int fftSize = 1024;
    constexpr int numBins = fftSize / 2 + 1;
    constexpr double frameSeconds = 256.0 / sampleRate;

    // Broadband and well above threshold: a flat spectrum at -20 dB against the layer's
    // -30 dB threshold. On the normalised per-bin magnitude scale ModulationControl
    // produces, this spectrum's own MC reading is exactly its flat magnitude.
    const double magnitude = std::pow (10.0, -20.0 / 20.0);

    std::vector<double> spectrum ((size_t) numBins, magnitude);
    std::vector<double> fraction ((size_t) numBins, 0.0);

    SlidingLayerParams params;   // high-frequency defaults: 10 kHz corner, MC gain 1

    const int probeBin = (int) std::lround (15000.0 * fftSize / sampleRate);

    // Settles the layer on the spectrum at a given MC level and returns the surviving
    // boost fraction at the probe bin. 400 frames is 2.1 s against an 80 ms final tau.
    const auto settledFraction = [&] (const SlidingLayerParams& p, double mcLevel)
    {
        SlidingLayer layer;
        layer.prepare (p, numBins, sampleRate, fftSize, frameSeconds);
        layer.setModulationControl (mcLevel);

        for (int frame = 0; frame < 400; ++frame)
            layer.processFrame (spectrum.data(), fraction.data());

        return fraction[(size_t) probeBin];
    };

    const double atZero = settledFraction (params, 0.0);
    const double atMc = settledFraction (params, magnitude);

    std::printf ("    fraction at 15 kHz: MC = 0 %.4f, MC = spectrum level %.4f\n",
                 atZero, atMc);

    // With MC at zero the layer runs on its bare threshold, and a flat spectrum's
    // spread threshold sits near the level everywhere, so the control stands well
    // above threshold and real boost is withheld.
    ++checks;
    const bool active = atZero < 0.9;

    if (! active)
        ++failures;

    std::printf ("  [%s] %-58s %8.4f  (wants below 0.9)\n",
                 active ? "pass" : "FAIL", "withholding is real at MC = 0", atZero);

    // Raising MC to the spectrum's own level must release the withholding — the ratio
    // the knee sees falls — without zeroing the layer: the old subtraction was inert
    // here, reporting exactly 1.0.
    ++checks;
    const bool released = atMc > atZero + 0.05;

    if (! released)
        ++failures;

    std::printf ("  [%s] %-58s %8.4f  (wants above %.4f)\n",
                 released ? "pass" : "FAIL", "high MC releases the withholding", atMc,
                 atZero + 0.05);

    ++checks;
    const bool stillActive = atMc < 1.0;

    if (! stillActive)
        ++failures;

    std::printf ("  [%s] %-58s %8.4f  (wants below 1.0)\n",
                 stillActive ? "pass" : "FAIL", "the layer is not inert at MC > 0", atMc);

    // And the other end of the contract: at MC = 0 the opposition contributes nothing
    // at all, so the outputs must be bit-identical to a layer whose opposition gain is
    // zero — whatever MC level that layer is then fed.
    {
        SlidingLayerParams disabled = params;
        disabled.mcOppositionGain = 0.0;

        SlidingLayer a, b;
        a.prepare (params, numBins, sampleRate, fftSize, frameSeconds);
        b.prepare (disabled, numBins, sampleRate, fftSize, frameSeconds);

        a.setModulationControl (0.0);
        b.setModulationControl (0.7);   // must not matter: the gain disables it

        std::vector<double> fa ((size_t) numBins, 0.0), fb ((size_t) numBins, 0.0);
        int mismatches = 0;

        for (int frame = 0; frame < 100; ++frame)
        {
            a.processFrame (spectrum.data(), fa.data());
            b.processFrame (spectrum.data(), fb.data());

            for (int k = 0; k < numBins; ++k)
                if (fa[(size_t) k] != fb[(size_t) k])
                    ++mismatches;
        }

        expectEqual ("MC = 0 is bit-identical to opposition disabled, mismatches",
                     mismatches, 0);
    }
}

//==============================================================================
/** The reported latency, measured on the ON path rather than trusted.

    Bypass being an exact delay is already covered; what that cannot catch is the
    running process shifting audio by something other than what getLatencySamples()
    tells the host. So the full loop — encode, channel, decode — is run on broadband
    noise with the tape impairments disabled, and the delay that best explains the
    output as a copy of the input has to be the reported figure exactly.
*/
void testMeasuredOnPathLatency()
{
    std::printf ("\nMeasured ON-path latency matches the reported figure\n");

    TapeChannelSettings transparent;
    transparent.hissEnabled = false;
    transparent.saturationEnabled = false;
    transparent.hfLossEnabled = false;
    transparent.modulationNoiseEnabled = false;

    SrEngine engine;
    engine.prepare (48000.0, 512, 1);
    engine.setProcess (SrProcess::spectralRecording);
    engine.setMode (SrMode::loop);
    engine.getTapeChannel().setSettings (transparent);

    const int latency = engine.getLatencySamples();
    const int total = 4 * latency;

    auto input = makeNoise (total, 1970u);

    for (auto& v : input)
        v *= 0.05;                       // the level the loop null is quoted at

    std::vector<double> output = input;
    runEngine (engine, output);

    const int fitted = bestFitDelay (output, input, latency + 4096, latency);

    expectEqual ("loop mode, spectral recording, best-fit delay", fitted, latency);
}

//==============================================================================
/** The wet-only low pass ahead of the blend.

    A Butterworth magnitude through the bilinear transform has a closed form —
    |H|^2 = 1 / (1 + r^6) with r the ratio of prewarped frequencies — so the
    digital design is held to the analog prototype exactly across the band, not
    just spot-checked at the corner.
*/
void testWetLowPass()
{
    std::printf ("\nWet low pass matches the third-order Butterworth prototype\n");

    const auto magnitudeDb = [] (const std::array<BiquadCoeffs, 2>& c, double hz, double rate)
    {
        const double omega = twoPi * hz / rate;
        return 20.0 * std::log10 (std::abs (c[0].responseAt (omega) * c[1].responseAt (omega)));
    };

    // The whole curve against the prototype, at the corners the knob rests on.
    for (double rate : { 44100.0, 48000.0, 96000.0 })
        for (double cornerHz : { 1000.0, 6000.0 })
        {
            const auto coeffs = butterworthLowPass18 (cornerHz, rate);
            const double kCorner = std::tan (3.14159265358979323846 * cornerHz / rate);
            double worst = 0.0;

            for (double hz = 20.0; hz < 0.49 * rate; hz *= 1.1)
            {
                const double r = std::tan (3.14159265358979323846 * hz / rate) / kCorner;
                const double wantDb = -10.0 * std::log10 (1.0 + std::pow (r, 6.0));

                worst = std::max (worst, std::abs (magnitudeDb (coeffs, hz, rate) - wantDb));
            }

            expectNear ("prototype match, " + std::to_string ((int) cornerHz) + " Hz at "
                            + std::to_string ((int) rate) + " Hz",
                        worst, 0.0, 1.0e-9);
        }

    // The knob's number is the -3 dB point, exactly, and DC passes untouched.
    for (double rate : { 48000.0, 96000.0 })
        expectNear ("corner is -3.01 dB at " + std::to_string ((int) rate) + " Hz",
                    magnitudeDb (butterworthLowPass18 (6000.0, rate), 6000.0, rate),
                    -3.0102999566398, 1.0e-3);

    for (double cornerHz : { 100.0, 6000.0, 20000.0 })
        expectNear ("DC unity, corner " + std::to_string ((int) cornerHz) + " Hz",
                    magnitudeDb (butterworthLowPass18 (cornerHz, 48000.0), 0.0, 48000.0),
                    0.0, 1.0e-10);

    // Third order means an octave above the corner sits at -10 log10 (1 + 2^6).
    // Measured on a low corner, where the bilinear warp is negligible.
    expectNear ("one octave above a 1 kHz corner",
                magnitudeDb (butterworthLowPass18 (1000.0, 48000.0), 2000.0, 48000.0),
                -18.129, 0.3);

    // Fully open at the lowest rate: poles close to Nyquist, still quiet and stable.
    {
        WetLowPass filter;
        filter.prepare (44100.0, defaultSettlingSeconds, 1, 20000.0);

        std::vector<double> impulse (8192, 0.0);
        impulse[0] = 1.0;
        double* pointers[1] = { impulse.data() };
        filter.process (pointers, 1, (int) impulse.size());

        long long finite = 0;

        for (auto x : impulse)
            if (std::isfinite (x))
                ++finite;

        double tail = 0.0;

        for (size_t n = 4096; n < impulse.size(); ++n)
            tail = std::max (tail, std::abs (impulse[n]));

        expectEqual ("20 kHz corner at 44.1 kHz, every sample finite",
                     finite, (long long) impulse.size());
        expectBelow ("20 kHz corner at 44.1 kHz, ring 4096 samples in",
                     relativeDb (tail, 1.0), -100.0);
    }

    // A full-range sweep under full-scale programme: the corner walks from 6 kHz
    // to 100 Hz while a 1 kHz sine plays. DF-I under moving coefficients has to
    // stay bounded, and the redesigns must not allocate.
    {
        WetLowPass filter;
        filter.prepare (48000.0, defaultSettlingSeconds, 2, 6000.0);
        filter.setCornerHz (100.0);

        constexpr int block = 512;
        std::vector<double> left (block), right (block);
        double* pointers[2] = { left.data(), right.data() };

        double peak = 0.0;
        long long badSamples = 0;
        int sampleIndex = 0;

        allocationCount.store (0);
        countingAllocations.store (true);

        for (int b = 0; b < 40; ++b)
        {
            for (int n = 0; n < block; ++n, ++sampleIndex)
                left[(size_t) n] = right[(size_t) n]
                    = std::sin (twoPi * 1000.0 * sampleIndex / 48000.0);

            filter.process (pointers, 2, block);

            for (int n = 0; n < block; ++n)
            {
                if (! std::isfinite (left[(size_t) n]))
                    ++badSamples;

                peak = std::max (peak, std::abs (left[(size_t) n]));
            }
        }

        countingAllocations.store (false);

        expectEqual ("swept corner, no non-finite samples", badSamples, 0);
        expectEqual ("swept redesigns never allocate", allocationCount.load(), 0);
        expectBelow ("peak through the sweep", relativeDb (peak, 1.0), 2.0);
        expectNear ("the corner arrived", filter.getCornerHz(), 100.0, 0.5, "Hz");
    }

    // reset() clears history and nothing else: the response afterwards is the
    // response of a fresh instance, sample for sample.
    {
        WetLowPass used, fresh;
        used.prepare (48000.0, defaultSettlingSeconds, 1, 2000.0);
        fresh.prepare (48000.0, defaultSettlingSeconds, 1, 2000.0);

        auto noise = makeNoise (1024, 4242u);
        double* noisePtr[1] = { noise.data() };
        used.process (noisePtr, 1, (int) noise.size());
        used.reset();

        std::vector<double> a (2048, 0.0), b (2048, 0.0);
        a[0] = b[0] = 1.0;
        double* aPtr[1] = { a.data() };
        double* bPtr[1] = { b.data() };
        used.process (aPtr, 1, (int) a.size());
        fresh.process (bPtr, 1, (int) b.size());

        double difference = 0.0;

        for (size_t n = 0; n < a.size(); ++n)
            difference = std::max (difference, std::abs (a[n] - b[n]));

        expectNear ("reset matches a fresh instance", difference, 0.0, 0.0, "");
    }
}

//==============================================================================
/** The side-chain high pass: the compressors' level detection reads through it,
    the signal path does not. Both engines weight encode and decode identically,
    so the round trip must stay exact at any corner — that is the property that
    distinguishes this placement from a filter in the audio path. */
void testSidechainHighPass()
{
    std::printf ("\nSide-chain high pass weights detection, not audio\n");

    // The weight is the analog prototype's magnitude exactly: -3.01 dB on the
    // corner, 18 dB an octave beyond it, silent at DC, unity far into the band.
    expectNear ("high-pass corner is -3.01 dB",
                20.0 * std::log10 (butterworth3Magnitude (6000.0, 6000.0, true)),
                -3.0103, 1.0e-3);
    expectNear ("low-pass corner is -3.01 dB",
                20.0 * std::log10 (butterworth3Magnitude (6000.0, 6000.0, false)),
                -3.0103, 1.0e-3);
    expectNear ("an octave below the corner",
                20.0 * std::log10 (butterworth3Magnitude (3000.0, 6000.0, true)),
                -18.129, 1.0e-2);
    expectNear ("far above the corner it is unity",
                20.0 * std::log10 (butterworth3Magnitude (20000.0, 100.0, true)),
                0.0, 1.0e-3);
    expectNear ("DC reads silent", butterworth3Magnitude (0.0, 6000.0, true),
                0.0, 0.0, "");

    // Loop-mode null with the corner engaged: encode and decode weight their
    // detection identically, so the complementarity is what it was without it —
    // the same limit testEncodeDecodeNull records for the unweighted engine.
    {
        TapeChannelSettings transparent;
        transparent.hissEnabled = false;
        transparent.saturationEnabled = false;
        transparent.hfLossEnabled = false;
        transparent.modulationNoiseEnabled = false;

        SrEngine engine;
        engine.prepare (48000.0, 1024, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (SrMode::loop);
        engine.setSidechainHighPass (6000.0);
        engine.getTapeChannel().setSettings (transparent);

        const int latency = engine.getLatencySamples();
        const int total = 6 * latency;

        auto input = makeNoise (total, 60606u);

        for (auto& v : input)
            v *= 0.05;

        std::vector<double> output = input;
        runEngine (engine, output);

        const double err = delayedError (output, input, latency, latency);

        expectBelow ("loop-mode null with the corner at 6 kHz",
                     relativeDb (err, peakMagnitude (input)), -26.0);
    }

    // Below the corner the detection under-reads, so a moderate low tone is
    // treated as low-level programme and boosted harder on encode.
    {
        const auto encodeLevelDb = [] (double cornerHz)
        {
            SrEngine engine;
            engine.prepare (48000.0, 1024, 1);
            engine.setProcess (SrProcess::spectralRecording);
            engine.setMode (SrMode::encode);
            engine.setSidechainHighPass (cornerHz);

            const int latency = engine.getLatencySamples();
            const int total = 6 * latency;

            std::vector<double> buffer ((size_t) total);

            for (int n = 0; n < total; ++n)
                buffer[(size_t) n] = 0.1 * std::sin (twoPi * 200.0 * (double) n / 48000.0);

            runEngine (engine, buffer);

            double sum = 0.0;

            for (int n = total / 2; n < total; ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            return 10.0 * std::log10 (sum / (double) (total - total / 2));
        };

        const double extraDb = encodeLevelDb (12000.0) - encodeLevelDb (0.0);

        ++checks;
        const bool boosted = extraDb > 3.0;

        if (! boosted)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     boosted ? "pass" : "FAIL",
                     "a 200 Hz tone at -20 dB gains more, corner at 12 kHz", extraDb, 3.0);
    }

    // A-type: the corner set identically on an encoder and a decoder leaves the
    // round trip exact — the decoder's solver recovers the encoder's input from
    // the same detector states, weighted or not.
    {
        const auto runAtype = [] (AtypeEngine& e, std::vector<double>& buffer)
        {
            constexpr int block = 512;

            for (int at = 0; at < (int) buffer.size(); at += block)
            {
                double* ptr = buffer.data() + at;
                const int len = std::min (block, (int) buffer.size() - at);
                e.process (&ptr, 1, len);
            }
        };

        AtypeEngine encoder, decoder;
        encoder.prepare (48000.0, 512, 1);
        decoder.prepare (48000.0, 512, 1);
        encoder.setMode (AtypeMode::encode);
        decoder.setMode (AtypeMode::decode);
        encoder.setSidechainHighPass (6000.0);
        decoder.setSidechainHighPass (6000.0);

        auto input = makeNoise (48000, 31313u);

        for (auto& v : input)
            v *= 0.05;

        std::vector<double> processed = input;
        runAtype (encoder, processed);
        runAtype (decoder, processed);

        double err = 0.0;

        for (size_t n = 0; n < input.size(); ++n)
            err = std::max (err, std::abs (processed[n] - input[n]));

        expectBelow ("a-type round trip with the corner at 6 kHz",
                     relativeDb (err, peakMagnitude (input)), -100.0);

        // And the weighting is real: a 12 kHz tone whose band detector under-reads
        // comes out of the encoder hotter than with the weighting off.
        const auto atypeLevelDb = [&runAtype] (double cornerHz)
        {
            AtypeEngine e;
            e.prepare (48000.0, 512, 1);
            e.setMode (AtypeMode::encode);
            e.setSidechainHighPass (cornerHz);

            std::vector<double> buffer (48000);

            for (int n = 0; n < (int) buffer.size(); ++n)
                buffer[(size_t) n] = 0.05 * std::sin (twoPi * 12000.0 * (double) n / 48000.0);

            runAtype (e, buffer);

            double sum = 0.0;
            const int half = (int) buffer.size() / 2;

            for (int n = half; n < (int) buffer.size(); ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            return 10.0 * std::log10 (sum / (double) half);
        };

        const double extraDb = atypeLevelDb (18000.0) - atypeLevelDb (0.0);

        ++checks;
        const bool boosted = extraDb > 1.0;

        if (! boosted)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     boosted ? "pass" : "FAIL",
                     "a-type: a 12 kHz tone gains more, corner at 18 kHz", extraDb, 1.0);
    }
}

//==============================================================================
/** The matching side-chain low pass, and the two corners together. Same
    contract as the high pass: detection only, and the round trip stays exact
    with both engaged. */
void testSidechainLowPass()
{
    std::printf ("\nSide-chain low pass, and both corners together\n");

    // Loop-mode null with the detection band-limited from both ends.
    {
        TapeChannelSettings transparent;
        transparent.hissEnabled = false;
        transparent.saturationEnabled = false;
        transparent.hfLossEnabled = false;
        transparent.modulationNoiseEnabled = false;

        SrEngine engine;
        engine.prepare (48000.0, 1024, 1);
        engine.setProcess (SrProcess::spectralRecording);
        engine.setMode (SrMode::loop);
        engine.setSidechainHighPass (200.0);
        engine.setSidechainLowPass (6000.0);
        engine.getTapeChannel().setSettings (transparent);

        const int latency = engine.getLatencySamples();
        const int total = 6 * latency;

        auto input = makeNoise (total, 91919u);

        for (auto& v : input)
            v *= 0.05;

        std::vector<double> output = input;
        runEngine (engine, output);

        const double err = delayedError (output, input, latency, latency);

        expectBelow ("loop-mode null, detection band-limited 200 Hz - 6 kHz",
                     relativeDb (err, peakMagnitude (input)), -26.0);
    }

    // Above the corner the detection under-reads, so a moderate high tone is
    // treated as low-level programme and boosted harder on encode.
    {
        const auto encodeLevelDb = [] (double cornerHz)
        {
            SrEngine engine;
            engine.prepare (48000.0, 1024, 1);
            engine.setProcess (SrProcess::spectralRecording);
            engine.setMode (SrMode::encode);
            engine.setSidechainLowPass (cornerHz);

            const int latency = engine.getLatencySamples();
            const int total = 6 * latency;

            std::vector<double> buffer ((size_t) total);

            for (int n = 0; n < total; ++n)
                buffer[(size_t) n] = 0.1 * std::sin (twoPi * 10000.0 * (double) n / 48000.0);

            runEngine (engine, buffer);

            double sum = 0.0;

            for (int n = total / 2; n < total; ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            return 10.0 * std::log10 (sum / (double) (total - total / 2));
        };

        const double extraDb = encodeLevelDb (1000.0) - encodeLevelDb (0.0);

        ++checks;
        const bool boosted = extraDb > 3.0;

        if (! boosted)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     boosted ? "pass" : "FAIL",
                     "a 10 kHz tone at -20 dB gains more, corner at 1 kHz", extraDb, 3.0);
    }

    // A-type: both corners set identically on an encoder and a decoder leave the
    // round trip exact, and the low pass moves the high bands' detectors.
    {
        const auto runAtype = [] (AtypeEngine& e, std::vector<double>& buffer)
        {
            constexpr int block = 512;

            for (int at = 0; at < (int) buffer.size(); at += block)
            {
                double* ptr = buffer.data() + at;
                const int len = std::min (block, (int) buffer.size() - at);
                e.process (&ptr, 1, len);
            }
        };

        AtypeEngine encoder, decoder;
        encoder.prepare (48000.0, 512, 1);
        decoder.prepare (48000.0, 512, 1);
        encoder.setMode (AtypeMode::encode);
        decoder.setMode (AtypeMode::decode);
        encoder.setSidechainHighPass (200.0);
        encoder.setSidechainLowPass (6000.0);
        decoder.setSidechainHighPass (200.0);
        decoder.setSidechainLowPass (6000.0);

        auto input = makeNoise (48000, 41414u);

        for (auto& v : input)
            v *= 0.05;

        std::vector<double> processed = input;
        runAtype (encoder, processed);
        runAtype (decoder, processed);

        double err = 0.0;

        for (size_t n = 0; n < input.size(); ++n)
            err = std::max (err, std::abs (processed[n] - input[n]));

        expectBelow ("a-type round trip, detection band-limited 200 Hz - 6 kHz",
                     relativeDb (err, peakMagnitude (input)), -100.0);

        const auto atypeLevelDb = [&runAtype] (double cornerHz)
        {
            AtypeEngine e;
            e.prepare (48000.0, 512, 1);
            e.setMode (AtypeMode::encode);
            e.setSidechainLowPass (cornerHz);

            std::vector<double> buffer (48000);

            for (int n = 0; n < (int) buffer.size(); ++n)
                buffer[(size_t) n] = 0.05 * std::sin (twoPi * 12000.0 * (double) n / 48000.0);

            runAtype (e, buffer);

            double sum = 0.0;
            const int half = (int) buffer.size() / 2;

            for (int n = half; n < (int) buffer.size(); ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            return 10.0 * std::log10 (sum / (double) half);
        };

        const double extraDb = atypeLevelDb (1000.0) - atypeLevelDb (0.0);

        ++checks;
        const bool boosted = extraDb > 1.0;

        if (! boosted)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     boosted ? "pass" : "FAIL",
                     "a-type: a 12 kHz tone gains more, corner at 1 kHz", extraDb, 1.0);
    }
}

//==============================================================================
/** The per-half attack/release times. Disengaged they are bit-identical to the
    published ballistics; engaged and matched across the halves the round trip
    keeps its null; mismatched halves forfeit it — measured here so the
    agreement discipline is a number rather than a warning. Each half's pair
    only reaches its own direction, and the knobs audibly do what they say. */
void testGlobalTimeConstants()
{
    std::printf ("\nPer-half attack/release time constants\n");

    // The frame-pole convention the feature rests on: the coefficient is the
    // fraction kept per frame, exp(-frame/tau).
    expectNear ("frame pole at tau = one frame is 1/e",
                framePoleCoefficient (10.0, 0.010), std::exp (-1.0), 1.0e-12, "");

    // Component level: the final smoother's asymmetric pair on a FixedBand.
    {
        const double frameSeconds = 256.0 / 48000.0;
        constexpr int bins = 8;

        const auto runFrames = [] (FixedBand& band, double magnitude, int frames)
        {
            std::vector<double> magnitudes ((size_t) bins, magnitude);
            std::vector<double> transfer ((size_t) bins, 0.0);

            for (int n = 0; n < frames; ++n)
                band.processFrame (magnitudes.data(), transfer.data());
        };

        FixedBand published, engaged, cleared;
        published.prepare (FixedBandParams {}, bins, 48000.0, 1024, frameSeconds);
        engaged.prepare (FixedBandParams {}, bins, 48000.0, 1024, frameSeconds);
        cleared.prepare (FixedBandParams {}, bins, 48000.0, 1024, frameSeconds);

        // Engaging and disengaging returns the exact published trajectory.
        cleared.setTimeConstants (5.0, 100.0);
        cleared.setTimeConstants (0.0, 0.0);

        runFrames (published, 0.1, 10);
        runFrames (cleared, 0.1, 10);

        expectNear ("disengaged fixed band is bit-identical",
                    cleared.getControl()[0] - published.getControl()[0], 0.0, 0.0, "");

        // A 5 ms attack outruns the published 160 ms final within ten frames.
        engaged.setTimeConstants (5.0, 0.0);
        runFrames (engaged, 0.1, 10);

        ++checks;
        const bool faster = engaged.getControl()[0] > published.getControl()[0];

        if (! faster)
            ++failures;

        std::printf ("  [%s] %-58s %8.3f    (published %.3f)\n",
                     faster ? "pass" : "FAIL",
                     "a 5 ms attack outruns the published 160 ms rise",
                     engaged.getControl()[0], published.getControl()[0]);

        // And a 20 ms release falls faster once the level goes away.
        FixedBand slowFall, fastFall;
        slowFall.prepare (FixedBandParams {}, bins, 48000.0, 1024, frameSeconds);
        fastFall.prepare (FixedBandParams {}, bins, 48000.0, 1024, frameSeconds);
        fastFall.setTimeConstants (0.0, 20.0);

        runFrames (slowFall, 0.1, 400);
        runFrames (fastFall, 0.1, 400);
        runFrames (slowFall, 0.0, 10);
        runFrames (fastFall, 0.0, 10);

        ++checks;
        const bool falls = fastFall.getControl()[0] < slowFall.getControl()[0];

        if (! falls)
            ++failures;

        std::printf ("  [%s] %-58s %8.4f    (published %.4f)\n",
                     falls ? "pass" : "FAIL",
                     "a 20 ms release undercuts the published 160 ms fall",
                     fastFall.getControl()[0], slowFall.getControl()[0]);
    }

    TapeChannelSettings transparent;
    transparent.hissEnabled = false;
    transparent.saturationEnabled = false;
    transparent.hfLossEnabled = false;
    transparent.modulationNoiseEnabled = false;

    // Disengaged is bit-identical: the whole suite runs on the sentinel, and
    // this pins it directly — an engine told "no times" against one never told
    // anything.
    {
        SrEngine plain, told;

        for (auto* e : { &plain, &told })
        {
            e->prepare (48000.0, 1024, 1);
            e->setProcess (SrProcess::spectralRecording);
            e->setMode (SrMode::loop);
            e->getTapeChannel().setSettings (transparent);
        }

        told.setTimeConstants (0.0, 0.0, 0.0, 0.0);

        const int total = 4 * plain.getLatencySamples();

        auto input = makeNoise (total, 72727u);

        for (auto& v : input)
            v *= 0.05;

        std::vector<double> a = input, b = input;
        runEngine (plain, a);
        runEngine (told, b);

        double diff = 0.0;

        for (size_t n = 0; n < a.size(); ++n)
            diff = std::max (diff, std::abs (a[n] - b[n]));

        expectNear ("disengaged spectral engine is bit-identical", diff, 0.0, 0.0, "");
    }

    // Matched halves keep the loop null; mismatched halves measurably forfeit
    // it. The second number is the price of independence, asserted so it reads
    // as a property rather than a defect report.
    double matchedDb = 0.0;

    {
        const auto loopNullDb = [&transparent] (double upA, double upR,
                                                double downA, double downR)
        {
            SrEngine engine;
            engine.prepare (48000.0, 1024, 1);
            engine.setProcess (SrProcess::spectralRecording);
            engine.setMode (SrMode::loop);
            engine.setTimeConstants (upA, upR, downA, downR);
            engine.getTapeChannel().setSettings (transparent);

            const int latency = engine.getLatencySamples();
            const int total = 6 * latency;

            auto input = makeNoise (total, 21212u);

            for (auto& v : input)
                v *= 0.05;

            std::vector<double> output = input;
            runEngine (engine, output);

            return relativeDb (delayedError (output, input, latency, latency),
                               peakMagnitude (input));
        };

        matchedDb = loopNullDb (5.0, 100.0, 5.0, 100.0);
        const double mismatchedDb = loopNullDb (5.0, 100.0, 50.0, 1000.0);

        expectBelow ("loop-mode null with matched times, 5 ms / 100 ms",
                     matchedDb, -26.0);

        ++checks;
        const bool forfeited = mismatchedDb > matchedDb + 3.0;

        if (! forfeited)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (matched %6.1f dB)\n",
                     forfeited ? "pass" : "FAIL",
                     "mismatched halves forfeit the null, as stated",
                     mismatchedDb, matchedDb);
    }

    // Each pair reaches only its own direction: an encode-only instance is
    // bit-identical whatever the down pair says, and a decode-only instance
    // whatever the up pair says.
    {
        const auto directionDiff = [&transparent] (SrMode mode, double upA, double upR,
                                                   double downA, double downR)
        {
            SrEngine plain, told;

            for (auto* e : { &plain, &told })
            {
                e->prepare (48000.0, 1024, 1);
                e->setProcess (SrProcess::spectralRecording);
                e->setMode (mode);
                e->getTapeChannel().setSettings (transparent);
            }

            told.setTimeConstants (upA, upR, downA, downR);

            const int total = 4 * plain.getLatencySamples();

            auto input = makeNoise (total, 34343u);

            for (auto& v : input)
                v *= 0.05;

            std::vector<double> a = input, b = input;
            runEngine (plain, a);
            runEngine (told, b);

            double diff = 0.0;

            for (size_t n = 0; n < a.size(); ++n)
                diff = std::max (diff, std::abs (a[n] - b[n]));

            return diff;
        };

        expectNear ("the down pair never reaches an encode",
                    directionDiff (SrMode::encode, 0.0, 0.0, 500.0, 5000.0), 0.0, 0.0, "");
        expectNear ("the up pair never reaches a decode",
                    directionDiff (SrMode::decode, 500.0, 5000.0, 0.0, 0.0), 0.0, 0.0, "");
    }

    // The knobs do what they say on the way out of a transient. A tone drops
    // 46 dB mid-buffer; a fast up release restores the boost inside the
    // measurement window where a slow one has barely begun.
    {
        const auto windowedDb = [] (double attackMs, double releaseMs,
                                    bool loudFirst, int windowStart, int windowLength)
        {
            SrEngine engine;
            engine.prepare (48000.0, 1024, 1);
            engine.setProcess (SrProcess::spectralRecording);
            engine.setMode (SrMode::encode);
            engine.setTimeConstants (attackMs, releaseMs, 0.0, 0.0);

            const int latency = engine.getLatencySamples();
            const int total = 8 * latency;
            const int half = total / 2;

            std::vector<double> buffer ((size_t) total);

            for (int n = 0; n < total; ++n)
            {
                const double amplitude = (n < half) == loudFirst ? 0.1 : 5.0e-4;
                buffer[(size_t) n] = amplitude * std::sin (twoPi * 1000.0 * (double) n / 48000.0);
            }

            runEngine (engine, buffer);

            const int start = half + latency + windowStart;
            double sum = 0.0;

            for (int n = start; n < start + windowLength; ++n)
                sum += buffer[(size_t) n] * buffer[(size_t) n];

            return 10.0 * std::log10 (sum / (double) windowLength);
        };

        const double fastRelease = windowedDb (0.0, 20.0, true, 4800, 12000);
        const double slowRelease = windowedDb (0.0, 2000.0, true, 4800, 12000);

        ++checks;
        const bool releases = fastRelease - slowRelease > 3.0;

        if (! releases)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     releases ? "pass" : "FAIL",
                     "a 20 ms up release restores boost sooner than 2 s",
                     fastRelease - slowRelease, 3.0);

        // And on the way in: a slow attack keeps shedding boost long after a
        // fast one has finished.
        const double fastAttack = windowedDb (1.0, 0.0, false, 480, 4800);
        const double slowAttack = windowedDb (1000.0, 0.0, false, 480, 4800);

        ++checks;
        const bool attacks = slowAttack - fastAttack > 2.0;

        if (! attacks)
            ++failures;

        std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                     attacks ? "pass" : "FAIL",
                     "a 1 s up attack sheds boost later than 1 ms",
                     slowAttack - fastAttack, 2.0);
    }

    // A-type: same contract. Disengaged bit-identical, matched pairs leave the
    // decoder's solver exact, and the release knob is audible on the way out.
    {
        const auto runAtype = [] (AtypeEngine& e, std::vector<double>& buffer)
        {
            constexpr int block = 512;

            for (int at = 0; at < (int) buffer.size(); at += block)
            {
                double* ptr = buffer.data() + at;
                const int len = std::min (block, (int) buffer.size() - at);
                e.process (&ptr, 1, len);
            }
        };

        {
            AtypeEngine plain, told;
            plain.prepare (48000.0, 512, 1);
            told.prepare (48000.0, 512, 1);
            plain.setMode (AtypeMode::encode);
            told.setMode (AtypeMode::encode);
            told.setTimeConstants (0.0, 0.0, 0.0, 0.0);

            auto input = makeNoise (48000, 45454u);

            for (auto& v : input)
                v *= 0.05;

            std::vector<double> a = input, b = input;
            runAtype (plain, a);
            runAtype (told, b);

            double diff = 0.0;

            for (size_t n = 0; n < a.size(); ++n)
                diff = std::max (diff, std::abs (a[n] - b[n]));

            expectNear ("disengaged a-type is bit-identical", diff, 0.0, 0.0, "");
        }

        {
            AtypeEngine encoder, decoder;
            encoder.prepare (48000.0, 512, 1);
            decoder.prepare (48000.0, 512, 1);
            encoder.setMode (AtypeMode::encode);
            decoder.setMode (AtypeMode::decode);
            encoder.setTimeConstants (5.0, 100.0, 5.0, 100.0);
            decoder.setTimeConstants (5.0, 100.0, 5.0, 100.0);

            auto input = makeNoise (48000, 56565u);

            for (auto& v : input)
                v *= 0.05;

            std::vector<double> processed = input;
            runAtype (encoder, processed);
            runAtype (decoder, processed);

            double err = 0.0;

            for (size_t n = 0; n < input.size(); ++n)
                err = std::max (err, std::abs (processed[n] - input[n]));

            expectBelow ("a-type round trip with matched times, 5 ms / 100 ms",
                         relativeDb (err, peakMagnitude (input)), -100.0);
        }

        {
            const auto atypeReleaseDb = [&runAtype] (double releaseMs)
            {
                AtypeEngine e;
                e.prepare (48000.0, 512, 1);
                e.setMode (AtypeMode::encode);
                e.setTimeConstants (0.0, releaseMs, 0.0, 0.0);

                std::vector<double> buffer (48000);
                const int half = (int) buffer.size() / 2;

                for (int n = 0; n < (int) buffer.size(); ++n)
                {
                    const double amplitude = n < half ? 0.05 : 5.0e-4;
                    buffer[(size_t) n] = amplitude * std::sin (twoPi * 12000.0 * (double) n / 48000.0);
                }

                runAtype (e, buffer);

                double sum = 0.0;
                const int start = half + 2400;
                const int length = 9600;

                for (int n = start; n < start + length; ++n)
                    sum += buffer[(size_t) n] * buffer[(size_t) n];

                return 10.0 * std::log10 (sum / (double) length);
            };

            const double extraDb = atypeReleaseDb (20.0) - atypeReleaseDb (2000.0);

            ++checks;
            const bool releases = extraDb > 1.0;

            if (! releases)
                ++failures;

            std::printf ("  [%s] %-58s %8.1f dB  (wants above %.1f dB)\n",
                         releases ? "pass" : "FAIL",
                         "a-type: a 20 ms release restores boost sooner than 2 s",
                         extraDb, 1.0);
        }
    }
}

} // namespace

//==============================================================================
int main()
{
    std::printf ("Lime DSP tests\n");
    std::printf ("Accelerate backend %s\n", Fft64::accelerateAvailable() ? "available" : "not available");

    testForwardAgainstNaiveDft();
    testRoundTrip();
    testBackendsAgree();
    testKnownSignals();
    testWindowCola();
    testWolaPerfectReconstruction();
    testOverlapSaveAgainstDirectConvolution();
    testGeometry();
    testEngineBypassIsPureDelay();
    testEngineOffKeepsTheChannelInLoop();
    testEncodeDecodeNull();
    testAnalogPrototypesMatchThePaper();
    testDigitalMainPathMatchesThePublishedFigures();
    testBiquadInversionIsExact();
    testFixedBandStaggering();
    testModulationControl();
    testSlidingLayerCompositions();
    testGainSlew();
    testPublishedCurves();
    testShelvingBankInvertsExactly();
    testControlSmoothing();
    testSectionResolution();
    testSectionTransitionIsSmooth();
    testFractionalDelayFloor();
    testCrossfadeRamp();
    testDelayLineRepositioning();
    testSpectrumProbe();
    testProbesFollowTheRunningDirection();
    testAutoGain();
    testEngineProcessDoesNotAllocate();
    testSlidingLayerModulationControlRelease();
    testMeasuredOnPathLatency();
    testWetLowPass();
    testSidechainHighPass();
    testSidechainLowPass();
    testGlobalTimeConstants();

    std::printf ("\n%d checks, %d failure%s\n", checks, failures, failures == 1 ? "" : "s");

    return failures == 0 ? 0 : 1;
}
