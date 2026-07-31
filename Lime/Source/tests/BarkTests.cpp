/*
  ==============================================================================

    Standalone tests for Source/dsp/Bark.

    Bark depends on nothing but the standard library, so this builds on its own:

        clang++ -std=c++20 -O2 -Wall -Wextra BarkTests.cpp ../dsp/Bark.cpp \
                -o barktests && ./barktests

    The masking threshold is checked against a brute-force O(n^2) evaluation of
    maskingSpreadDb over every ordered pair, so the envelope method is held to
    numerical equality rather than to a tolerance.

  ==============================================================================
*/

#include "../dsp/Bark.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <string>
#include <vector>

//==============================================================================
// Allocation counter, so compute() can be shown to be real-time safe. Installed
// globally because a container allocating anywhere below compute() must be
// caught, not only one belonging to SpreadThreshold.
namespace
{
    std::size_t allocationCount = 0;
}

void* operator new (std::size_t size)
{
    ++allocationCount;

    if (void* p = std::malloc (size == 0 ? 1 : size))
        return p;

    throw std::bad_alloc();
}

void* operator new[] (std::size_t size)                       { return operator new (size); }
void operator delete (void* p) noexcept                       { std::free (p); }
void operator delete[] (void* p) noexcept                     { std::free (p); }
void operator delete (void* p, std::size_t) noexcept          { std::free (p); }
void operator delete[] (void* p, std::size_t) noexcept        { std::free (p); }

//==============================================================================
using namespace lime;

namespace
{

int failures = 0;
int checks = 0;

void expect (const std::string& what, bool ok, const std::string& detail = {})
{
    ++checks;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-58s %s\n", ok ? "pass" : "FAIL", what.c_str(), detail.c_str());
}

void expectBelow (const std::string& what, double measured, double limit)
{
    ++checks;
    const bool ok = measured <= limit;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-58s %12.4g  (limit %8.3g)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measured, limit);
}

void expectNear (const std::string& what, double measured, double target, double tolerance)
{
    ++checks;
    const bool ok = std::abs (measured - target) <= tolerance;

    if (! ok)
        ++failures;

    std::printf ("  [%s] %-58s %12.6g  (target %8.4g +/- %.3g)\n",
                 ok ? "pass" : "FAIL", what.c_str(), measured, target, tolerance);
}

//==============================================================================
void testBarkMapping()
{
    std::printf ("\nBark mapping\n");

    // The stated landmarks. Traunmuller's uncorrected form gives 0.7715 at
    // 100 Hz, 13.010 at 2 kHz and 19.679 at 6 kHz.
    expectNear ("bark at 100 Hz", barkFromHz (100.0), 0.7715, 0.001);
    expectNear ("bark at 2 kHz", barkFromHz (2000.0), 13.0104, 0.001);
    expectNear ("bark at 6 kHz", barkFromHz (6000.0), 19.6793, 0.001);

    // The looser statements the module is meant to satisfy.
    expectNear ("bark at 100 Hz is roughly 1", barkFromHz (100.0), 1.0, 0.3);
    expectNear ("bark at 2 kHz is roughly 13", barkFromHz (2000.0), 13.0, 0.5);
    expectNear ("bark at 6 kHz is roughly 20", barkFromHz (6000.0), 20.0, 0.5);

    expectNear ("bark at 0 Hz", barkFromHz (0.0), -0.53, 1.0e-15);

    // Round trip both ways over the audio band and beyond it.
    double worstHz = 0.0, worstBark = 0.0;

    for (int i = 0; i <= 4000; ++i)
    {
        const double hz = 10.0 * std::pow (10.0, 3.4 * (double) i / 4000.0); // 10 Hz .. 25 kHz
        const double back = hzFromBark (barkFromHz (hz));

        worstHz = std::max (worstHz, std::abs (back - hz) / hz);
    }

    for (int i = 0; i <= 4000; ++i)
    {
        const double bark = -0.5 + 26.0 * (double) i / 4000.0;
        const double back = barkFromHz (hzFromBark (bark));

        worstBark = std::max (worstBark, std::abs (back - bark));
    }

    expectBelow ("hz -> bark -> hz relative error", worstHz, 1.0e-13);
    expectBelow ("bark -> hz -> bark absolute error", worstBark, 1.0e-12);

    // Monotone and finite everywhere, including at and past the asymptote.
    bool monotone = true;
    double previous = -1.0e300;

    for (int i = 0; i <= 20000; ++i)
    {
        const double hz = 24000.0 * (double) i / 20000.0;
        const double z = barkFromHz (hz);

        monotone = monotone && (z > previous);
        previous = z;
    }

    expect ("barkFromHz strictly increasing, 0 .. 24 kHz", monotone);

    const double atAsymptote = hzFromBark (barkAsymptote);
    const double pastAsymptote = hzFromBark (barkAsymptote + 5.0);

    expect ("hzFromBark finite at the asymptote",
            std::isfinite (atAsymptote) && atAsymptote > 0.0,
            "= " + std::to_string (atAsymptote));
    expect ("hzFromBark finite past the asymptote",
            std::isfinite (pastAsymptote) && pastAsymptote > 0.0);
}

//==============================================================================
void testSpreadFunction()
{
    std::printf ("\nSpread function\n");

    // Zero at zero distance, exactly, at every masker position and level.
    double worstAtZero = 0.0;

    for (double z = 0.0; z <= 24.0; z += 0.25)
        for (double level = -80.0; level <= 10.0; level += 5.0)
            worstAtZero = std::max (worstAtZero, std::abs (maskingSpreadDb (z, z, level)));

    expectBelow ("spread is exactly 0 dB at zero distance", worstAtZero, 0.0);

    // Never positive, and monotone away from the masker in both directions.
    bool nonPositive = true, monotoneDown = true, monotoneUp = true;

    for (double z = 1.0; z <= 22.0; z += 0.5)
    {
        for (double level = -80.0; level <= 10.0; level += 10.0)
        {
            double previousDown = 0.0, previousUp = 0.0;

            for (int i = 1; i <= 400; ++i)
            {
                const double distance = 0.02 * (double) i;

                const double down = maskingSpreadDb (z, z - distance, level);
                const double up = maskingSpreadDb (z, z + distance, level);

                nonPositive = nonPositive && down <= 0.0 && up <= 0.0;
                monotoneDown = monotoneDown && down < previousDown;
                monotoneUp = monotoneUp && up < previousUp;

                previousDown = down;
                previousUp = up;
            }
        }
    }

    expect ("spread is never positive", nonPositive);
    expect ("spread decreases monotonically below the masker", monotoneDown);
    expect ("spread decreases monotonically above the masker", monotoneUp);

    // Downward slope is the published constant, exactly.
    double worstDownSlope = 0.0;

    for (double z = 1.0; z <= 22.0; z += 0.5)
        for (double distance = 0.1; distance <= 8.0; distance += 0.1)
            worstDownSlope = std::max (worstDownSlope,
                                       std::abs (maskingSpreadDb (z, z - distance, -20.0)
                                                 / distance + 27.0));

    expectBelow ("downward slope is -27 dB/Bark", worstDownSlope, 1.0e-12);

    // Upward slope: inside the bounds, and strictly shallower as level rises.
    bool inBounds = true, shallowerWithLevel = true;
    double widest = 0.0, narrowest = 100.0;

    for (double z = 0.0; z <= 24.0; z += 0.5)
    {
        double previous = 1.0e300;

        for (double level = -120.0; level <= 40.0; level += 0.5)
        {
            const double slope = upwardSpreadSlopeDbPerBark (z, level);

            inBounds = inBounds
                    && slope > barkSpread::upwardSlopeShallowestDbPerBark
                    && slope < barkSpread::upwardSlopeSteepestDbPerBark;
            shallowerWithLevel = shallowerWithLevel && slope < previous;
            previous = slope;

            widest = std::max (widest, slope);
            narrowest = std::min (narrowest, slope);
        }
    }

    expect ("upward slope magnitude lies strictly within 6 .. 15 dB/Bark", inBounds,
            "range " + std::to_string (narrowest) + " .. " + std::to_string (widest));
    expect ("upward slope gets shallower as masker level rises", shallowerWithLevel);

    // The same statement expressed through maskingSpreadDb itself: a louder
    // masker attenuates less at a fixed distance above it.
    bool louderSpreadsFurther = true;
    double previousSpread = -1.0e300;

    for (double level = -100.0; level <= 20.0; level += 0.5)
    {
        const double spread = maskingSpreadDb (12.0, 16.0, level);

        louderSpreadsFurther = louderSpreadsFurther && spread > previousSpread;
        previousSpread = spread;
    }

    expect ("louder maskers spread further upward", louderSpreadsFurther);

    // Slope over the range SR itself works in.
    std::printf ("       upward slope at 1 kHz: %5.2f dB/Bark at -62 dB, "
                 "%5.2f at -30 dB, %5.2f at 0 dB\n",
                 upwardSpreadSlopeDbPerBark (barkFromHz (1000.0), -62.0),
                 upwardSpreadSlopeDbPerBark (barkFromHz (1000.0), -30.0),
                 upwardSpreadSlopeDbPerBark (barkFromHz (1000.0), 0.0));

    // Continuity in all three arguments: the largest step over a fine grid must
    // scale with the grid, so compare against a generous multiple of it.
    constexpr double step = 1.0e-6;
    double worstProbeStep = 0.0, worstMaskerStep = 0.0, worstLevelStep = 0.0;

    for (double d = -6.0; d <= 6.0; d += 0.001)
    {
        const double a = maskingSpreadDb (10.0, 10.0 + d, -30.0);
        const double b = maskingSpreadDb (10.0, 10.0 + d + step, -30.0);
        worstProbeStep = std::max (worstProbeStep, std::abs (b - a));

        const double c = maskingSpreadDb (10.0 + d, 14.0, -30.0);
        const double e = maskingSpreadDb (10.0 + d + step, 14.0, -30.0);
        worstMaskerStep = std::max (worstMaskerStep, std::abs (e - c));
    }

    for (double level = -200.0; level <= 100.0; level += 0.001)
    {
        const double a = maskingSpreadDb (10.0, 15.0, level);
        const double b = maskingSpreadDb (10.0, 15.0, level + step);
        worstLevelStep = std::max (worstLevelStep, std::abs (b - a));
    }

    expectBelow ("continuous in probe position", worstProbeStep, 40.0 * step);
    expectBelow ("continuous in masker position", worstMaskerStep, 40.0 * step);
    expectBelow ("continuous in masker level", worstLevelStep, 40.0 * step);
}

//==============================================================================
/** Brute-force reference: for each bin, the largest masking level from any other
    bin, computed straight from maskingSpreadDb. */
std::vector<double> bruteForceThreshold (const std::vector<double>& levelsDb,
                                         const std::vector<double>& barks)
{
    const int n = (int) levelsDb.size();
    std::vector<double> out ((size_t) n, barkSpread::minimumLevelDb);

    for (int j = 0; j < n; ++j)
    {
        double best = barkSpread::minimumLevelDb;

        for (int i = 0; i < n; ++i)
        {
            if (i == j)
                continue;

            const double level = std::min (std::max (levelsDb[(size_t) i],
                                                     barkSpread::minimumLevelDb),
                                           barkSpread::maximumLevelDb);

            best = std::max (best, level + maskingSpreadDb (barks[(size_t) i],
                                                            barks[(size_t) j], level));
        }

        out[(size_t) j] = best;
    }

    return out;
}

double worstDifference (const std::vector<double>& a, const std::vector<double>& b)
{
    double worst = 0.0;

    for (size_t i = 0; i < a.size(); ++i)
        worst = std::max (worst, std::abs (a[i] - b[i]));

    return worst;
}

void testSpreadThresholdAgainstBruteForce()
{
    std::printf ("\nSpreadThreshold versus brute force\n");

    struct Case { int fftSize; double sampleRate; const char* name; };

    const Case cases[] = { {  512, 48000.0, "512 @ 48k"  },
                           { 1024, 44100.0, "1024 @ 44.1k" },
                           { 2048, 96000.0, "2048 @ 96k" },
                           { 8192, 48000.0, "8192 @ 48k" } };   // the production LF size

    std::mt19937_64 rng (12345u);
    std::uniform_real_distribution<double> dist (-120.0, 0.0);

    for (const auto& c : cases)
    {
        const int numBins = c.fftSize / 2 + 1;

        SpreadThreshold spread;
        spread.prepare (numBins, c.sampleRate, c.fftSize);

        std::vector<double> barks (spread.getBinBarks(), spread.getBinBarks() + numBins);
        std::vector<double> levels ((size_t) numBins);
        std::vector<double> got ((size_t) numBins);

        double worst = 0.0;

        // Random spectra.
        for (int trial = 0; trial < 4; ++trial)
        {
            for (auto& v : levels)
                v = dist (rng);

            spread.compute (levels.data(), got.data());
            worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));
        }

        // Single tone: the case where the widest Bark distances actually decide
        // the result, so any windowing or slope quantisation would show up here.
        for (int tone : { 1, numBins / 8, numBins / 2, numBins - 2 })
        {
            std::fill (levels.begin(), levels.end(), -140.0);
            levels[(size_t) tone] = 0.0;

            spread.compute (levels.data(), got.data());
            worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));
        }

        // Two loud tones of different level, which is what makes the upper
        // envelope non-trivial: the quieter one has the steeper upward slope.
        std::fill (levels.begin(), levels.end(), -160.0);
        levels[(size_t) (numBins / 10)] = -60.0;
        levels[(size_t) (numBins / 4)] = -6.0;

        spread.compute (levels.data(), got.data());
        worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));

        // Monotone ramps in both directions. A rising ramp means every masker is
        // louder, and so shallower, than the one before it, which is the case the
        // dominated-line shortcut must never take.
        for (int direction : { 1, -1 })
        {
            for (int k = 0; k < numBins; ++k)
            {
                const double t = (double) k / (double) (numBins - 1);
                levels[(size_t) k] = direction > 0 ? -120.0 + 120.0 * t : -120.0 * t;
            }

            spread.compute (levels.data(), got.data());
            worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));
        }

        // A pink-ish floor with a harmonic series on top: the shape the sliding
        // layer will actually be handed.
        for (int k = 0; k < numBins; ++k)
            levels[(size_t) k] = -60.0 - 10.0 * std::log10 (std::max (k, 1));

        for (int harmonic = 1; harmonic <= 12; ++harmonic)
        {
            const int bin = harmonic * numBins / 64;

            if (bin < numBins)
                levels[(size_t) bin] += 45.0;
        }

        spread.compute (levels.data(), got.data());
        worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));

        // Digital silence expressed as -inf, which must not produce a NaN.
        std::fill (levels.begin(), levels.end(), -std::numeric_limits<double>::infinity());
        levels[(size_t) (numBins / 3)] = -12.0;

        spread.compute (levels.data(), got.data());

        bool finite = true;

        for (auto v : got)
            finite = finite && std::isfinite (v);

        expect (std::string ("no NaN or infinity from -inf input levels, ") + c.name, finite);

        worst = std::max (worst, worstDifference (got, bruteForceThreshold (levels, barks)));

        expectBelow (std::string ("exact all-pairs maximum, ") + c.name, worst, 1.0e-9);
    }

    // Degenerate sizes must not crash or read out of range.
    for (int n : { 0, 1, 2, 3 })
    {
        SpreadThreshold spread;
        spread.prepare (n, 48000.0, 8);

        std::vector<double> levels ((size_t) std::max (n, 1), -20.0);
        std::vector<double> got ((size_t) std::max (n, 1), 0.0);

        spread.compute (levels.data(), got.data());

        expect ("numBins = " + std::to_string (n) + " is handled",
                spread.getNumBins() == n);
    }
}

//==============================================================================
void testSingleToneSkirt()
{
    std::printf ("\nSingle-tone skirt\n");

    constexpr int fftSize = 8192;
    constexpr double sampleRate = 48000.0;
    constexpr int numBins = fftSize / 2 + 1;
    constexpr double toneHz = 1000.0;
    constexpr double toneDb = -20.0;

    SpreadThreshold spread;
    spread.prepare (numBins, sampleRate, fftSize);

    const double* barks = spread.getBinBarks();

    // A plausible noise floor rather than silence, so the skirt has to stand out
    // of something.
    std::vector<double> levels ((size_t) numBins, -140.0);
    std::vector<double> got ((size_t) numBins);

    const int toneBin = (int) std::lround (toneHz * (double) fftSize / sampleRate);
    levels[(size_t) toneBin] = toneDb;

    spread.compute (levels.data(), got.data());

    const double toneBark = barks[(size_t) toneBin];

    // Locate the bins nearest to a given Bark offset from the tone.
    const auto binAtOffset = [&] (double offset)
    {
        const double target = toneBark + offset;
        int best = 0;

        for (int k = 1; k < numBins; ++k)
            if (std::abs (barks[(size_t) k] - target) < std::abs (barks[(size_t) best] - target))
                best = k;

        return best;
    };

    // The expected slopes: the tone's own level fixes the upward one.
    const double upSlope = upwardSpreadSlopeDbPerBark (toneBark, toneDb);

    std::printf ("       tone at bin %d (%.1f Hz, %.3f Bark), upward slope %.2f dB/Bark\n",
                 toneBin, (double) toneBin * sampleRate / (double) fftSize, toneBark, upSlope);

    for (double offset : { 0.5, 1.0, 2.0, 4.0 })
    {
        const int below = binAtOffset (-offset);
        const int above = binAtOffset (offset);

        const double barkBelow = toneBark - barks[(size_t) below];
        const double barkAbove = barks[(size_t) above] - toneBark;

        const double expectedBelow = toneDb - barkSpread::downwardSlopeDbPerBark * barkBelow;
        const double expectedAbove = toneDb - upSlope * barkAbove;

        expectNear ("skirt " + std::to_string (offset) + " Bark below the tone",
                    got[(size_t) below], expectedBelow, 1.0e-9);
        expectNear ("skirt " + std::to_string (offset) + " Bark above the tone",
                    got[(size_t) above], expectedAbove, 1.0e-9);

        // The asymmetry is the difference of the two slopes times the distance,
        // allowing for the bins not landing exactly on the requested offsets.
        const double expectedGap = (barkSpread::downwardSlopeDbPerBark - upSlope) * offset;

        expect ("skirt is higher above than below at " + std::to_string (offset) + " Bark",
                got[(size_t) above] - got[(size_t) below] > 0.9 * expectedGap,
                std::to_string (got[(size_t) above] - got[(size_t) below]) + " dB gap, expected "
                    + std::to_string (expectedGap));
    }

    // The tone's own bin is masked only by its neighbours, so it must sit well
    // below the tone: this is what makes it a detectable dominant.
    expect ("the tone's own bin is far below its own level",
            got[(size_t) toneBin] < toneDb - 100.0,
            "threshold " + std::to_string (got[(size_t) toneBin]) + " dB");

    // Monotone away from the tone, over the region where the tone rather than the
    // noise floor sets the threshold. Beyond that the floor takes over and the
    // threshold creeps back up, because the Bark spacing between adjacent bins
    // narrows towards Nyquist - which is right, not a defect.
    constexpr double toneDominates = -130.0;

    bool fallsBelow = true, fallsAbove = true;

    for (int k = 1; k < toneBin; ++k)
        if (got[(size_t) k] > toneDominates)
            fallsBelow = fallsBelow && got[(size_t) k] >= got[(size_t) k - 1];

    for (int k = toneBin + 2; k < numBins && got[(size_t) k] > toneDominates; ++k)
        fallsAbove = fallsAbove && got[(size_t) k] <= got[(size_t) k - 1];

    expect ("skirt rises monotonically towards the tone from below", fallsBelow);
    expect ("skirt falls monotonically away from the tone above", fallsAbove);

    // A bin one Bark above the tone is masked; a bin one Bark below is not, at
    // least not by anything like as much. Print the shape for the record.
    std::printf ("       skirt, dB relative to the tone:");

    for (double offset : { -4.0, -2.0, -1.0, 1.0, 2.0, 4.0 })
        std::printf ("  %+.0f Bk %6.1f", offset, got[(size_t) binAtOffset (offset)] - toneDb);

    std::printf ("\n");
}

//==============================================================================
void testRealTimeSafety()
{
    std::printf ("\nReal-time safety\n");

    constexpr int fftSize = 4096;
    constexpr int numBins = fftSize / 2 + 1;

    SpreadThreshold spread;
    spread.prepare (numBins, 48000.0, fftSize);

    std::mt19937_64 rng (777);
    std::uniform_real_distribution<double> dist (-120.0, 0.0);

    std::vector<double> levels ((size_t) numBins);
    std::vector<double> out ((size_t) numBins);

    for (auto& v : levels)
        v = dist (rng);

    // Warm up once outside the counted region, then count.
    spread.compute (levels.data(), out.data());

    const std::size_t before = allocationCount;

    for (int i = 0; i < 200; ++i)
        spread.compute (levels.data(), out.data());

    const std::size_t during = allocationCount - before;

    expect ("compute() performs no allocation", during == 0,
            std::to_string (during) + " allocations over 200 calls");

    // Free functions too, since the sliding layer calls them directly.
    const std::size_t beforeFunctions = allocationCount;
    double accumulator = 0.0;

    for (int i = 0; i < 10000; ++i)
    {
        const double z = 0.0026 * (double) i;
        accumulator += barkFromHz (z * 1000.0) + hzFromBark (z)
                     + maskingSpreadDb (z, z + 1.0, -30.0)
                     + upwardSpreadSlopeDbPerBark (z, -30.0);
    }

    // Read the counter into a local before calling expect(): argument evaluation
    // order is unspecified, and building the message string allocates.
    const std::size_t duringFunctions = allocationCount - beforeFunctions;

    expect ("free functions perform no allocation",
            duringFunctions == 0 && std::isfinite (accumulator),
            std::to_string (duringFunctions) + " allocations");
}

} // namespace

//==============================================================================
int main()
{
    std::printf ("Bark tests\n");

    testBarkMapping();
    testSpreadFunction();
    testSpreadThresholdAgainstBruteForce();
    testSingleToneSkirt();
    testRealTimeSafety();

    std::printf ("\n%d checks, %d failures\n", checks, failures);

    return failures == 0 ? 0 : 1;
}
