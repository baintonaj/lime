/*
  ==============================================================================

    Bark — critical-band mapping, two-slope spread function, and the exact
    all-pairs masking threshold.

  ==============================================================================
*/

#include "Bark.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lime
{

namespace
{
    constexpr double barkNumerator = 26.81;
    constexpr double barkOffset = 0.53;
    constexpr double barkCornerHz = 1960.0;

    constexpr double negativeInfinity = -std::numeric_limits<double>::infinity();

    double clampLevelDb (double db) noexcept
    {
        // NaN falls through unchanged, which is intended: a NaN level is a caller
        // bug and hiding it behind a clamp would make it harder to find.
        return std::min (std::max (db, barkSpread::minimumLevelDb),
                         barkSpread::maximumLevelDb);
    }

    /** Terhardt's raw upper-slope magnitude, before saturation, for a masker
        frequency in Hz and a level on this module's own dB scale. Split into a
        level-independent base so SpreadThreshold can precompute it per bin. */
    double rawSlopeBase (double maskerHz) noexcept
    {
        const double hz = std::max (maskerHz, barkSpread::minimumSlopeFrequencyHz);

        return barkSpread::terhardtConstantDb
             + barkSpread::terhardtFrequencyTermHzDb / hz
             - barkSpread::terhardtLevelCoefficient * barkSpread::referenceLevelSpl;
    }

    double rawSlope (double base, double clampedLevelDb) noexcept
    {
        return base - barkSpread::terhardtLevelCoefficient * clampedLevelDb;
    }

    /** Smooth, strictly increasing map of the raw slope magnitude onto the open
        interval between the two published bounds.

        x / sqrt(1 + x^2) is chosen over a logistic for two reasons: it costs one
        square root rather than an exponential, which matters because this runs
        once per bin per hop, and it is odd, so scaling its argument by the
        half-range gives unit gain at the interval midpoint with no further
        constants. Unit gain there is the point of the whole construction - it
        means the saturation is invisible wherever Terhardt's line already lies
        inside the bounds. */
    double saturateSlope (double raw) noexcept
    {
        constexpr double centre = 0.5 * (barkSpread::upwardSlopeSteepestDbPerBark
                                       + barkSpread::upwardSlopeShallowestDbPerBark);
        constexpr double halfRange = 0.5 * (barkSpread::upwardSlopeSteepestDbPerBark
                                          - barkSpread::upwardSlopeShallowestDbPerBark);

        const double x = (raw - centre) / halfRange;

        return centre + halfRange * x / std::sqrt (1.0 + x * x);
    }
}

//==============================================================================
double barkFromHz (double hz) noexcept
{
    return barkNumerator * hz / (barkCornerHz + hz) - barkOffset;
}

double hzFromBark (double bark) noexcept
{
    // z = 26.81*f/(1960+f) - 0.53  =>  f = 1960*(z + 0.53) / (26.81 - z - 0.53).
    // The denominator vanishes at the asymptote; back off by an amount that is
    // large enough to stay well clear of it in double precision and small enough
    // that the mapping is still monotone up to it. 1.0e-9 sits about six orders
    // of magnitude above the rounding granularity of the subtraction (ULP near
    // 26.28 is ~4e-15) yet corresponds to roughly 5e13 Hz, so any value in that
    // wide window is equivalent: every affected input already maps so far above
    // the audio band that callers treat the result identically.
    constexpr double smallestDenominator = 1.0e-9;

    const double numerator = barkCornerHz * (bark + barkOffset);
    const double denominator = std::max (barkAsymptote - bark, smallestDenominator);

    return numerator / denominator;
}

//==============================================================================
double upwardSpreadSlopeDbPerBark (double maskerBark, double maskerLevelDb) noexcept
{
    const double base = rawSlopeBase (hzFromBark (maskerBark));

    return saturateSlope (rawSlope (base, clampLevelDb (maskerLevelDb)));
}

double maskingSpreadDb (double maskerBark, double probeBark, double maskerLevelDb) noexcept
{
    const double distance = probeBark - maskerBark;

    // Both branches are zero at zero distance, so the pair is continuous there
    // whichever way the comparison falls.
    if (distance <= 0.0)
        return barkSpread::downwardSlopeDbPerBark * distance;

    return -upwardSpreadSlopeDbPerBark (maskerBark, maskerLevelDb) * distance;
}

//==============================================================================
void SpreadThreshold::prepare (int numBinsIn, double sampleRate, int fftSize)
{
    numBins = std::max (0, numBinsIn);

    binBark.assign ((size_t) numBins, 0.0);
    slopeBase.assign ((size_t) numBins, 0.0);

    // 4n nodes is the standard bound for a segment tree over [0, n-1] addressed
    // from the root as 1.
    envelope.assign ((size_t) std::max (4 * numBins, 1), Line {});

    // Every node now stamped 0, and the first compute() bumps to 1, so the tree starts
    // empty without having to say so twice.
    generation = 0;

    if (numBins == 0 || sampleRate <= 0.0 || fftSize <= 0)
        return;

    const double binWidthHz = sampleRate / (double) fftSize;

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * binWidthHz;

        binBark[(size_t) k] = barkFromHz (hz);
        slopeBase[(size_t) k] = rawSlopeBase (hz);
    }
}

SpreadThreshold::Line SpreadThreshold::lineAt (int node) const noexcept
{
    const Line& stored = envelope[(size_t) node];

    if (stored.generation == generation)
        return stored;

    return Line { negativeInfinity, 0.0, generation };
}

void SpreadThreshold::insertLine (Line line) const noexcept
{
    line.generation = generation;

    int node = 1, low = 0, high = numBins - 1;

    for (;;)
    {
        const int mid = (low + high) / 2;
        Line& resident = envelope[(size_t) node];

        // Left over from an earlier frame, so it holds nothing. Normalising it here is
        // what makes bumping the counter equivalent to clearing the tree, and it only
        // ever costs a store on the O(log n) nodes an insert would visit anyway.
        if (resident.generation != generation)
            resident = Line { negativeInfinity, 0.0, generation };

        const bool betterAtLow = line.at (binBark[(size_t) low]) > resident.at (binBark[(size_t) low]);
        const bool betterAtMid = line.at (binBark[(size_t) mid]) > resident.at (binBark[(size_t) mid]);

        // The line that wins the midpoint stays; the other can only win on the
        // half where the two disagree, so at most one child is ever visited.
        if (betterAtMid)
            std::swap (line, resident);

        if (low == high)
            return;

        if (betterAtLow != betterAtMid)
        {
            node = 2 * node;
            high = mid;
        }
        else
        {
            node = 2 * node + 1;
            low = mid + 1;
        }
    }
}

SpreadThreshold::Line SpreadThreshold::queryEnvelope (int bin) const noexcept
{
    const double z = binBark[(size_t) bin];

    int node = 1, low = 0, high = numBins - 1;
    Line best = lineAt (1);
    double bestValue = best.at (z);

    while (low != high)
    {
        const int mid = (low + high) / 2;

        if (bin <= mid)
        {
            node = 2 * node;
            high = mid;
        }
        else
        {
            node = 2 * node + 1;
            low = mid + 1;
        }

        const Line candidate = lineAt (node);
        const double value = candidate.at (z);

        if (value > bestValue)
        {
            best = candidate;
            bestValue = value;
        }
    }

    return best;
}

void SpreadThreshold::compute (const double* binLevelsDb, double* thresholdDbOut) const noexcept
{
    if (numBins <= 0)
        return;

    // Descending pass. Contribution of masker i to probe j < i is
    // (level[i] - slope*bark[i]) + slope*bark[j], so one running maximum of the
    // first bracket, taken downwards, covers every masker above every probe.
    {
        const double slope = barkSpread::downwardSlopeDbPerBark;
        double running = negativeInfinity;

        for (int j = numBins - 1; j >= 0; --j)
        {
            thresholdDbOut[j] = running + slope * binBark[(size_t) j];

            running = std::max (running,
                                clampLevelDb (binLevelsDb[j]) - slope * binBark[(size_t) j]);
        }
    }

    // Ascending pass. Each masker is one straight line in Bark and the threshold
    // is the upper envelope of the lines belonging to lower bins, so bin j is
    // queried before its own line joins the tree.
    // Empties the tree without touching it: see the note on `generation`.
    ++generation;

    for (int j = 0; j < numBins; ++j)
    {
        const double bark = binBark[(size_t) j];
        const Line incumbent = queryEnvelope (j);

        thresholdDbOut[j] = std::max (std::max (thresholdDbOut[j], incumbent.at (bark)),
                                      barkSpread::minimumLevelDb);

        const double level = clampLevelDb (binLevelsDb[j]);
        const double slope = saturateSlope (rawSlope (slopeBase[(size_t) j], level));
        const Line line { level + slope * bark, -slope };

        // A line that is already below the incumbent at its own Bark position and
        // decays at least as fast can never overtake it, so it need not be
        // inserted at all. Exact, not a heuristic, and it fires often: quiet bins
        // have both the lower level and the steeper slope.
        if (line.at (bark) <= incumbent.at (bark) && line.slope <= incumbent.slope)
            continue;

        insertLine (line);
    }
}

} // namespace lime
