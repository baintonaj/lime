/*
  ==============================================================================

    Bark — critical-band scale and the psychoacoustic masking spread function.

    Lime's sliding layer derives its control from the masking spread threshold
    evaluated at every bin — the skirt everything else in the spectrum casts
    there (plan decision 2). Local-maximum detection enters only as a metering
    count of dominants, not as a selection step. So the layer needs two things:
    a warping of the frequency axis onto critical bands, and a rule for how much
    a bin at one critical-band position masks a bin at another. Both live here.

    Nothing in this file depends on JUCE, on the transforms, or on any other Lime
    header, so it compiles and is tested standalone. Double throughout.

  ==============================================================================
*/

#pragma once

#include <cstdint>
#include <vector>

namespace lime
{

//==============================================================================
/** Value the Bark mapping approaches as frequency goes to infinity. */
inline constexpr double barkAsymptote = 26.81 - 0.53;

/** Frequency in Hz to critical-band rate in Bark, Traunmuller (1990):

        z = 26.81 * f / (1960 + f) - 0.53

    Traunmuller's low- and high-end corrections are deliberately not applied. The
    plain form already reproduces the textbook band numbers closely enough for
    peak selection (0.77 Bark at 100 Hz, 13.01 at 2 kHz, 19.68 at 6 kHz), and
    being a Mobius function of f it is strictly monotone and exactly invertible
    over the whole positive frequency axis, which the piecewise-corrected form is
    not. Exact invertibility is what lets the spread function be specified in
    Bark and evaluated in Hz without accumulating a mapping error.

    Expects hz >= 0; f = 0 maps to -0.53 rather than to zero, which is the
    formula's own value and is harmless because only differences in Bark are ever
    used.
*/
double barkFromHz (double hz) noexcept;

/** Exact inverse of barkFromHz for bark < barkAsymptote. There is no positive
    frequency at or above the asymptote, so such arguments are clamped to just
    below it: the result stays finite and monotone instead of dividing by zero.
*/
double hzFromBark (double bark) noexcept;

//==============================================================================
/** Constants of the two-slope spread function. Named rather than buried so the
    sliding layer and the tests can refer to the same numbers.
*/
namespace barkSpread
{
    /** Slope on the low-frequency side of a masker. Level independent, as
        published. */
    inline constexpr double downwardSlopeDbPerBark = 27.0;

    /** Bounds on the magnitude of the high-frequency-side slope. The published
        range is "about 6 to 15 dB per Bark" depending on masker level. */
    inline constexpr double upwardSlopeShallowestDbPerBark = 6.0;
    inline constexpr double upwardSlopeSteepestDbPerBark = 15.0;

    /** Terhardt's upper-slope law, S2 = -(24 + 230/f - 0.2*L) dB/Bark, with f in
        Hz and L in dB SPL. */
    inline constexpr double terhardtConstantDb = 24.0;
    inline constexpr double terhardtFrequencyTermHzDb = 230.0;
    inline constexpr double terhardtLevelCoefficient = 0.2;

    /** Sound pressure level in dB assumed for a level of 0 dB at this module's
        interface. Terhardt's law is written in dB SPL; Lime's bin levels are dB
        relative to reference level, so the two need tying together somewhere and
        this is the single place it happens.

        90 dB SPL is chosen because it puts SR's own thresholds in sensible
        places: the -30 dB first threshold lands at 60 dB SPL, and the -62 dB
        third threshold lands at 28 dB SPL, close to the audibility floor of a
        quiet room, which is precisely the region SR was designed to protect.
        Over that -62..0 dB span the upper slope moves from about 14.4 down to
        about 7.4 dB/Bark, so the whole published range is exercised by the
        signal levels SR actually deals with.
    */
    inline constexpr double referenceLevelSpl = 90.0;

    /** Lowest frequency used in the 230/f term. Bins below this exist at every
        transform size and the term would otherwise run away towards DC. */
    inline constexpr double minimumSlopeFrequencyHz = 20.0;

    /** Levels are clamped into this range before use. The floor doubles as the
        value SpreadThreshold reports for a bin that nothing masks, so it is far
        below anything a real spectrum reaches; the ceiling exists only so that
        an infinite or absurd input cannot turn into a NaN slope. */
    inline constexpr double minimumLevelDb = -1000.0;
    inline constexpr double maximumLevelDb = 300.0;
}

//==============================================================================
/** Magnitude in dB per Bark of the spread towards higher frequencies, for a
    masker at `maskerBark` whose own level is `maskerLevelDb`.

    Terhardt (1979), as reprinted in Painter & Spanias, "Perceptual Coding of
    Digital Audio", Proc. IEEE 88(4) 2000, gives

        S2 = -24 - 230/f + 0.2 * L    dB/Bark

    which is the standard published form of the result that loud maskers spread
    further upward. Taken literally it leaves the stated 6..15 dB/Bark range at
    both ends (about 24 dB/Bark for an inaudible masker, under 6 dB/Bark above
    roughly 95 dB SPL), so the raw value is passed through a smooth algebraic
    saturation onto the open interval (6, 15). The saturation has unit slope at
    the midpoint of the interval, so wherever Terhardt's line already lies well
    inside the bounds it is reproduced to a fraction of a dB per Bark; outside
    them it bends over instead of being clipped, which keeps the result strictly
    monotone in level and differentiable everywhere.

    Returns a positive magnitude, strictly inside the two bound constants, and
    strictly decreasing in `maskerLevelDb`.
*/
double upwardSpreadSlopeDbPerBark (double maskerBark, double maskerLevelDb) noexcept;

/** Attenuation in dB of a masker's contribution at a probe position, relative to
    the masker's own level: the masking level at the probe is
    `maskerLevelDb + maskingSpreadDb (...)`.

    Two straight lines in Bark meeting at the masker, so the result is zero at
    zero distance, never positive, continuous in both Bark arguments and in the
    level, and monotone away from the masker in both directions.
*/
double maskingSpreadDb (double maskerBark, double probeBark, double maskerLevelDb) noexcept;

//==============================================================================
/** Per-bin masking threshold: for every bin, the highest masking level produced
    by all the *other* bins.

    Cost is O(numBins * log numBins) per call, and the result is the exact
    all-pairs maximum rather than an approximation over a bounded window. Two
    passes achieve that:

    - Descending (probes below a masker). The slope is a constant, so a masker's
      contribution splits as (level - slope*maskerBark) + slope*probeBark. The
      first bracket does not depend on the probe, so one running maximum taken
      from the top bin downwards covers every masker above every probe. O(n).

    - Ascending (probes above a masker). Here the slope depends on the masker's
      own level and frequency, so no such split exists; each masker contributes a
      different straight line in Bark and the threshold is the upper envelope of
      those lines. That envelope is maintained in a Li Chao tree keyed on the bin
      Bark values, which are fixed by prepare(). Inserting the masker for bin j
      only after querying bin j is what excludes a bin from masking itself.
      O(log n) per bin. A masker already below the envelope at its own position
      and decaying at least as fast can never overtake it and is not inserted at
      all, which is exact and, measured on random and on tonal spectra, cuts the
      cost by a factor of about four.

    A bounded window was rejected rather than not considered. The shallowest
    upward slope is 6 dB/Bark, so a 60 dB window is 10 Bark wide, and at the top
    of an 8192-point transform at 48 kHz 10 Bark is over two thousand bins - the
    window is not small where it matters, and the envelope method is both cheaper
    and exact.

    Levels in and out are dB on the same scale, `numBins` entries each, and the
    two buffers may not overlap. Bins that nothing masks report
    barkSpread::minimumLevelDb. Negative infinity is an acceptable input level.
    Before prepare() the object has no bins and compute() leaves its output
    untouched rather than writing anything.

    compute() allocates nothing; prepare() does all of it. It is const but writes
    to internal scratch, so, like Fft64, give every concurrent STFT path its own
    instance.
*/
class SpreadThreshold
{
public:
    /** @param numBinsIn  number of bins to process, normally fftSize / 2 + 1
        @param sampleRate sample rate in Hz, > 0
        @param fftSize    transform size the bins came from, > 0; bin k is
                          centred at k * sampleRate / fftSize
    */
    void prepare (int numBinsIn, double sampleRate, int fftSize);

    void compute (const double* binLevelsDb, double* thresholdDbOut) const noexcept;

    int getNumBins() const noexcept { return numBins; }

    /** Bark of each bin centre, as used internally. */
    const double* getBinBarks() const noexcept { return binBark.data(); }

private:
    //==============================================================================
    /** A masker's contribution as a straight line in Bark, f(z) = intercept +
        slope*z, so that the ascending pass is an upper-envelope problem. */
    struct Line
    {
        double intercept = 0.0;
        double slope = 0.0;

        /** The frame this node was last written in. A node stamped with anything but
            the current frame holds nothing — see the note on `generation`. */
        std::uint64_t generation = 0;

        double at (double z) const noexcept { return intercept + slope * z; }
    };

    /** The line at a node, or the empty line if it was last written in an earlier frame. */
    Line lineAt (int node) const noexcept;

    void insertLine (Line line) const noexcept;

    /** The line of the current envelope that is highest at a bin's Bark
        position. Returning the line rather than only its value lets compute()
        recognise, and skip inserting, a masker the envelope already dominates. */
    Line queryEnvelope (int bin) const noexcept;

    int numBins = 0;
    std::vector<double> binBark;        // Bark of each bin centre, strictly increasing
    std::vector<double> slopeBase;      // Terhardt raw slope magnitude at 0 dB bin level

    mutable std::vector<Line> envelope; // Li Chao tree over binBark

    /** Which frame the tree currently holds.

        The tree has to start each frame empty, and clearing it was doing that literally:
        4n nodes of 24 bytes, which for the long window is 131 KB memset per stage per
        frame — around 270 MB/s across a stereo instance in loop mode, enough to evict L2
        every hop. Since an insert only ever touches O(log n) nodes, almost all of that
        clearing was of nodes nobody would have read.

        So the counter is bumped instead and a node whose stamp does not match reads as
        empty. 64 bits, so it cannot wrap inside any session. */
    mutable std::uint64_t generation = 0;
};

} // namespace lime
