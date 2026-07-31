/*
  ==============================================================================

    SlidingLayer — the sliding-band half of each staggered stage.

    In the analog, a stage's sliding band is one variable filter whose corner
    slides with level, maintaining noise reduction above a dominant component in
    the high-frequency stages and below it in the low-frequency stages. Its
    partner, the fixed band, supplies uniform but reduced boost on the other side.
    Action substitution combines them: "the action of one circuit augments that of
    the other, and, as the occasion requires, may be substituted for that of the
    other."

    The per-bin reinterpretation changes what is left for this layer to do. A
    per-bin fixed band already boosts every bin except a dominant's own, which is
    *better* than the analog by the paper's own least-treatment criterion. What it
    cannot do is anything across bins: it produces a one-bin notch at each
    dominant, where the analog produces a smooth skirt. Section 2.2 says why that
    matters — fixed-band behaviour is preferred on the stop-band side precisely so
    that a gain or frequency response error in the channel causes "no undue
    exaggeration of the error at other, nondominant signal frequencies", which a
    razor-thin notch cannot deliver.

    So this layer supplies the cross-bin profile, shaped by the masking spread
    function. That is not a liberty: the paper states the same thing in
    psychoacoustic terms — "The only region of the spectrum that is not boosted in
    gain is the region that is controlled by masking, where audible low-level
    information does not exist."

    Two compositions are provided because they are genuinely different designs and
    the choice is being settled by measurement rather than argument. See
    SlidingComposition.

  ==============================================================================
*/

#pragma once

#include "Bark.h"

#include <vector>

namespace lime
{

enum class SlidingComposition
{
    /** phi_stage = phi_fixed * sigma_slide.

        The layer withholds boost across the masking skirt around each dominant,
        multiplying the fixed band's per-bin result. Applying the paper's
        augmenting law to a withholding operation would restore boost exactly where
        the skirt meant to remove it, so restriction is expressed as a product.
        Several skirts still combine with each other by 1 - prod(1 - F_i), where
        each genuinely is a provider. Buys the smooth composite filter and the
        channel-error tolerance of section 2.2, at the cost of boosting less than
        the theoretical maximum. */
    restriction,

    /** phi_stage = 1 - (1 - phi_fixed)(1 - phi_sliding), Eq. (1) applied literally,
        with *both* circuits driven by the spread level so that each sees the masked
        region as occupied and compresses there.

        Closer to the stated intent that the compressor "should try to keep all
        signal components fully boosted at all times", and it preserves the paper's
        own composition law. The reduction profile is thinner, so less of section
        2.2's tolerance is bought. */
    actionSubstitution
};

//==============================================================================
struct SlidingLayerParams
{
    /** High-frequency stages slide upward and maintain action above a dominant;
        low-frequency stages slide downward and maintain it below. */
    bool highFrequency = true;

    /** Matches the partner fixed band, so the two are on the same scale. */
    double lowLevelGainDb = 8.0;
    double thresholdDb = -30.0;
    double precedingGainDb = 0.0;
    double kneeExponent = 2.0;

    /** Control weighting corner: about 10 kHz single-pole high-pass in the
        high-frequency stages (section 3.3, "different for each stage"), 80 Hz
        single-pole low-pass in the low-frequency stages (section 3.4). */
    double controlCornerHz = 10000.0;

    /** Two-stage smoothing: 5 ms then 80 ms at high frequencies, 7.5 ms then
        150 ms at low. Deliberately slower than the fixed band, which the paper
        says the fixed band's fast frequency-independent recovery permits. */
    double firstTauMs = 5.0;
    double finalTauMs = 80.0;

    /** End-stop on how much boost the layer may withhold.

        Section 2.3: once "an overall gain of about unity is obtained — when the
        variable filter cutoff frequency is about two to three octaves above the
        dominant signal frequency — there is no reason for further sliding of the
        variable filter." Capping the withheld fraction is the direct expression of
        that, and it is what stops a loud dominant from stripping boost across the
        whole spectrum. */
    double maximumWithheldFraction = 0.85;

    /** How strongly the modulation control level opposes the layer's own control.

        The MC level arrives normalised to the per-bin magnitude scale (see
        ModulationControl) and is added onto the compression-law reference, so the
        ratio the knee sees is control against threshold-plus-MC — section 3.3's
        end-stop expressed as the ratio it describes. At zero the opposition is
        disabled and the layer runs on its bare threshold. */
    double mcOppositionGain = 1.0;

    SlidingComposition composition = SlidingComposition::restriction;
};

//==============================================================================
class SlidingLayer
{
public:
    void prepare (const SlidingLayerParams& params, int numBins,
                  double sampleRate, int fftSize, double frameSeconds);

    void reset();

    int getNumBins() const noexcept { return numBins; }
    SlidingComposition getComposition() const noexcept { return params.composition; }

    /** MC1 for the high-frequency stages, MC4 for the low, on the normalised
        per-bin magnitude scale ModulationControl produces. Fed in opposition to
        the control signal, as section 3.2 specifies: it raises the reference the
        compression law measures against, so the ratio between the rectified
        control and the MC level forms the sliding end-stop. Smoothed inside the
        frame loop by the same two poles as the control itself, so the opposition
        is as well-behaved as what it opposes. */
    void setModulationControl (double level) noexcept { modulationControl = std::max (0.0, level); }

    /** The user's attack/release override on the final smoothers — the 80/150 ms
        poles the control and the MC opposition share. Both take the pair
        together: the opposition has to stay as smooth as what it opposes (see
        processFrame), so letting the per-bin control run user times while the
        MC reference kept the published pole would transiently mis-form the
        end-stop ratio. Engaged, each smoother runs the attack time while its
        input is rising and the release time while it falls; zero for either
        value falls back to that direction's published tau, and both zero is
        bit-identical to the published symmetric smoother. The 5/7.5 ms first
        poles keep published values. Allocation free. */
    void setTimeConstants (double attackMs, double releaseMs);

    /** Advances one frame.

        @param magnitudes    this stage's input magnitude per bin
        @param fractionOut   in restriction mode, the surviving boost fraction in
                             [0,1] to multiply the fixed band by; in action
                             substitution mode, this layer's own boost fraction
        @param spreadLevelOut  the per-bin spread level in linear magnitude, which
                               action substitution mode also needs to drive the
                               fixed band. May be null.
    */
    void processFrame (const double* magnitudes, double* fractionOut,
                       double* spreadLevelOut = nullptr);

    /** Number of bins selected as dominants in the last frame, for metering. */
    int getNumDominants() const noexcept { return numDominants; }

    /** Per-bin masking threshold from the last frame, in dB, for the display. */
    const double* getSpreadThresholdDb() const noexcept { return thresholdDb.data(); }

private:
    SlidingLayerParams params;

    int numBins = 0;
    int numDominants = 0;

    double thresholdLinear = 1.0;
    double firstCoeff = 0.0, finalCoeff = 0.0;
    double modulationControl = 0.0;

    // The final smoothers' asymmetric pair. Both sit at finalCoeff until a user
    // time engages, so the default path is the published symmetric smoother.
    double finalRiseCoeff = 0.0, finalFallCoeff = 0.0;
    double currentFrameSeconds = 0.0;

    // The MC level smoothed by the same two poles as the per-bin control; see
    // processFrame.
    double mcFirstSmoothed = 0.0, mcSmoothed = 0.0;

    SpreadThreshold spread;

    std::vector<double> controlWeight;   // |H(f_k)| of the control weighting
    std::vector<double> levelDb, thresholdDb, spreadLinear;
    std::vector<double> firstSmoothed, control;
};

} // namespace lime
