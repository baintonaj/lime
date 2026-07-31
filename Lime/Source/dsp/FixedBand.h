/*
  ==============================================================================

    FixedBand — one staggered stage's fixed band, per STFT bin.

    In the analog there is a single wideband fixed band per stage. Here every bin
    carries its own, which is the core of the per-bin reinterpretation.

    That change has a consequence worth stating, because it decides how modulation
    control is applied. MC3 and MC6 exist mainly to stop a *wideband* fixed band
    reducing its gain in response to energy outside its pass band. A per-bin fixed
    band only ever sees its own bin, so it is immune to that by construction. What
    remains necessary is the other half of modulation control, the level-dependent
    end-stop that section 2.3 describes: "after an appropriate degree of limiting
    has taken place at a given frequency, then it is unnecessary to continue the
    limiting when the signal level rises even further". That is applied per bin,
    via FixedBandParams::endStopFraction.

    Faithful features retained from Fig. 9 and Fig. 10:

      - Feedback control. Section 3.3 takes the control signal from the output of
        the variable gain element, not its input, which gives the soft 2:1 law the
        paper quotes for the intermediate level region.
      - Dual control paths into a maximum selector. The main path has the
        modulation control signal subtracted; the pass-band path does not. When
        modulation control over-suppresses the main path under complex signals, the
        pass-band path is "phased in, via the maximum selector circuit".
      - Two-stage smoothing: 15 ms on each control path, then 160 ms (HF) or
        300 ms (LF) after the selector.

  ==============================================================================
*/

#pragma once

#include <vector>

namespace lime
{

struct FixedBandParams
{
    /** Threshold in dB referred to the *system input*, not to this stage's input.
        Section 2.5 quotes approximately -30, -48 and -62 dB. */
    double thresholdDb = -30.0;

    /** Total low-level gain of the stages ahead of this one in the cascade. The
        analog achieves the same offset by raising "the ac and dc circuit gains
        for the mid- and low-level stages"; here it is explicit, so that a
        subthreshold signal receives full boost from every stage rather than
        driving the later ones over their own thresholds. */
    double precedingGainDb = 0.0;

    /** Low-level boost contributed by this stage. Section 1 gives "somewhat over
        8 dB" per stage, for 24 dB at high frequencies and 16 dB at low. */
    double lowLevelGainDb = 8.0;

    /** The per-bin end-stop: the fraction of the band's own input level at which
        the differential component stops falling.

        This is the surviving half of modulation control. Section 2.3: once enough
        limiting has happened at a given frequency, "the level of the side chain
        signal can be allowed to rise as a function of a further increase in signal
        level, whereby it stabilizes at some significant fraction of the main path
        signal level."

        Subtracting fraction * |x| from the control input produces exactly that.
        With |s| = g|x|, the effective control becomes |x|(g - fraction), so once g
        falls to `fraction` the control stops growing, the gain floors there, and
        |s| tracks |x| at a fixed ratio.

        Measured caveat, worth recording because it is easy to expect otherwise:
        this does *not* set a hard floor on the steady-state gain. The pass-band
        control path is deliberately not opposed by the end-stop — that is exactly
        what "phased in, via the maximum selector circuit" means — so as soon as the
        end-stop pulls the main path below the pass-band path, the selector switches
        over and the control keeps growing. Sweeping this parameter moves the
        steady-state curve by hundredths of a dB.

        Its real effect is on modulation immunity: it stops gain being reduced more
        than a given signal condition warrants, which shows up in the probe-tone
        behaviour of Figs. 14-16 rather than in a steady single-tone sweep. It is
        left enabled for that reason. Zero disables it. */
    double endStopFraction = 0.026;

    /** Weight of the feedforward term in the control signal, relative to the
        feedback term.

        Pure feedback control gives |s| = |x|/(1 + |s|/T), which tends to
        sqrt(T|x|) — the 2:1 law section 4.1 quotes for the intermediate region,
        but one that keeps *rising* with input, so the boost only decays as a
        square root and is still 3.7 dB at reference level. The 1967 paper says how
        the real characteristic is built: "such a characteristic is formed by
        deriving the compressor control voltage from a combination of feedforward
        and feedback signals". The feedforward term grows with the input itself,
        which limits the differential component to a roughly constant amplitude and
        collapses the boost to near unity above the action region — the "bilinear
        characteristic" with no action in the top 25 dB.

        The combination is normalised so the low-level knee stays at thresholdDb
        regardless of this weight.

        0.1 was chosen by measurement. This parameter trades the two figures the
        paper quotes for high-level behaviour against each other, and no value
        satisfies both:

            weight   total boost   at reference   steepest ratio
             0.0        24.00 dB      2.20 dB         1.67:1
             0.1        24.00 dB      0.91 dB         1.77:1
             0.5        24.00 dB      0.15 dB         1.89:1
             1.0        24.00 dB      0.06 dB         1.98:1

        Section 1 explicitly quantifies "a further dynamic action of about 1 dB
        above the reference level", which 0.1 reproduces at 0.91 dB. Section 4.1's
        "about 2:1" is the looser of the two statements, and 1.77:1 falls short of
        it because the thresholds are spaced 14 to 18 dB apart, so three stages of
        8 dB cannot compound more steeply than that without eliminating the 1 dB of
        high-level action. Matching the quantified figure was preferred to matching
        the approximate one. */
    double feedforwardWeight = 0.1;

    /** Exponent on the control ratio in g = 1 / (1 + (control/T)^n).

        n = 1 is the plain feedback law. Raising it sharpens the knee and steepens
        the rate at which boost is shed: above threshold the boost falls at roughly
        n dB per dB of input. Section 4.1 quotes "a compression ratio of about 2:1"
        for the intermediate level region, which n = 1 alone does not reach because
        the thresholds are spaced 14 to 18 dB apart and each stage can only shed its
        own 8 dB — measured, n = 1 peaks at 1.67:1.

        n = 2 was chosen by measurement: it gives a steepest local ratio of 1.99:1,
        holds the total boost at exactly 24.00 dB (HF) and 16.00 dB (LF), and pulls
        the onset of action up from about -80 dB towards the -62 dB the paper quotes.
        The cost is that the boost remaining at reference level falls to 0.06 dB
        rather than the "further dynamic action of about 1 dB" mentioned in
        section 1; sharpening and softening the knee trade these two figures against
        each other and cannot satisfy both. Near-unity at reference is the more
        important of the two, being the core dual-path premise that high-level
        signals travel the main path unprocessed, and section 4.1 attributes about
        1 dB of wideband compensation to the stage signal combiner rather than to the
        compressors, which plausibly accounts for the difference. */
    double kneeExponent = 2.0;

    double mainTauMs = 15.0;
    double passBandTauMs = 15.0;
    double finalTauMs = 160.0;

    /** Pass-band control weighting: a 1.6 kHz single-pole high-pass in the
        high-frequency stages, a 400 Hz single-pole low-pass in the low-frequency
        stages. */
    double passBandCornerHz = 1600.0;
    bool passBandIsHighPass = true;
};

//==============================================================================
class FixedBand
{
public:
    /** @param bandWeight  per-bin fraction of lowLevelGainDb this band contributes,
                             in [0,1]. The two frequency halves use power-complementary
                             single-pole weights, which sum to exactly one, so their dB
                             contributions add smoothly across the 800 Hz crossover.
                             Null means a uniform weight of one. */
    void prepare (const FixedBandParams& params, int numBins,
                  double sampleRate, int fftSize, double frameSeconds,
                  const double* bandWeight = nullptr);

    /** Recomputes the per-bin transfer scale — 10^(gain*weight/20) - 1 at full
        boost — for a new band weighting, which is what resampling the split onto a
        different number of sections amounts to. Allocation free, so it is safe to
        call from the audio thread: the vector is already sized by prepare(). */
    void setBandWeight (const double* bandWeight);

    const double* getTransferScale() const noexcept { return transferScale.data(); }

    void reset();

    int getNumBins() const noexcept { return numBins; }

    /** Advances one frame.

        @param magnitudes    this stage's input magnitude for each bin
        @param transferOut   per-bin side chain transfer F, >= 0, such that the
                             stage's output spectrum is input * (1 + F)
    */
    void processFrame (const double* magnitudes, double* transferOut);

    /** Smoothed per-bin control signal, for metering and the spectral display. */
    const double* getControl() const noexcept { return control.data(); }

private:
    FixedBandParams params;

    int numBins = 0;
    double thresholdLinear = 1.0;
    double controlNormalisation = 1.0;   // keeps the knee at thresholdDb

    std::vector<double> transferScale;   // per-bin F at full boost

    double mainCoeff = 0.0, passCoeff = 0.0, finalCoeff = 0.0;

    std::vector<double> passWeight;  // |H_passband(f_k)|
    std::vector<double> mainSmoothed, passSmoothed, control;
};

} // namespace lime
