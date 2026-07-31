/*
  ==============================================================================

    SrSideChain — the staggered stage cascade for one frequency half.

    Three stages above the 800 Hz crossover, two below. Section 2.5: "In the
    series-connected staggered action format there is a compounding of the actions
    of the individual stages; the transfer functions of the several stages are
    multiplied, whereby the dB characteristics add." The stages are cascaded, each
    seeing the previous one's output, and the half's total transfer is

        total = prod (1 + F_i)

    which is applied to the spectrum directly on encode and divided out on decode.

    ## Why the two halves multiply rather than sum

    The plan's first sketch summed the two halves' side chain contributions into the
    main path. That cannot be inverted exactly. Summing gives v = (I + H + L)u with H
    and L computed at different analysis resolutions, and there is no closed form for
    the inverse of a sum of two different-resolution filters; iterating diverges,
    because the side chain contribution reaches roughly fifteen times the signal at
    full boost.

    Cascading multiplicatively makes each half exactly invertible by per-bin division
    in its own resolution, and it is what section 2.5 describes anyway. The cost is
    that the halves' boosts add in dB, so each is shaped by a per-bin band weight
    chosen so the sum lands where Fig. 13 puts it: the underlying profile is the
    power-complementary single-pole pair at 800 Hz, |HP|^2 + |LP|^2 = 1 exactly, so
    the total runs smoothly from 16 dB below the crossover through 20 dB at it to
    24 dB above.

    ## Why the split is a resolution rather than a frequency

    The analog machine has exactly two band-defining networks and meets them at 800 Hz
    (reference hardware manual, section 5.6 and Fig. 5.10). That is a limit of what can be built
    out of op-amps, not of the process, and a transform is not subject to it: the
    allocation between the two halves can be resolved as finely as the bins allow.

    So the profile above is not evaluated per bin directly. It is sampled onto a number
    of log-spaced sections, and that count is the control — a power of two from 4 up to
    one section per bin, shown on the panel as its exponent. The corner itself is a
    constant again, at the value the reference system uses.

    Two things about that are load-bearing rather than incidental:

      - The section grid is defined in **hertz**, not in bins. The two halves run at
        different resolutions, so a grid in bins would quantise the profile differently
        in each and |HP|^2 + |LP|^2 would stop being 1. In hertz both halves evaluate
        the identical function of frequency and the sum stays exact.

      - For the same reason the count is **not clamped against a half's own bin count**.
        Sections finer than a half can resolve cost nothing and change nothing; clamping
        them per half would break the complementarity the whole cascade rests on.

    ## Why the control state is one frame behind

    Every stage derives its control from the *previous* frame's magnitudes. That
    single choice is what makes the decoder exact rather than approximate: encoder
    and decoder both build their control state from the same |u| history, so they
    compute bit-identical transfers and the division undoes the multiplication. The
    cost is one frame of control lag, 5.33 ms at 48 kHz, against time constants of
    15 to 300 ms.

  ==============================================================================
*/

#pragma once

#include "Bark.h"
#include "GainSlew.h"
#include "ShelvingBank64.h"
#include "SpectrumProbe.h"
#include "ModulationControl.h"
#include "SrStage.h"
#include "Wola64.h"

#include <complex>
#include <memory>
#include <vector>

namespace lime
{

/** How the computed gain profile is actually applied to the signal. */
enum class GainRealisation
{
    /** Multiply the bins, divide them out on decode. Simple, and it reproduces every
        published curve, but it is only approximately invertible: a modified spectrum
        with abrupt bin-to-bin gain changes is not the transform of any signal. */
    stft,

    /** EXPERIMENTAL, and currently worse than `stft`. Do not select it expecting an
        improvement.

        Use the transform for analysis only and apply the profile through a Bark-spaced
        bank of direct-form-I peaking biquads. The bank itself is exact: measured in
        isolation it inverts to -235 dB with a static profile and -234 dB with its gains
        randomised every block, against -41 dB for the STFT. Integrating it is the part
        that is unsolved, and three attempts are recorded here because each failed for a
        different and instructive reason:

          1. Encoder analysing its input, decoder analysing its own input. The two build
             different state, so nothing cancels. Measured +1 to -18 dB.
          2. Both analysing the *encoded* signal, signal delayed to line up with the
             analysis. Better, -31 to -41 dB, but still short: the decoder undoes a block
             that was encoded with coefficients from one window earlier while using its
             own current ones, so the coefficient sequences are offset by the analysis
             latency.
          3. Both analysing the encoded signal, no signal delay, so both use equally
             stale coefficients. This diverges outright, +6 dB: with no delay between the
             applied gain and the level that gain is derived from, the loop runs away —
             gain up raises the analysed level, which pulls the gain down.

        What the three together say is that the coefficient sequences must match sample
        for sample *and* the loop from applied gain back to measured level must be broken.
        Satisfying both needs the control to be derived from something neither side has to
        infer — most likely a level estimate transmitted or reconstructed independently of
        the gain, or per-band envelope followers outside the gain loop, as the analog has.

        The trade it would buy is resolution of the *applied* gain: 24 bands rather than
        one per bin, with analysis staying per-bin. See ShelvingBank64 and DESIGN.md. */
    filterBank
};

enum class SrDirection
{
    encode,   // multiply the spectrum by the half's transfer
    decode    // divide it out
};

//==============================================================================
class SrSideChain
{
public:
    /** Where the underlying power-complementary pair meets, which the reference system
        fixes at 800 Hz and this no longer exposes. See the note at the top of the file. */
    static constexpr double referenceCrossoverHz = 800.0;

    /** The band the section grid spans, and the range of the count.

        20 Hz to 20 kHz because that is the band the profile has anything to say about, and
        because it is what the panel's plot already draws. Bins outside it take the outermost
        section, which costs nothing: the pole is 0.0006 at 20 Hz and 0.9984 at 20 kHz, so
        both ends are already within a hundredth of a decibel of their asymptotes.

        The count runs 4 to 16384 in powers of two, so it is exactly 2^e for an exponent e of
        2 to 14, which is what the panel exposes. Fineness is close to free and coarseness is
        not: a staircase spanning 8 dB in N steps has 8/N dB edges, and a gain profile that
        jumps between neighbouring bins is what setEnvelopeSmoothingBark measured as
        uninvertible.

        The floor is 4 rather than 1, and what it excludes is worth recording because both
        exclusions were deliberate and neither was obvious:

          - **1 section** is a single weight across the whole band. It has no edges anywhere,
            so it inverts better than any other setting — measured at -67.5 dB against the
            fine grid's -46.6 — but it has no frequency-dependent allocation either, so there
            is nothing left of the process it belongs to. A control position that nulls best
            by not working is a trap, not an option.

          - **2 sections** is the analog machine's own arrangement: one crossover at 800 Hz,
            three stages above and two below. Dropping it costs the ability to reproduce the
            reference system exactly, which in a project about fidelity is a real loss and is
            not pretended otherwise. It goes because it is also the worst setting in the
            range by a wide margin — the entire 8 dB arrives as one step, and it nulls at
            -24.0 dB, 22 dB behind everything above it.

        So the range starts where the control starts being useful rather than where the
        arithmetic starts. Where one section per bin lands depends on the rate, because the
        long-window half's bin count does:

            rate        long window   bins    one section per bin at
            44.1/48 kHz     4096      2049            2048
            88.2/96 kHz     8192      4097            4096
            176.4/192 kHz  16384      8193            8192

        The top of the range is 16384, which is past per-bin at every rate the plugin runs
        at. It is there to make the position count odd, so that the default sits straight up
        on the knob — the same reason the old crossover's taper was skewed to put 800 Hz at
        twelve o'clock. That is a panel reason rather than a process one and is recorded as
        such; the position is inert, for the reason the next paragraph gives, so the cost of
        having it is a few kilobytes of table and nothing else.

        So at 48 kHz the last two positions are asking for a grid finer than the transform
        has bins to put it on. That costs nothing and is deliberately not clamped: several
        sections landing in one bin simply means the bin takes the weight of whichever
        contains it, and the profile is a function of frequency either way. Clamping the
        count to a half's own bin count is the one thing that would break it — the two
        halves would quantise differently and |HP|^2 + |LP|^2 would stop being 1. See the
        note at the top of the file.

        The panel prints the count it was asked for rather than the count the rate can
        resolve, which is the honest way round: the setting has to mean the same thing in
        every session, because an encode and its decode must agree on it. */
    static constexpr double sectionLoHz = 20.0;
    static constexpr double sectionHiHz = 20000.0;

    static constexpr int minSections = 4;
    static constexpr int maxSections = 16384;

    /** The exponent the panel actually shows, so that the mapping lives in one place rather
        than being reinvented by the parameter layout. sections == 1 << exponent. */
    static constexpr int minExponent = 2;    // 4 sections
    static constexpr int maxExponent = 14;   // 16384 sections

    /** 256, which is not the top of the range and is chosen rather than conceded.

        The null plateaus at about 64 sections — past that the 8/N dB steps the staircase
        puts into the profile are already below whatever else limits the round trip, and
        every finer setting measures the same. Measured at 48 kHz, 256 nulls at -46.6 dB,
        which is the finest grid's figure to the tenth, so nothing is given up on the round
        trip by not going further.

        It also lands on the continuous profile the count replaced: 18.68 dB at 500 Hz
        against the finest grid's 18.67, a hundredth of a decibel out. That is worth
        recording because it is better than it had to be. The knee argument alone only
        justifies 128, which is the first power of two clear of it, and 128 reads 18.83 dB
        there — 0.16 dB out, inaudible but real. Doubling to 256 removes even that.

        What settles it between the two is the panel rather than the process. An odd number
        of detents puts the default straight up on the knob, the range runs 2 to 14 in the
        exponent, and 8 is its centre. So the arithmetic that makes the pointer vertical also
        happens to pick the setting that traces the reference characteristic most closely,
        which is a coincidence and is flagged as one — if the range moves, this reasoning has
        to be redone rather than assumed to still hold. */
    static constexpr int defaultSections = 256;

    /** The scratch the filter-bank realisation analyses through is sized here, so
        process() never allocates. Callers that know their host's block size pass it;
        the default covers any block a standalone test is likely to use. */
    static constexpr int defaultMaxBlockSamples = 8192;

    /** @param highFrequency  true for the three-stage half above 800 Hz
        @param fftOrder       window is 1 << fftOrder
        @param slewFrameScale the gain slew's limits are quoted per frame (see
                              GainSlewParams), so a half running a longer hop has to be
                              allowed proportionally more per frame or it would be
                              limiting four times as hard for the same rate of change.
                              1.0 for the half the figures were quoted for.
        @param maxBlockSamples the largest block process() will be handed, so every
                              scratch buffer is sized once here rather than grown on
                              the audio thread.
    */
    void prepare (bool highFrequency, int fftOrder, int hop, double sampleRate,
                  SlidingComposition composition = SlidingComposition::restriction,
                  double slewFrameScale = 1.0,
                  int maxBlockSamples = defaultMaxBlockSamples);

    void reset();

    void setDirection (SrDirection d) noexcept { direction = d; }

    /** Sets how finely the band-defining split between the two halves is resolved.

        Rounded down to a power of two and clamped to [minSections, maxSections], so a
        caller can hand it anything. At 2 the process is the analog machine's: one crossover
        at 800 Hz, three stages above it and two below. At the default it is one section per
        bin and the allocation is effectively continuous.

        It changes the encode characteristic, so an encode and its decode have to agree on
        it — the same discipline as tape alignment. In Both mode they agree by construction.

        Allocation free and cheap enough to call from the audio thread, and a no-op unless
        the value actually moved, so a static control costs nothing. */
    void setSections (int count);

    int getSections() const noexcept { return sections; }

    /** High-passes what the compressors' level detection reads; the signal path
        never passes through it.

        Third-order Butterworth magnitude, sampled per bin onto the detection
        levels ahead of the stages. Below the corner a bin reads quieter than it
        is, so the compressors treat it as low-level programme; above the corner
        detection is untouched. Zero or negative disables the weighting, which is
        the prepared default.

        Like setSections it changes the encode characteristic, so an encode and
        its decode have to agree on it — and for the same reason it is pushed to
        all four side chains by the engine, never set on one alone. Allocation
        free and a no-op unless the corner actually moved. */
    void setDetectorHighPass (double cornerHz);

    double getDetectorHighPass() const noexcept { return detectorHighPassHz; }

    /** The matching low pass on the detection: above the corner a bin reads
        quieter than it is. The two corners multiply into one per-bin weight, so
        together they band-limit what the compressors hear without either knowing
        about the other. Same contract as setDetectorHighPass in every respect. */
    void setDetectorLowPass (double cornerHz);

    double getDetectorLowPass() const noexcept { return detectorLowPassHz; }

    void setGainRealisation (GainRealisation r) noexcept { realisation = r; }
    GainRealisation getGainRealisation() const noexcept { return realisation; }

    /** Smooths the composite per-bin transfer across frequency, in bins.

        Built to test the hypothesis that the single-bin notch a per-bin fixed band
        creates at each dominant was what limited the encode/decode null, on the
        reasoning that sharp spectral gain implies a long impulse response. Measured,
        that is wrong, and the numbers are kept here because the result is useful:

            smoothing (bins)   broadband null   1 kHz tone null
                   0               -48.2 dB         -77.0 dB
                   1               -49.7            -58.2
                   4               -50.3            -33.2
                  16               -51.3            -21.4
                  32               -52.0            -18.4

        Smoothing barely helps broadband material and severely harms tonal material.
        Smoothing a deep single-bin notch in dB spreads it into a wide, deep trough,
        which is *more* total spectral modification than the notch it replaced, so
        consistency gets worse rather than better. Sharpness was not the problem;
        depth and extent are.

        Left in place, defaulting to zero, because it is a genuine control over the
        composite filter's smoothness and section 2.2 does care about that for
        channel-error tolerance — just not for the null. */
    void setSpectralSmoothingBins (double bins) noexcept { smoothingBins = bins; }

    /** Width, in Bark, of the critical-band smoothing applied to the level envelope
        that drives the stage controls.

        This is what makes the process invertible, and it is not a cosmetic choice.
        Measured on a bare pair of WOLAs with a time-invariant gain, so that control
        state and divergence are excluded entirely:

            gain profile              depth    residual
            smooth tilt across band   24 dB    -125.9 dB
            single-bin notch          24 dB     -15.9 dB
            notch every 32 bins       24 dB      +4.0 dB

        A per-bin fixed band notches single bins, and a notched spectrum is simply not
        the transform of any signal — no decoder can undo it. A gradual profile is
        undone to -126 dB. So the per-bin *units* are fine; what cannot survive is a
        gain that jumps tens of dB between neighbouring bins.

        Smoothing the level envelope before the control sees it fixes this at the
        source, and keeps every bin's own threshold, time constants and weighting.
        Smoothing the finished gain instead does not work — see
        setSpectralSmoothingBins. */
    void setEnvelopeSmoothingBark (double bark)
    {
        envelopeBark = bark;
        rebuildEnvelopeCoefficients();
    }

    double getEnvelopeSmoothingBark() const noexcept { return envelopeBark; }

    /** Enables the per-frame slew limit on the composite transfer, which stands in
        for the analog overshoot suppressors. Applied identically in both directions,
        so it does not disturb complementarity. */
    void setGainSlewEnabled (bool enabled) noexcept { slewEnabled = enabled; }
    bool isGainSlewEnabled() const noexcept { return slewEnabled; }

    /** Floor on the envelope smoothing width, in bins.

        A Bark width alone is not enough in the short-window half. At 1024 points and
        48 kHz the bins are 46.9 Hz apart, so one Bark near 100 Hz spans barely two
        bins and smoothing over it does almost nothing — which is why the measured
        null in that half was indifferent to the Bark setting while improving steadily
        with window length. Imposing a minimum width in bins smooths where the Bark
        scale cannot.

        Measured, this does not help, and the negative result is worth recording. The
        short-window half sits at -31.2 dB for every width from 0 to 64 bins, while the
        long-window half improves steadily (-63.5 to -67.8 dB). So profile sharpness is
        not what limits the short-window half — a third mechanism is, and it is not yet
        identified. Defaults to zero. See DESIGN.md. */
    void setMinimumSmoothingBins (double bins)
    {
        minimumSmoothingBins = bins;
        rebuildEnvelopeCoefficients();
    }

    double getMinimumSmoothingBins() const noexcept { return minimumSmoothingBins; }
    double getSpectralSmoothingBins() const noexcept { return smoothingBins; }
    SrDirection getDirection() const noexcept { return direction; }

    /** Signal-path latency of this half under the current realisation.

        The STFT realisation delays the signal by one analysis window. The
        filter-bank realisation applies its gain in the time domain and uses the
        transform for analysis only, so its signal path has no delay at all —
        deliberately, see the note in process() — and reporting the window here
        anyway would make a host compensate for a delay that does not exist. */
    int getLatencySamples() const noexcept
    {
        if (realisation == GainRealisation::filterBank)
            return 0;

        return wola != nullptr ? wola->getLatencySamples() : 0;
    }
    int getNumStages() const noexcept { return (int) stages.size(); }
    int getNumBins() const noexcept { return wola != nullptr ? wola->getNumBins() : 0; }

    /** Dominants seen across the stages in the last frame, for metering. */
    int getNumDominants() const noexcept { return lastDominants; }

    /** Total boost in dB this half would apply at each bin under full low-level
        conditions, for the display and for tests. */
    const double* getLowLevelBoostDb() const noexcept { return lowLevelBoostDb.data(); }

    /** Hands this half a probe to publish its input level and composite transfer
        into, for the editor's display. The probe is owned by the caller and must
        outlive this side chain — the engine keeps two engine-lifetime probes for
        exactly that reason, so the editor's pointers stay valid across a
        re-prepare. Null (the default) publishes nothing. */
    void setProbe (SpectrumProbe* p) noexcept { probe = p; }

    /** Streams a block through the half. Input and output may alias. */
    void process (const double* input, double* output, int numSamples);

private:
    int sections = defaultSections;

    void processFrame (std::complex<double>* bins, int numBins);

    /** One scratch-sized piece of the filter-bank realisation; process() splits
        oversized blocks so this never outgrows the buffers prepare() sized. */
    void processFilterBankChunk (const double* input, double* output, int numSamples);

    /** Fills bandWeight and lowLevelBoostDb from the current section count and pushes the
        weighting into the stages. Allocation free once prepare() has sized the vectors. */
    void rebuildBandWeights();

    /** Fills detectorWeight from the current detector high-pass corner. Allocation free
        once prepare() has sized the vector. */
    void rebuildDetectorWeights();

    /** Rebuilds the per-bin envelope smoothing coefficients. Cheap, and only reached from
        prepare() and the two setters that can change them, so nothing in the frame loop
        has to work out what they are. */
    void rebuildEnvelopeCoefficients();

    std::unique_ptr<Wola64> wola;
    std::vector<SrStage> stages;
    ModulationControl modulation;
    GainSlew slew;
    ShelvingBank64 bank;
    SpectrumProbe* probe = nullptr;

    bool isHighFrequency = true;
    SrDirection direction = SrDirection::encode;
    GainRealisation realisation = GainRealisation::stft;
    int lastDominants = 0;

    std::vector<double> bandWeight, lowLevelBoostDb;
    std::vector<double> binHz;   // bin centre frequencies, so the weights can be rebuilt

    // The user's side-chain filters, combined into one per-bin weight on the
    // detection levels.
    std::vector<double> detectorWeight;
    double detectorHighPassHz = 0.0;
    double detectorLowPassHz = 0.0;

    // Each bin's position along the log section grid, in [0,1], and the weight of each
    // section. Split this way so that rebuilding costs one pow per section rather than one
    // per bin: the positions depend only on the sample rate and the window, so they are
    // built once, and the section table is at most 2048 long however many bins there are.
    std::vector<double> binLogPosition, sectionWeight;
    std::vector<double> magnitudes, running, transfer, total, mcFeed, pending;
    std::vector<double> bandTargetDb, analyserScratch, discardScratch;
    std::vector<int> bandBin;

    /** MC comes from the second-stage adder output, so like the stage control it is
        one frame behind. Both directions therefore see the same value. */
    double slidingMc = 0.0;
    double smoothingBins = 0.0;
    double magnitudeScale = 1.0;   // bin magnitude -> amplitude
    double currentSampleRate = 48000.0;
    double envelopeBark = 1.0;
    double minimumSmoothingBins = 0.0;   // measured not to help; see the setter
    bool slewEnabled = true;
    bool haveState = false;

    std::vector<double> envelope, forwardCoeff;

    // The envelope smoother's per-bin coefficient, its complement, and the peak
    // normalisation that goes with it. Tables rather than expressions: see
    // rebuildEnvelopeCoefficients.
    std::vector<double> envelopeCoeff, envelopeOneMinusCoeff, envelopePeakNormalisation;
};

} // namespace lime
