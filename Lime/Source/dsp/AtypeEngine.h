/*
  ==============================================================================

    AtypeEngine — the four-band process noise reduction, four bands, time domain.

    A-type ("An Audio Noise Reduction System", JAES 15, 1967) is the
    ancestor of SR and is offered here as the switchable alternative the Cat. No.
    300 module's SR/Off/A selector implies. It is built directly in the time
    domain rather than on the STFT core: A-type is a four-band analog design with
    wide, gentle bands and per-band time constants of tens of milliseconds, so
    biquads and sample recursions realise it more faithfully — and with no
    latency — than any transform would.

    Dual-path differential topology (section 1 of the paper). The main path
    passes the signal unaltered; a side chain produces a small differential
    component, added on encode and subtracted on decode:

        encode:  y = (1 + G(x)) * x
        decode:  z = y - z * G(z)    =>    z = y / (1 + G(z))

    The decoder is a feedback configuration, and that is what makes the pair
    exactly complementary: z = x whenever the two side chains are identical. Here
    the encoder's sample map is continuous, piecewise linear and strictly
    increasing in the current sample, so the decoder inverts it in closed form
    rather than approximately — see solveForInput() in the implementation. The
    null is therefore limited by floating point, not by the control smoothing.

    Bands (section 2 and Fig. 4): 80 Hz low-pass, 80 Hz to 3 kHz band-pass,
    3 kHz high-pass, 9 kHz high-pass. Bands 1, 3 and 4 are conventional 12 dB per
    octave filters; band 2 is the exact complement of bands 1 and 3, so the three
    sum to unity at every frequency and phase. Bands 3 and 4 overlap, which is
    what makes them contribute roughly equally above 5 kHz.

    Channels are processed independently. The paper notes that the topology's
    tolerance of gain error "enables the noise reduction system to operate
    without control signal interconnections", so there is no stereo linking.

    No JUCE, no float, and no allocation after prepare().

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace lime
{

enum class AtypeMode
{
    encode,   // differential component added: the signal you would print to tape
    decode    // differential component subtracted, feedback topology
};

//==============================================================================
/** The band-defining filters and the low-level gain target they are fitted to.

    Mirrors srNetworks in FilterDesign.h: the published constants, the design
    functions, and the acceptance curve all live together so a test can reach
    them without instantiating an engine.
*/
namespace atypeBands
{
    // Section 2, stated verbatim.
    inline constexpr double band1CornerHz = 80.0;     // low-pass
    inline constexpr double band3CornerHz = 3000.0;   // high-pass
    inline constexpr double band4CornerHz = 9000.0;   // high-pass
    inline constexpr double butterworthQ = 0.70710678118654752440;

    // "The compression thresholds are 40 dB below peak operating level."
    inline constexpr double thresholdDb = -40.0;

    // "The clipping level is chosen so that the overshoot is limited to 2 dB for
    // peak-amplitude step inputs."
    inline constexpr double overshootLimitDb = 2.0;

    // Section 3. The steady-state attack is deliberately slow so that gain
    // modulation stays inaudible; it is shortened in proportion to the size of
    // the amplitude transition, down to the fast limit below. Recovery has to
    // restore the noise reduction action within about 0.1 s, the interval over
    // which residual auditory masking holds — that is a requirement on the
    // completion of recovery, so the time constant is a third of it.
    inline constexpr double attackSeconds = 0.1;
    inline constexpr double fastestAttackSeconds = 0.001;
    inline constexpr double recoverySeconds = 0.0333333333333333333;

    // Control is taken from a combination of peak and average value, which the
    // paper offers as a sufficient proxy for rms and which keeps the system
    // tracking when the channel has badly non-linear phase. The average is
    // scaled by pi/2 so that a sine reads its own amplitude in both halves, and
    // the two are then mixed equally.
    inline constexpr double peakAttackSeconds = 0.0005;
    inline constexpr double peakReleaseSeconds = 0.1;
    inline constexpr double averageSeconds = 0.03;
    inline constexpr double peakWeight = 0.5;
    inline constexpr double averageScale = 1.57079632679489661923;   // pi/2

    // The low-level acceptance curve: "an output which is uniformly 10 dB above
    // the input up to about 5 kHz, rising smoothly to 15 dB at 15 kHz".
    inline constexpr double flatGainDb = 10.0;
    inline constexpr double flatUpToHz = 5000.0;
    inline constexpr double topGainDb = 15.0;
    inline constexpr double topHz = 15000.0;

    /** 12 dB/octave Butterworth low-pass, bilinear with the corner prewarped so
        it lands on `cornerHz`. Corners are clamped just below Nyquist so the
        design stays valid at sample rates that do not reach 9 kHz. */
    BiquadCoeffs lowPass12 (double cornerHz, double sampleRate);

    /** 12 dB/octave Butterworth high-pass, same conventions. */
    BiquadCoeffs highPass12 (double cornerHz, double sampleRate);

    /** The target of the low-level gain fit: flatGainDb up to flatUpToHz, then a
        smooth (log-linear in dB) rise to topGainDb at topHz, flat above. */
    double targetGainDb (double hz);

    /** Solves the four side-chain band gains so that the summed low-level
        response, 1 + sum g_i B_i, reproduces targetGainDb().

        Four scalars are the only freedom taken: the band corners and slopes are
        the paper's. The fit has to be numerical because the bands sum as complex
        responses — band 4's 12 dB/octave high-pass is 180 degrees out below its
        corner, so its contribution subtracts there and no closed-form choice of
        gains makes the total both flat and 5 dB up at 15 kHz. Deterministic
        Nelder-Mead, run once per sample rate change, exactly as
        srNetworks::mainPathBiquadsForRate solves for antisaturation.

        Band order is the paper's: index 0 is band 1.
    */
    std::array<double, 4> solveBandGains (double sampleRate,
                                          const BiquadCoeffs& band1,
                                          const BiquadCoeffs& band3,
                                          const BiquadCoeffs& band4);
}

//==============================================================================
class AtypeEngine
{
public:
    static constexpr int numBands = 4;

    AtypeEngine();
    ~AtypeEngine();

    AtypeEngine (const AtypeEngine&) = delete;
    AtypeEngine& operator= (const AtypeEngine&) = delete;

    /** `maxBlockSize` is accepted only so this matches the other engines'
        prepare(); the engine is a sample recursion and keeps no scratch. */
    void prepare (double sampleRate, int maxBlockSize, int numChannels);
    void reset();

    void setMode (AtypeMode m) noexcept { mode = m; }
    AtypeMode getMode() const noexcept { return mode; }

    /** Peak operating level as an amplitude; thresholds sit thresholdDb below
        it. Full scale by default. Takes effect on the next prepare(). */
    void setPeakOperatingLevel (double amplitude) noexcept;
    double getPeakOperatingLevel() const noexcept { return peakLevel; }

    /** The user's side-chain high pass: high-passes what the four compressors'
        level detection reads, without touching the signal path or the band
        filters. Each band's detector input is scaled by the third-order
        Butterworth magnitude at the band's geometric centre — the per-band
        equivalent of the per-bin weighting the spectral process applies. The
        encoder and the decoder feed their detectors the same band signals, so an
        identical corner on both leaves the round trip exact. Zero or negative
        disables it, which is the default. Allocation free and a no-op unless the
        corner moved. */
    void setSidechainHighPass (double cornerHz) noexcept;
    double getSidechainHighPass() const noexcept { return detectorHighPassHz; }

    /** The matching low pass on the detection: the two corners multiply into one
        weight per band, so together they band-limit what the compressors hear.
        Same contract as setSidechainHighPass in every respect. */
    void setSidechainLowPass (double cornerHz) noexcept;
    double getSidechainLowPass() const noexcept { return detectorLowPassHz; }

    /** The user's attack/release times, one pair per direction.

        The up pair replaces the steady-state attack rate and the recovery time
        while the engine encodes; the down pair does the same while it decodes.
        The peak and average detector constants and the 1 ms fast-attack cap
        stay published — a user attack below 1 ms only steepens the
        transition-proportional scaling, it never beats the cap. Zero for any
        value keeps that direction's published constant, and all four zero is
        bit-identical to the published ballistics. The round trip is exact only
        while an encoder and its decoder run the same pair; setting them apart
        is a creative choice that forfeits the exact null. Allocation free and
        a no-op unless a time actually moved. */
    void setTimeConstants (double upAttackMs, double upReleaseMs,
                           double downAttackMs, double downReleaseMs) noexcept;

    /** Processes in place. `channelData` must hold `numChannels` pointers to at
        least `numSamples` doubles. */
    void process (double* const* channelData, int numChannels, int numSamples);

    /** Zero: nothing in this engine looks ahead or transforms. */
    int getLatencySamples() const noexcept { return 0; }

    bool isPrepared() const noexcept { return prepared; }
    double getSampleRate() const noexcept { return currentSampleRate; }

    //==========================================================================
    // Design readout. The engine is the authority on what it actually runs, so
    // the tests measure against these rather than re-deriving the design.

    /** Band index 0..3 in the paper's numbering. Band 2 (index 1) has no filter
        of its own — it is formed as 1 - band1 - band3 — so only 0, 2 and 3 have
        coefficients. */
    const BiquadCoeffs& getBandCoeffs (int band) const noexcept;

    double getBandGain (int band) const noexcept { return gains[(size_t) band]; }
    const std::array<double, 4>& getBandGains() const noexcept { return gains; }

    /** Band-defining response at a frequency in Hz, band 2 included. */
    std::complex<double> bandResponseAt (int band, double hz) const;

    /** Summed low-level response 1 + sum g_i B_i: the encode transfer function
        below every threshold, where all four compressors sit at unity gain. */
    std::complex<double> lowLevelResponseAt (double hz) const;
    double lowLevelGainDb (double hz) const;

    double getThresholdLevel() const noexcept { return threshold; }
    double getClipLevel() const noexcept { return clipLevel; }

private:
    struct BandState
    {
        double peak = 0.0;      // peak follower on the rectified band signal
        double average = 0.0;   // averaged rectified band signal, scaled to read amplitude
        double gain = 1.0;      // linear limiter's smoothed gain
    };

    struct Channel;

    /** Solves the encoder's sample map for the sample that produced `y`. */
    double solveForInput (Channel& c, double y) const noexcept;

    void updateControl (BandState& s, double bandSample) const noexcept;

    void rebuildDetectorWeights() noexcept;

    /** Derives the per-direction effective ballistics from the user times. */
    void rebuildEffectiveTimes() noexcept;

    std::vector<std::unique_ptr<Channel>> channels;

    // Bands 1, 3 and 4; band 2 is derived.
    std::array<BiquadCoeffs, 3> bandCoeffs {};
    std::array<double, 4> gains { 0.0, 0.0, 0.0, 0.0 };

    // The user's side-chain filters, combined into one per-band weight on the
    // detector inputs.
    std::array<double, 4> detectorWeight { 1.0, 1.0, 1.0, 1.0 };
    double detectorHighPassHz = 0.0;
    double detectorLowPassHz = 0.0;

    double currentSampleRate = 48000.0;
    double peakLevel = 1.0;
    double threshold = 0.01;
    double clipLevel = 0.028;

    double peakAttackCoef = 1.0, peakReleaseCoef = 1.0, averageCoef = 1.0;
    double attackRate = 1.0, attackCoefMax = 1.0, recoveryCoef = 1.0;

    // The user's per-direction attack/release times; zero keeps the published
    // constant for that direction.
    double upAttackMs = 0.0, upReleaseMs = 0.0;
    double downAttackMs = 0.0, downReleaseMs = 0.0;

    // Effective ballistics per direction, [0] encode and [1] decode, derived
    // from the user times here and in prepare. Equal to the published
    // attackRate / recoveryCoef while the times are zero.
    std::array<double, 2> effectiveAttackRate { 1.0, 1.0 };
    std::array<double, 2> effectiveRecoveryCoef { 1.0, 1.0 };

    AtypeMode mode = AtypeMode::encode;
    bool prepared = false;
};

} // namespace lime
