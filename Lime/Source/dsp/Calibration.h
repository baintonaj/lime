/*
  ==============================================================================

    Calibration — the line-up signals and Auto Compare.

    Section 5 of the reference hardware manual, and section 4 of the 1987 paper,
    describe three things an installation needs before any noise reduction is
    meaningful: a broadband signal for checking recorder response, a tone for
    checking level, and a way of hearing tape against source without having to
    trust a meter.

      - calibration noise: pink noise interrupted by 20 ms "nicks" every 2 s, recorded
        15 dB below reference level. The attenuation is deliberate. At reference level
        the noise would drive the tape into saturation at both frequency
        extremes, so what came back would show the saturation characteristic
        rather than the recorder's response.
      - calibration tone: 850 Hz, modulated *in frequency* up 10 percent for 30 ms
        every 750 ms, at reference level. The manual is specific that the modulation
        is FM and not AM so that meters read a steady level throughout, which is
        why the oscillator here keeps phase continuous and amplitude untouched.
      - Auto Compare: the decoder alternates tape return and an internally
        generated *uninterrupted* pink noise in 4 s segments. The two are told
        apart by ear: the tape carries the nicks, the internal reference does
        not. The hardware lights red for the internal reference and green for
        tape, so the active source is exposed here for an indicator to follow.

    Level convention: reference level is amplitude 1.0, meaning the
    reference tone is a sine of *peak* amplitude 1.0 and therefore RMS
    1/sqrt(2). Every "dB relative to reference level" figure in this file is an RMS
    ratio against that, so calibration noise at -15 dB has RMS
    10^(-15/20)/sqrt(2) = 0.12574, i.e. -18.0 dBFS. Nothing here can clip: the
    noise source is bounded and the tone is exactly unit amplitude.

    No JUCE, no float, and no allocation after prepare().

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"

#include <atomic>
#include <cstdint>
#include <vector>

namespace lime
{

//==============================================================================
/** The published constants. Every figure is from section 5 of the reference hardware manual. */
namespace calibration
{
    /** reference level as a normalised peak amplitude. */
    inline constexpr double referenceAmplitude = 1.0;

    /** RMS of a sine at reference level, which is what the dB figures below are
        relative to. */
    inline constexpr double referenceRms = 0.70710678118654752440;

    inline constexpr double noiseBelowReferenceDb = -15.0;

    /** 10^(-15/20), spelled out because std::pow is not constexpr. */
    inline constexpr double noiseGainFromReference = 0.177827941003892280;

    /** RMS of calibration noise between the nicks. */
    inline constexpr double noiseRms = referenceRms * noiseGainFromReference;

    inline constexpr double nickSeconds = 0.020;
    inline constexpr double nickPeriodSeconds = 2.0;

    inline constexpr double toneHz = 850.0;
    inline constexpr double toneModulationFactor = 1.10;
    inline constexpr double toneModulatedHz = 935.0;      // 850 * 1.10
    inline constexpr double toneBurstSeconds = 0.030;
    inline constexpr double tonePeriodSeconds = 0.750;

    inline constexpr double autoCompareSegmentSeconds = 4.0;

    inline constexpr std::uint64_t defaultSeed = 0x5352'4e4f'4953'4531ull;   // "SRNOISE1"
}

//==============================================================================
/** Deterministic uniform white noise on [-1, 1).

    xoshiro256++ over a splitmix64-expanded seed. Written out rather than taken
    from <random> so the sequence is fixed by this file alone: the standard
    pins down mt19937_64's output but not, for instance, normal_distribution's,
    and the tests assert bit-exact reproducibility.

    Uniform rather than Gaussian, which differs from the analog hardware's zener
    source, for two reasons: the peak is bounded, so a calibration signal can be
    guaranteed not to clip, and no transcendental is needed, so the sequence is
    reproducible on any libm. The pinking filter that follows has several
    thousand samples of memory, so its output is very close to Gaussian anyway.
*/
class WhiteNoise64
{
public:
    explicit WhiteNoise64 (std::uint64_t seed = calibration::defaultSeed) noexcept;

    /** Reseeds and rewinds. */
    void setSeed (std::uint64_t seed) noexcept;

    /** Rewinds to the start of the current seed's sequence. */
    void reset() noexcept;

    double nextSample() noexcept;

    /** RMS of the sequence, 1/sqrt(3). */
    static constexpr double rms() noexcept { return 0.57735026918962576451; }

private:
    std::uint64_t seedValue = calibration::defaultSeed;
    std::uint64_t state[4] {};
};

//==============================================================================
/** The pinking filter: a cascade of first-order shelving sections giving a
    -3 dB/octave magnitude.

    Method. A pole/zero pair spaced a fixed ratio apart in log frequency
    produces -6 dB/octave between them and 0 dB above, so a staggered chain of
    them averages to any slope between the two. Poles are laid out four per
    decade from 4 Hz to 0.47*sampleRate; the *depth* of each section's shelf is
    then solved for by damped Gauss-Newton (Levenberg-Marquardt) against an
    exact -3 dB/octave target over 20 Hz to min(20 kHz, 0.45*sampleRate), with
    an overall gain as a free parameter.

    Solving numerically rather than placing the zeros by the asymptotic rule
    matters more than it looks. The bilinear transform maps digital frequency f
    onto analog 2*fs*tan(pi*f/fs), so a cascade that is exactly pink in the
    analog variable is pink in the *warped* frequency, which is 0.7 dB out by
    10 kHz and 2.2 dB out by 16 kHz at 48 kHz. Asymptotic corrections for that
    over- or under-shoot by a few tenths of a dB; the fit removes it entirely.

    Measured: 0.017 dB peak-to-peak or better over 20 Hz to min(20 kHz,
    0.45*sampleRate), at every rate from 8 kHz to 384 kHz. The fit runs once per
    prepare() and takes a few milliseconds. Below the lowest pole the spectrum
    necessarily flattens, which is what keeps the total power finite.
*/
class PinkFilter64
{
public:
    void prepare (double sampleRate);
    void reset() noexcept;

    double processSample (double x) noexcept { return cascade.processSample (x); }

    /** RMS gain for white input: sqrt of the mean of |H|^2 over the band.
        Computed by quadrature at prepare() time, so the output level of a
        generator using this filter is exact rather than measured. */
    double getNoiseGain() const noexcept { return noiseGain; }

    /** Worst deviation from -3 dB/octave over the fitted band, peak to peak. */
    double getWorstDeviationDb() const noexcept { return worstDeviationDb; }

    /** Magnitude in dB at a frequency in Hz. Design-time helper, not real time. */
    double responseDbAt (double hz) const;

    int getNumSections() const noexcept { return (int) sections.size(); }

    double getFitLowHz()  const noexcept { return fitLowHz; }
    double getFitHighHz() const noexcept { return fitHighHz; }

private:
    std::vector<BiquadCoeffs> sections;
    BiquadCascade64 cascade;

    double currentSampleRate = 48000.0;
    double noiseGain = 1.0;
    double worstDeviationDb = 0.0;
    double fitLowHz = 20.0, fitHighHz = 20000.0;
};

//==============================================================================
/** Uninterrupted pink noise: the internal reference Auto Compare switches
    against, and the source calibration noise is gated from.

    Defaults to the same RMS as calibration noise, 15 dB below reference level,
    because Auto Compare is only meaningful if the two sides match in level.
*/
class PinkNoiseGenerator
{
public:
    void prepare (double sampleRate);
    void reset();

    /** Reseeds and rewinds the white source immediately. The filter state is
        left alone; call reset() as well for a signal identical to a fresh
        instance. */
    void setSeed (std::uint64_t seed) noexcept;

    void setRmsLevel (double rms) noexcept;
    double getRmsLevel() const noexcept { return targetRms; }

    double nextSample() noexcept { return scale * filter.processSample (white.nextSample()); }

    /** Overwrites `out`. */
    void process (double* out, int numSamples) noexcept;

    const PinkFilter64& getFilter() const noexcept { return filter; }

private:
    void updateScale() noexcept;

    WhiteNoise64 white;
    PinkFilter64 filter;

    double targetRms = calibration::noiseRms;
    double scale = 1.0;
};

//==============================================================================
/** calibration noise: the pink source gated off for 20 ms every 2 s.

    The noise source and its filter keep running through the nick and only the
    output is muted, exactly as the hardware gates a continuously running
    generator. That way the noise either side of a nick is one uninterrupted
    process, with no filter transient to be mistaken for a recorder fault.

    The nick occupies the last 20 ms of each 2 s period, so the signal opens
    with noise and the first nick begins at 1.980 s.

    The gate is hard-edged. The click is the point: it is what an operator, and
    the Auto Compare indicator, uses to tell tape from the internal reference.
*/
class CalibrationNoiseGenerator
{
public:
    void prepare (double sampleRate);
    void reset();

    void setSeed (std::uint64_t seed) noexcept;

    /** Overwrites `out`. */
    void process (double* out, int numSamples) noexcept;

    /** True if the *next* sample generated falls inside a nick. */
    bool isInNick() const noexcept { return phase >= nickPeriodSamples - nickLengthSamples; }

    int getNickLengthSamples() const noexcept { return nickLengthSamples; }
    int getNickPeriodSamples() const noexcept { return nickPeriodSamples; }

    const PinkNoiseGenerator& getNoise() const noexcept { return noise; }

private:
    PinkNoiseGenerator noise;

    int nickLengthSamples = 0;
    int nickPeriodSamples = 1;
    int phase = 0;
};

//==============================================================================
/** calibration tone: 850 Hz, frequency modulated up 10 percent to 935 Hz for 30 ms
    every 750 ms, at reference level.

    One phase accumulator whose increment switches between the two frequencies.
    Phase is therefore continuous across every transition and the amplitude term
    is never touched, so there is no envelope ripple at all — which is the whole
    reason the manual gives for the modulation being FM rather than AM: a meter
    must read one steady level through the burst.

    The burst occupies the last 30 ms of each 750 ms period, matching the nick
    placement in CalibrationNoiseGenerator, so the first burst begins at 0.720 s.
*/
class CalibrationToneGenerator
{
public:
    void prepare (double sampleRate);
    void reset();

    /** Peak amplitude; 1.0, i.e. reference level, by default. */
    void setAmplitude (double newAmplitude) noexcept { amplitude = newAmplitude; }
    double getAmplitude() const noexcept { return amplitude; }

    /** Overwrites `out`. */
    void process (double* out, int numSamples) noexcept;

    /** True if the *next* sample generated falls inside a modulation burst. */
    bool isModulated() const noexcept { return phase >= periodSamples - burstSamples; }

    /** Frequency the next sample will be generated at, in Hz. */
    double getCurrentFrequencyHz() const noexcept;

    int getBurstSamples()  const noexcept { return burstSamples; }
    int getPeriodSamples() const noexcept { return periodSamples; }

private:
    double currentSampleRate = 48000.0;
    double amplitude = calibration::referenceAmplitude;

    double angle = 0.0;              // radians
    double baseIncrement = 0.0;
    double modulatedIncrement = 0.0;

    int burstSamples = 0;
    int periodSamples = 1;
    int phase = 0;
};

//==============================================================================
/** Decides whether the incoming signal is calibration noise.

    The manual only asks for enough to stop Auto Compare running on material it
    cannot make sense of, and warns that the circuit "may deliver spurious and
    unexpected noises if it receives signals other than calibration noise (including
    silence); this is normal and harmless". So this is deliberately a plain
    pattern matcher rather than a classifier, and it errs towards not firing.

    Three conditions, all of which must hold:

      - Level. The slow RMS must exceed -60 dBFS. Silence, and anything close to
        it, is rejected here.
      - Noise-like. Any sinusoid satisfies x[n+1] + x[n-1] = 2*cos(w)*x[n]
        exactly, so the residual of that one-coefficient predictor, relative to
        the signal it is predicting, separates noise from tones by a very wide
        margin. Measured: -7.0 dB for pink noise, -28.2 dB for pink noise through
        three cascaded 5 kHz one-poles, -43.7 dB through three cascaded 2 kHz
        one-poles, and -141 dB or the numerical floor for a sine at any
        frequency. The threshold is -55 dB, so a tape return can be in very poor
        shape and still pass while a tone has no chance. Note that this rejects
        tones rather than proving broadband content: narrowband noise passes.
        That is deliberate — see below.
      - Nicks. A gap is a stretch where the 0.7 ms mean square falls more than
        20 dB below the 400 ms mean square. Gaps lasting 8 to 45 ms count as
        nicks, and two consecutive nick-to-nick intervals must both fall between
        1.7 and 2.3 s.

    The nick pattern does nearly all the work; the other two conditions only stop
    the sequence running on material where it would be meaningless. Consequences
    worth knowing: three nicks are needed, so detection arms about 6 s after the
    noise starts, and one missed nick, or 2.6 s without one, disarms it. A steady
    tone and uninterrupted pink noise both fail on the absence of gaps, not on
    the noise-like test.

    The gate's decay means a detected gap runs about 3 ms short of the true 20 ms
    and starts about 3 ms late. The offset is the same for every nick, so the
    measured intervals are unaffected.
*/
class CalibrationNoiseDetector
{
public:
    void prepare (double sampleRate);
    void reset();

    void push (const double* samples, int numSamples) noexcept;

    bool isCalibrationNoisePresent() const noexcept { return present; }

    /** Diagnostics, for tests and for a panel that wants to show why the
        sequence is not running. */
    int getNickCount() const noexcept { return nickCount; }
    int getValidIntervals() const noexcept { return validIntervals; }
    double getLevelDbFs() const noexcept;

    /** Sinusoid-predictor residual relative to the signal, in dB. Near 0 dB for
        white noise, about -7 dB for pink noise, below -120 dB for a pure tone. */
    double getNoiseLikenessDb() const noexcept;

private:
    double fastCoeff = 0.0, slowCoeff = 0.0, predictorCoeff = 0.0;
    double fastPower = 0.0, slowPower = 0.0;

    // Accumulators for the one-coefficient sinusoid predictor.
    double midPower = 0.0, sumPower = 0.0, crossPower = 0.0;
    double history1 = 0.0, history2 = 0.0;

    long long position = 0;
    long long gapStart = 0;
    long long lastNickStart = -1;

    int minGapSamples = 0, maxGapSamples = 0;
    int minIntervalSamples = 0, maxIntervalSamples = 0;
    int timeoutSamples = 0;

    int nickCount = 0;
    int validIntervals = 0;
    bool inGap = false;
    bool present = false;
};

//==============================================================================
enum class AutoCompareSource
{
    tape,        // green on the hardware; carries the nicks
    reference    // red on the hardware; internal, uninterrupted
};

//==============================================================================
/** Auto Compare: 4 s of tape, 4 s of internal reference, repeating.

    Hard switched, not crossfaded. The hardware has nothing to crossfade with,
    and a fade would blur the level step that the operator is listening for.
    Both inputs are read every sample so that neither side is ever asked to
    start from cold at a boundary.

    Starts on tape, so the first thing heard after arming is what is actually on
    the tape. `setRunning (false)`, and `reset()`, park the output on tape and
    rewind the cycle; drive the former from CalibrationNoiseDetector, since section 5
    has the sequence run only while calibration noise is present.
*/
class AutoCompare
{
public:
    void prepare (double sampleRate);
    void reset();

    void setRunning (bool shouldRun) noexcept;
    bool isRunning() const noexcept { return running.load (std::memory_order_relaxed); }

    /** Writes the currently selected source to every channel of `out`; `out` may
        alias `tape`. One call per block, whatever the channel count: the segment
        clock advances once per *sample*, so a stereo pair switches source on the
        same sample and the segments stay the configured length however many
        channels there are. */
    void process (const double* const* tape, const double* reference,
                  double* const* out, int numChannels, int numSamples) noexcept;

    /** Single-channel convenience over the multi-channel form. */
    void process (const double* tape, const double* reference, double* out, int numSamples) noexcept;

    AutoCompareSource getCurrentSource() const noexcept;

    /** True while the internal reference is selected: the red indicator. */
    bool isShowingReference() const noexcept { return getCurrentSource() == AutoCompareSource::reference; }

    double getSecondsIntoSegment() const noexcept;

    int getSegmentSamples() const noexcept { return segmentSamples; }

private:
    double currentSampleRate = 48000.0;
    int segmentSamples = 1;

    // Written on the audio thread, read by the editor's indicator on the message
    // thread; atomic so the front panel never sees a torn value. The audio thread
    // works on a local copy and stores once per block.
    std::atomic<int> phase { 0 };
    std::atomic<bool> running { false };
};

} // namespace lime
