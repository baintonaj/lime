/*
  ==============================================================================

    TapeChannel — the analog tape impairments process B was designed to defeat.

    Loop mode runs encode -> this -> decode inside one instance, so this is the
    thing in the middle. It is not a tape emulation for its own sake: it models
    only the four defects the papers actually name, so that the noise reduction
    can be heard working and measured.

      1. Hiss. "Statistical fluctuations in magnetic domain orientation in the
         tape coating" — broadband, but spectrally shaped by the replay
         equalisation rather than white, which is why the papers insist noise be
         measured weighted.
      2. Gradual saturation. The papers contrast analog tape's "gradual or soft
         clipping characteristic" with digital hard clipping, and stress that
         headroom is worst at the frequency extremes: the 3180 us low-frequency
         and 35 us high-frequency record equalisations boost those regions onto
         the tape, so they reach the coating's limit first. That asymmetry is the
         entire reason section 2.7 needs antisaturation networks. Modelled the
         way the machine causes it — record equalisation, one smooth saturator,
         exactly inverse replay equalisation — rather than as three separately
         distorted bands. The emphasis pair is an exact biquad inverse, so at low
         level the stage is transparent to numerical precision.
      3. High-frequency loss. Head and gap loss rolling off the top.
      4. Modulation noise. Noise that "appears only when a signal is present",
         spread over a wide range of spectrum around the signal and scaling with
         its level. Modelled multiplicatively, which is both the physical
         mechanism (coating irregularity modulating the recorded signal) and the
         only model that puts the noise around the signal rather than at
         baseband. A no-signal noise measurement misses it entirely, which is the
         papers' point.

    Double precision throughout, no JUCE, standard library only. Channels are
    processed independently with independent noise streams, matching the rest of
    Lime: there is no interchannel linking anywhere in this project.

    Each defect's level is normalised through its own shaping filter, so what is
    asked for in dB is what is injected. Hiss arises at replay, after every
    other defect, so its measured level is fully independent of the other
    settings. Modulation noise is not: it originates at the tape/head interface,
    so it is injected before the head loss and rides through that filter along
    with the signal it modulates — its delivered level and spectrum shift with
    the head-loss corner, which is the physically correct ordering rather than a
    coupling to normalise away.

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace lime
{

//==============================================================================
/** Fixed constants of the tape model. Everything here is either a published
    equalisation time constant or a documented modelling choice; nothing is
    exposed as a user parameter, in keeping with decision 14 of the plan.
*/
namespace tapeModel
{
    /** Digital amplitude taken to be reference level. Every level in
        TapeChannelSettings is relative to this, so 20 dB of peak headroom sits
        above reference — the same order as the analog chain's. */
    inline constexpr double referenceAmplitude = 0.1;

    /** Midband headroom above reference before the saturator's knee. Tape is
        usually quoted at 10 to 14 dB above reference level for a few percent
        distortion; the extremes lose far more than this, via the record
        equalisation below. */
    inline constexpr double midbandHeadroomDb = 12.0;

    /** 3180 us NAB low-frequency record equalisation, 1/(2*pi*T). Speed
        independent: both 15 and 30 ips use 3180 us. */
    inline constexpr double nabLowCornerHz = 50.0472975;

    /** 35 us CCIR high-frequency record equalisation at 15 ips. Halves in time,
        so doubles in frequency, at 30 ips — faster tape has better high
        frequency headroom, and this is where that comes from. */
    inline constexpr double ccirHighCornerHz = 4547.2841;

    /** Depth of the record equalisation shelves, and so the headroom lost at the
        extremes. Chosen to mirror section 4.1's antisaturation figures (10 dB at
        25 Hz and 15 kHz), since those networks exist precisely to pre-compensate
        this loss. */
    inline constexpr double recordLowBoostDb = 12.0;
    inline constexpr double recordHighBoostDb = 12.0;

    /** Hiss shaping: a low shelf for the replay low-frequency lift, a high shelf
        for the replay high-frequency lift, and a two-pole limit at the end of the
        chain's bandwidth. The shelf and bandwidth corners scale with tape speed;
        overall level does not, because it is what the user asks for. */
    inline constexpr double hissLowShelfHz = 150.0;
    inline constexpr double hissLowShelfDb = 9.0;
    inline constexpr double hissHighShelfDb = 5.0;
    inline constexpr double hissBandwidthHz = 15000.0;

    /** Modulation noise: the modulating process is band limited, so its sidebands
        form a skirt of roughly this half-width either side of every signal
        component. High-passed so the result reads as noise rather than as a slow
        level drift, which would be flutter and is a different defect. */
    inline constexpr double modulationSpreadHz = 1000.0;
    inline constexpr double modulationHighPassHz = 20.0;

    /** Tape speed reference and the range the model will accept. */
    inline constexpr double referenceSpeedIps = 15.0;
    inline constexpr double minSpeedIps = 3.75;
    inline constexpr double maxSpeedIps = 30.0;
}

//==============================================================================
/** The soft saturation curve, normalised so the knee sits at unit input.

    Odd, C-infinity, strictly increasing with unit slope at the origin, and
    asymptotic to +/-1: it compresses but never folds, and it leaves small
    signals alone to within a relative error of x^2/3. Exposed so tests can
    verify monotonicity of the nonlinearity directly rather than inferring it
    from the channel's output.
*/
inline double tapeSaturationCurve (double x) noexcept
{
    return std::tanh (x);
}

//==============================================================================
struct TapeChannelSettings
{
    double hissLevelDb = -65.0;          // RMS relative to reference level
    double saturationDb = 0.0;           // extra drive into the saturator
    double hfLossHz = 18000.0;           // -3 dB corner of the head/gap loss, at 15 ips
    double modulationNoiseDb = -55.0;    // RMS relative to the signal driving it
    double tapeSpeedIps = 15.0;          // 15 or 30 ips; shifts hiss shape and HF loss

    bool hissEnabled = true;
    bool saturationEnabled = true;
    bool hfLossEnabled = true;
    bool modulationNoiseEnabled = true;
};

//==============================================================================
class TapeChannel
{
public:
    TapeChannel();
    ~TapeChannel();

    TapeChannel (const TapeChannel&) = delete;
    TapeChannel& operator= (const TapeChannel&) = delete;

    /** Allocates all state. `maxBlockSize` is accepted for symmetry with the
        other engines and is not used: nothing here is block-sized, so process()
        is indifferent to how the host divides the stream. */
    void prepare (double sampleRate, int maxBlockSize, int numChannels);

    /** Clears filter state and restarts every noise stream from the current
        seed, so a reset channel reproduces its previous output exactly. */
    void reset();

    /** Redesigns the coefficients. Reallocates nothing once prepared, but does
        run the noise normalisation integrals, so call it on parameter change
        rather than per block. */
    void setSettings (const TapeChannelSettings& newSettings);

    const TapeChannelSettings& getSettings() const noexcept { return settings; }

    bool isPrepared() const noexcept { return prepared; }

    /** The head-loss corner actually in use: the requested corner scaled by tape
        speed and held below Nyquist. Before prepare() the Nyquist limit is
        unknown, so the speed-scaled value is reported unclamped. */
    double effectiveHfLossHz() const noexcept;

    /** Deterministic noise for tests. Reseeds immediately; each channel derives
        an independent stream from this value, and hiss and modulation noise draw
        from separate streams so enabling one cannot disturb the other. */
    void setSeed (uint64_t seed);

    /** Processes in place. `channelData` must hold `numChannels` pointers to at
        least `numSamples` doubles. Allocation free. */
    void process (double* const* channelData, int numChannels, int numSamples);

private:
    struct Channel;

    void updateCoefficients();
    void applyCoefficients();
    void seedChannels();

    std::vector<std::unique_ptr<Channel>> channels;
    TapeChannelSettings settings;

    // Designed once per settings or rate change, shared by every channel.
    std::vector<BiquadCoeffs> hissShapeCoeffs;
    std::vector<BiquadCoeffs> emphasisCoeffs;
    std::vector<BiquadCoeffs> headLossCoeffs;
    std::vector<BiquadCoeffs> modulationShapeCoeffs;

    double hissGain = 0.0;          // scales unit-variance shaped noise to the asked level
    double modulationGain = 0.0;    // ditto, relative to the signal it multiplies
    double saturationInScale = 1.0; // into the curve: drive / knee
    double saturationOutScale = 1.0;// out of it: knee / drive

    double currentSampleRate = 48000.0;
    uint64_t baseSeed = 0x9E3779B97F4A7C15ull;
    bool prepared = false;
};

} // namespace lime
