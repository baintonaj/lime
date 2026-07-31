/*
  ==============================================================================

    SrEngine — top level of the process B reinterpretation.

    Signal order, following Fig. 5.10 of the reference hardware manual:

        encode:  x -> spectral skewing -> HF stages -> LF stages -> antisaturation
        decode:  y -> antisaturation^-1 -> LF stages^-1 -> HF stages^-1 -> skewing^-1

    The fixed networks are double-precision biquad cascades and invert exactly by
    swapping numerator and denominator; the two staggered halves are STFT cascades
    at their own resolutions and invert by per-bin division. Nothing runs in
    parallel, so there is no alignment delay: latency is simply the sum of the two
    analysis windows.

    Channels are processed independently. The 1967 paper notes the dual-path
    topology's tolerance of gain error "enables the noise reduction system to
    operate without control signal interconnections", and the reference hardware quotes
    crosstalk below -100 dB, so there is no stereo linking.

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"
#include "ControlSmoother.h"
#include "DelayLine64.h"
#include "FilterDesign.h"
#include "SpectrumProbe.h"
#include "SrSideChain.h"
#include "TapeChannel.h"

#include <atomic>
#include <memory>
#include <vector>

namespace lime
{

enum class SrMode
{
    encode,   // the signal you would print to tape
    decode,   // the complementary inverse
    loop      // encode -> tape channel -> decode, all in one instance
};

enum class SrProcess
{
    off,
    atype,               // handled by AtypeEngine, not here
    spectralRecording
};

//==============================================================================
/** Window and hop geometry for a sample rate.

    the paper's constants are all specified in Hz and milliseconds, so every size scales
    with the sample rate to keep bin width constant in Hz and latency constant in
    milliseconds. The geometry is one power-of-two exponent relative to 48 kHz,
    clamped to [-2, +1]:

      exponent  rates            HF window  HF hop  LF window  LF hop  latency
        -2      <= ~17 kHz          256        64      1024      256     1280
        -1      ~17 to ~34 kHz      512       128      2048      512     2560
         0      ~34 to ~68 kHz     1024       256      4096     1024     5120
        +1      >= ~68 kHz         2048       512      8192     2048    10240

    The upper clamp is a cost decision: 176.4/192 kHz reuse the 96 kHz sizes, which
    coarsens bins to 11.7 Hz rather than 5.9 Hz. The 40 Hz skewing corner still lands
    near bin 3, which is adequate, and latency in milliseconds improves. The lower
    clamp stops the transforms becoming degenerate at rates nothing musical uses.

    Each window is a whole multiple of its own hop, so every delay in the engine is an
    integer number of frames and the encoder's and decoder's analysis framing stay
    aligned — which is what lets their control state match exactly. The LF hop is a whole
    multiple of the HF hop as well, so the two halves' frames still land together.

    The halves ran at one shared hop of 256 until it was measured. That put the long
    window at sixteen-fold overlap against the short window's four, and re-analysed every
    low-frequency sample sixteen times for almost nothing: a 4096-point window already
    smears across 85 ms, so the extra frames cannot resolve anything the window has
    already blurred. Since the LF half carries 2049 of the 2562 bins, that redundancy was
    around three quarters of the engine's cost. The sqrt-Hann pair is COLA-exact for any
    hop dividing size/2 (see WolaWindow), so a quarter of the frames reconstructs just as
    exactly.
*/
struct SrGeometry
{
    int rateExponent = 0;

    int hfOrder = 10;   // 1024 at 48 kHz
    int lfOrder = 12;   // 4096 at 48 kHz
    int hfHop = 256;
    int lfHop = 1024;

    /** The largest bin counts any rate can produce, from the exponent clamp above:
        +1 puts the HF window at 2048 and the LF window at 8192. The engine's display
        probes are allocated to these once, so no re-prepare ever reallocates them. */
    static constexpr int maxHfBins = (1 << 11) / 2 + 1;
    static constexpr int maxLfBins = (1 << 13) / 2 + 1;

    int hfSize() const noexcept { return 1 << hfOrder; }
    int lfSize() const noexcept { return 1 << lfOrder; }

    /** One pass through both halves. */
    int passLatencySamples() const noexcept { return hfSize() + lfSize(); }

    static SrGeometry forSampleRate (double sampleRate);
};

//==============================================================================
class SrEngine
{
public:
    SrEngine();
    ~SrEngine();

    SrEngine (const SrEngine&) = delete;
    SrEngine& operator= (const SrEngine&) = delete;

    void prepare (double sampleRate, int maxBlockSize, int numChannels);
    void reset();

    const SrGeometry& getGeometry() const noexcept { return geometry; }
    bool isPrepared() const noexcept { return prepared; }

    /** Constant regardless of mode.

        Loop mode genuinely costs two passes through the halves; encode and decode
        cost one. Reporting whichever the current mode needs would change the plugin's
        latency when a parameter moves, which hosts handle badly, so the larger figure
        is reported always and the single-pass modes are padded to match.

        The EXPERIMENTAL filter-bank realisation applies its gain in the time domain
        with no signal delay at all, so under it both the report and the padding are
        zero — reporting the window against a path that has none would misalign the
        wet path by a fifth of a second. */
    int getLatencySamples() const noexcept
    {
        return realisation == GainRealisation::filterBank
                 ? 0 : 2 * geometry.passLatencySamples();
    }

    void setMode (SrMode m) noexcept { mode = m; }
    void setProcess (SrProcess p) noexcept { processType = p; }

    SrMode getMode() const noexcept { return mode; }
    SrProcess getProcess() const noexcept { return processType; }

    void setSlidingComposition (SlidingComposition c) noexcept { composition = c; }

    /** Spectral smoothing width in bins, applied to both halves. See
        SrSideChain::setSpectralSmoothingBins. */
    void setSpectralSmoothing (double bins);

    /** Critical-band width of the level-envelope smoothing. See
        SrSideChain::setEnvelopeSmoothingBark. */
    void setEnvelopeSmoothing (double bark);

    /** Sets how finely the band-defining split between the two staggered halves is
        resolved, in both directions at once so an encode and a decode in the same instance
        cannot disagree. See SrSideChain::setSections for the range and why it exists. */
    void setSections (int count);

    /** Read by the editor while the audio thread writes; atomic so neither side
        can see a torn value. */
    int getSections() const noexcept { return sections.load (std::memory_order_relaxed); }

    /** The user's side-chain high pass: high-passes what every compressor's level
        detection reads, without touching the signal path. Pushed to all four side
        chains at once for the same reason setSections is — quoting the two
        directions different corners is exactly what would break the
        complementarity. Zero or negative disables it. See
        SrSideChain::setDetectorHighPass. */
    void setSidechainHighPass (double cornerHz);

    /** The matching low pass on the detection; the two corners band-limit what
        the compressors hear. Same fan-out and same discipline. See
        SrSideChain::setDetectorLowPass. */
    void setSidechainLowPass (double cornerHz);

    /** The user's attack/release times, one pair per direction.

        The up pair goes to both encode side chains and the down pair to both
        decode side chains — the two halves of one detection must always agree.
        The pairs themselves may differ, and that is the deliberate trade: the
        round trip is exact only while up and down run the same ballistics, so
        setting them apart is a creative choice that forfeits the exact null
        until they agree again. Zero for any value keeps that direction's
        published tau. See SrSideChain::setTimeConstants. */
    void setTimeConstants (double upAttackMs, double upReleaseMs,
                           double downAttackMs, double downReleaseMs);

    /** The three calibration trims that sit *inside* the loop, as linear gains.

        In loop mode the channel is inside the instance, so send out, return in and return
        out have nowhere else to be: send out scales what is printed to the channel,
        return in scales what comes back before the decoder reads it, and return out is
        the make-up after decoding. Applied only in loop mode; encode and decode take
        their pair from the wrapper, where the process has an outside.

        Smoothed, so any of the three can be swept while audio is running. Unity by
        default, and at unity the loop is exactly what it was before they existed.

        @param immediately  true from prepareToPlay, where a restored session's settings
                            have to be in force on the first sample rather than faded in. */
    void setLoopTrims (double sendOutGain, double returnInGain, double returnOutGain,
                       bool immediately = false) noexcept;

    /** How the gain profile is applied. See GainRealisation. */
    void setGainRealisation (GainRealisation);

    TapeChannel& getTapeChannel() noexcept { return tape; }

    /** Switches the encode probes on while an editor is watching.

        Publishing a frame costs two logarithms per bin in both encode side chains, which
        is not something a session with no window open should be paying for. Off until an
        editor says otherwise, and re-applied across prepare(). */
    void setProbesActive (bool shouldBeActive) noexcept;

    /** Live encode-curve snapshots for the editor.

        The probes live on the engine itself rather than inside the per-channel side
        chains, so the pointers stay valid for the engine's whole lifetime — a
        re-prepare rebuilds the channels but only re-wires these, and an editor that
        fetched them at any point can keep reading them safely. Channel 0 only, by
        design: the display shows one channel's encode curve, and the reference
        hardware processes channels without any control interconnection anyway.
        Null for any other channel index. */
    SpectrumProbe* getHfProbe (int channel) noexcept;
    SpectrumProbe* getLfProbe (int channel) noexcept;

    /** Dominants seen in the last block, summed across both halves and every
        channel. Metering; atomic so a display thread can sample it safely. */
    int getNumDominants() const noexcept { return lastDominants.load (std::memory_order_relaxed); }

    /** Processes in place. */
    void process (double* const* channelData, int numChannels, int numSamples);

private:
    struct Channel;

    void encodePass (Channel&, double* buffer, int numSamples);
    void decodePass (Channel&, double* buffer, int numSamples);

    std::vector<std::unique_ptr<Channel>> channels;
    SrGeometry geometry;
    TapeChannel tape;

    SrMode mode = SrMode::loop;
    SrProcess processType = SrProcess::spectralRecording;
    SlidingComposition composition = SlidingComposition::restriction;

    // The loop's internal trims. See setLoopTrims.
    ControlSmoother sendOutTrim, returnInTrim, returnOutTrim;

    // Engine-lifetime display probes, allocated once to the largest geometry any
    // rate can need; prepare() re-wires channel 0's encode side chains to them.
    // See getHfProbe for why they are not per-channel members.
    SpectrumProbe hfDisplayProbe { SrGeometry::maxHfBins };
    SpectrumProbe lfDisplayProbe { SrGeometry::maxLfBins };

    // Both written from the message thread (editor open/close) and read across
    // prepare(); atomic so the two never race.
    std::atomic<bool> probesActive { false };
    std::atomic<int> sections { SrSideChain::defaultSections };

    double currentSampleRate = 48000.0;
    int maxBlock = 0;
    std::atomic<int> lastDominants { 0 };
    double smoothing = 0.0;
    double envelopeBark = 1.0;
    double sidechainHighPassHz = 0.0;
    double sidechainLowPassHz = 0.0;

    // The user's per-direction attack/release times, cached across prepare.
    double upAttackMs = 0.0, upReleaseMs = 0.0;
    double downAttackMs = 0.0, downReleaseMs = 0.0;
    GainRealisation realisation = GainRealisation::stft;
    bool prepared = false;

    std::vector<double*> tapePointers;
};

} // namespace lime
