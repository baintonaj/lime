/*
  ==============================================================================

    SrEngine — geometry, per-channel state, and the encode and decode passes.

  ==============================================================================
*/

#include "SrEngine.h"

#include <algorithm>
#include <cmath>

namespace lime
{

//==============================================================================
SrGeometry SrGeometry::forSampleRate (double sampleRate)
{
    SrGeometry g;

    int exponent = 0;

    if (sampleRate > 0.0)
    {
        const double ratio = sampleRate / 48000.0;
        exponent = (int) std::lround (std::log2 (ratio));
        exponent = std::clamp (exponent, -2, 1);
    }

    g.rateExponent = exponent;
    g.hfOrder = 10 + exponent;
    g.lfOrder = 12 + exponent;
    g.hfHop = 1 << (8 + exponent);

    // Four-fold overlap on the long window, matching the short one, rather than sixteen.
    // See the note on SrGeometry.
    g.lfHop = 1 << (10 + exponent);

    return g;
}

namespace
{
    /** Applies a smoothed gain across a block.

        Samples outer, channels inner, deliberately. Advancing the smoother inside a
        per-channel loop would give the left and right of a stereo pair different
        trajectories through the same ramp, and they would drift apart for as long as it
        was moving. */
    void applyLoopTrim (ControlSmoother& smoother, double* const* channelData,
                        int numChannels, int numSamples)
    {
        if (! smoother.isMoving())
        {
            const auto gain = smoother.getCurrent();

            if (gain == 1.0)
                return;

            for (int ch = 0; ch < numChannels; ++ch)
            {
                auto* samples = channelData[ch];

                for (int n = 0; n < numSamples; ++n)
                    samples[n] *= gain;
            }

            return;
        }

        for (int n = 0; n < numSamples; ++n)
        {
            const auto gain = smoother.next();

            for (int ch = 0; ch < numChannels; ++ch)
                channelData[ch][n] *= gain;
        }
    }
}

//==============================================================================
struct SrEngine::Channel
{
    // The fixed networks. Skewing feeds the stages; antisaturation follows them.
    BiquadCascade64 skew, skewInverse;
    BiquadCascade64 antiSat, antiSatInverse;

    // Separate instances per direction, so loop mode can encode and decode in the
    // same block without either sharing control state with the other.
    SrSideChain hfEncode, lfEncode;
    SrSideChain hfDecode, lfDecode;

    /** Pads encode and decode up to the reported loop-mode latency, and carries the
        whole delay when the process is switched off. */
    DelayLine64 modeDelay;
};

//==============================================================================
SrEngine::SrEngine() = default;
SrEngine::~SrEngine() = default;

void SrEngine::prepare (double sampleRate, int maxBlockSize, int numChannels)
{
    currentSampleRate = sampleRate;
    maxBlock = std::max (1, maxBlockSize);
    geometry = SrGeometry::forSampleRate (sampleRate);

    // Spectral skewing alone drives the stages; the antisaturation networks are the
    // remainder of the fixed main path. Both are solved for this rate so the total
    // reproduces the section 4.1 figures.
    const auto skewProtos = srNetworks::skewing();
    const auto skewCoeffs = srNetworks::toBiquads (skewProtos, sampleRate);

    const auto mainCoeffs = srNetworks::mainPathBiquadsForRate (sampleRate);

    // Antisaturation is whatever the fitted main path adds beyond skewing: the fit
    // returns skewing first, then the two antisaturation sections.
    std::vector<BiquadCoeffs> antiSatCoeffs;

    for (size_t i = skewCoeffs.size(); i < mainCoeffs.size(); ++i)
        antiSatCoeffs.push_back (mainCoeffs[i]);

    channels.clear();
    channels.reserve ((size_t) std::max (0, numChannels));

    for (int ch = 0; ch < numChannels; ++ch)
    {
        auto c = std::make_unique<Channel>();

        c->skew.setCoeffs (skewCoeffs);
        c->skewInverse.setInverseCoeffs (skewCoeffs);
        c->antiSat.setCoeffs (antiSatCoeffs);
        c->antiSatInverse.setInverseCoeffs (antiSatCoeffs);

        // The slew limits are quoted per frame, so the half with the longer frame is
        // allowed proportionally more per frame to keep the same rate in dB per second.
        const double lfSlewScale = (double) geometry.lfHop / (double) geometry.hfHop;

        c->hfEncode.prepare (true,  geometry.hfOrder, geometry.hfHop, sampleRate, composition, 1.0, maxBlock);
        c->lfEncode.prepare (false, geometry.lfOrder, geometry.lfHop, sampleRate, composition, lfSlewScale, maxBlock);
        c->hfDecode.prepare (true,  geometry.hfOrder, geometry.hfHop, sampleRate, composition, 1.0, maxBlock);
        c->lfDecode.prepare (false, geometry.lfOrder, geometry.lfHop, sampleRate, composition, lfSlewScale, maxBlock);

        // The display probes are the engine's, with the engine's lifetime; only the
        // first channel's encode side chains publish into them. See getHfProbe.
        c->hfEncode.setProbe (ch == 0 ? &hfDisplayProbe : nullptr);
        c->lfEncode.setProbe (ch == 0 ? &lfDisplayProbe : nullptr);

        c->hfEncode.setDirection (SrDirection::encode);
        c->lfEncode.setDirection (SrDirection::encode);
        c->hfDecode.setDirection (SrDirection::decode);
        c->lfDecode.setDirection (SrDirection::decode);

        c->hfEncode.setGainRealisation (realisation);
        c->lfEncode.setGainRealisation (realisation);
        c->hfDecode.setGainRealisation (realisation);
        c->lfDecode.setGainRealisation (realisation);

        c->hfEncode.setEnvelopeSmoothingBark (envelopeBark);
        c->lfEncode.setEnvelopeSmoothingBark (envelopeBark);
        c->hfDecode.setEnvelopeSmoothingBark (envelopeBark);
        c->lfDecode.setEnvelopeSmoothingBark (envelopeBark);

        c->hfEncode.setSpectralSmoothingBins (smoothing);
        c->lfEncode.setSpectralSmoothingBins (smoothing);
        c->hfDecode.setSpectralSmoothingBins (smoothing);
        c->lfDecode.setSpectralSmoothingBins (smoothing);

        c->modeDelay.prepare (2 * geometry.passLatencySamples());

        channels.push_back (std::move (c));
    }

    if (! channels.empty())
    {
        hfDisplayProbe.prepare (channels.front()->hfEncode.getNumBins(), sampleRate, geometry.hfSize());
        lfDisplayProbe.prepare (channels.front()->lfEncode.getNumBins(), sampleRate, geometry.lfSize());
    }

    // The halves come out of prepare at the default resolution, so a section count set
    // before prepareToPlay is not lost. Same for the probes: a sample-rate change rebuilds
    // every side chain, and an editor that is already open is not going to ask twice.
    // The side-chain filters ride along for the same reason.
    setSections (sections.load (std::memory_order_relaxed));
    setProbesActive (probesActive.load (std::memory_order_relaxed));

    // Through the apply functions, not the setters: the setters early-out on
    // unchanged values, and the values have not changed — the freshly built
    // channels simply have not heard them yet.
    applyDetectorCorners();
    applyTimeConstants();

    // Unity, and snapped rather than faded: prepare is not a control movement.
    sendOutTrim.prepare (sampleRate, defaultSettlingSeconds, 1.0);
    returnInTrim.prepare (sampleRate, defaultSettlingSeconds, 1.0);
    returnOutTrim.prepare (sampleRate, defaultSettlingSeconds, 1.0);

    tape.prepare (sampleRate, maxBlock, numChannels);
    tapePointers.assign ((size_t) std::max (1, numChannels), nullptr);

    prepared = true;
}

void SrEngine::reset()
{
    for (auto& c : channels)
    {
        c->skew.reset();
        c->skewInverse.reset();
        c->antiSat.reset();
        c->antiSatInverse.reset();

        c->hfEncode.reset();
        c->lfEncode.reset();
        c->hfDecode.reset();
        c->lfDecode.reset();

        c->modeDelay.reset();
    }

    tape.reset();
    lastDominants.store (0, std::memory_order_relaxed);

    sendOutTrim.snapTo (sendOutTrim.getTarget());
    returnInTrim.snapTo (returnInTrim.getTarget());
    returnOutTrim.snapTo (returnOutTrim.getTarget());
}

void SrEngine::setSpectralSmoothing (double bins)
{
    smoothing = bins;

    for (auto& c : channels)
    {
        c->hfEncode.setSpectralSmoothingBins (bins);
        c->lfEncode.setSpectralSmoothingBins (bins);
        c->hfDecode.setSpectralSmoothingBins (bins);
        c->lfDecode.setSpectralSmoothingBins (bins);
    }
}

void SrEngine::setSections (int count)
{
    // The same count to all four, unclamped by anything either half knows about its own
    // bins — see the note at the top of SrSideChain.h. Quantising the two resolutions
    // differently is exactly what would break the complementarity.
    for (auto& c : channels)
    {
        c->hfEncode.setSections (count);
        c->lfEncode.setSections (count);
        c->hfDecode.setSections (count);
        c->lfDecode.setSections (count);
    }

    // Stored as what the halves actually rounded to rather than what was asked for,
    // so the engine reports the truth back to the panel.
    if (! channels.empty())
        count = channels.front()->hfEncode.getSections();

    sections.store (count, std::memory_order_relaxed);
}

void SrEngine::setEnvelopeSmoothing (double bark)
{
    envelopeBark = bark;

    for (auto& c : channels)
    {
        c->hfEncode.setEnvelopeSmoothingBark (bark);
        c->lfEncode.setEnvelopeSmoothingBark (bark);
        c->hfDecode.setEnvelopeSmoothingBark (bark);
        c->lfDecode.setEnvelopeSmoothingBark (bark);
    }
}

void SrEngine::setSidechainHighPass (double cornerHz)
{
    // The engine-level guard matches the chain-level ones: the wrapper calls
    // this every chunk, and an unmoved corner should cost a compare, not a
    // walk of the channel list. prepare() re-applies through the apply
    // functions below, which a guard here would otherwise defeat.
    if (cornerHz == sidechainHighPassHz)
        return;

    sidechainHighPassHz = cornerHz;
    applyDetectorCorners();
}

void SrEngine::setSidechainLowPass (double cornerHz)
{
    if (cornerHz == sidechainLowPassHz)
        return;

    sidechainLowPassHz = cornerHz;
    applyDetectorCorners();
}

void SrEngine::applyDetectorCorners()
{
    // The same corners to all four — an encode and a decode that disagreed on
    // what the detection reads would compute different gains, and the exactness
    // of the round trip rests on them never doing that. The chain setters
    // early-out per corner, so pushing both here costs nothing extra.
    for (auto& c : channels)
    {
        c->hfEncode.setDetectorHighPass (sidechainHighPassHz);
        c->lfEncode.setDetectorHighPass (sidechainHighPassHz);
        c->hfDecode.setDetectorHighPass (sidechainHighPassHz);
        c->lfDecode.setDetectorHighPass (sidechainHighPassHz);

        c->hfEncode.setDetectorLowPass (sidechainLowPassHz);
        c->lfEncode.setDetectorLowPass (sidechainLowPassHz);
        c->hfDecode.setDetectorLowPass (sidechainLowPassHz);
        c->lfDecode.setDetectorLowPass (sidechainLowPassHz);
    }
}

void SrEngine::setTimeConstants (double upAttack, double upRelease,
                                 double downAttack, double downRelease)
{
    if (upAttack == upAttackMs && upRelease == upReleaseMs
        && downAttack == downAttackMs && downRelease == downReleaseMs)
        return;

    upAttackMs = upAttack;
    upReleaseMs = upRelease;
    downAttackMs = downAttack;
    downReleaseMs = downRelease;

    applyTimeConstants();
}

void SrEngine::applyTimeConstants()
{
    // Each direction's pair goes to both of its side chains — the hf and lf
    // halves are one detection and must keep moving together. Up and down may
    // differ; the price, stated plainly, is the exact null: encode and decode
    // that run different ballistics compute different gains, and the round trip
    // is complementary again only once the pairs agree.
    for (auto& c : channels)
    {
        c->hfEncode.setTimeConstants (upAttackMs, upReleaseMs);
        c->lfEncode.setTimeConstants (upAttackMs, upReleaseMs);
        c->hfDecode.setTimeConstants (downAttackMs, downReleaseMs);
        c->lfDecode.setTimeConstants (downAttackMs, downReleaseMs);
    }
}

void SrEngine::setLoopTrims (double sendOutGain, double returnInGain,
                             double returnOutGain, bool immediately) noexcept
{
    if (immediately)
    {
        sendOutTrim.snapTo (sendOutGain);
        returnInTrim.snapTo (returnInGain);
        returnOutTrim.snapTo (returnOutGain);
        return;
    }

    sendOutTrim.setTarget (sendOutGain);
    returnInTrim.setTarget (returnInGain);
    returnOutTrim.setTarget (returnOutGain);
}

void SrEngine::setProbesActive (bool shouldBeActive) noexcept
{
    // Straight onto the engine-owned probes rather than through the channel list:
    // this is called from the message thread when an editor opens or closes, and
    // prepare() may be rebuilding the channels at that very moment. The probes
    // themselves never move, so there is nothing here to race with.
    hfDisplayProbe.setActive (shouldBeActive);
    lfDisplayProbe.setActive (shouldBeActive);

    probesActive.store (shouldBeActive, std::memory_order_relaxed);
}

SpectrumProbe* SrEngine::getHfProbe (int channel) noexcept
{
    return channel == 0 ? &hfDisplayProbe : nullptr;
}

SpectrumProbe* SrEngine::getLfProbe (int channel) noexcept
{
    return channel == 0 ? &lfDisplayProbe : nullptr;
}

void SrEngine::setGainRealisation (GainRealisation r)
{
    realisation = r;

    for (auto& c : channels)
    {
        c->hfEncode.setGainRealisation (r);
        c->lfEncode.setGainRealisation (r);
        c->hfDecode.setGainRealisation (r);
        c->lfDecode.setGainRealisation (r);
    }
}

//==============================================================================
void SrEngine::encodePass (Channel& c, double* buffer, int numSamples)
{
    c.skew.process (buffer, numSamples);
    c.hfEncode.process (buffer, buffer, numSamples);
    c.lfEncode.process (buffer, buffer, numSamples);
    c.antiSat.process (buffer, numSamples);

    // Accumulated across the channels of a block — process() zeroes it first — so
    // the reported figure is the instance's, not just the last channel's.
    lastDominants.fetch_add (c.hfEncode.getNumDominants() + c.lfEncode.getNumDominants(),
                             std::memory_order_relaxed);
}

void SrEngine::decodePass (Channel& c, double* buffer, int numSamples)
{
    // Exactly the encode order reversed, so each inverse meets the state its
    // forward counterpart built.
    c.antiSatInverse.process (buffer, numSamples);
    c.lfDecode.process (buffer, buffer, numSamples);
    c.hfDecode.process (buffer, buffer, numSamples);
    c.skewInverse.process (buffer, numSamples);
}

//==============================================================================
void SrEngine::process (double* const* channelData, int numChannels, int numSamples)
{
    if (! prepared || numSamples <= 0)
        return;

    const int usable = std::min (numChannels, (int) channels.size());

    if (usable <= 0)
        return;

    // Switched off still has to honour the reported latency, or engaging it would
    // shift the signal in time relative to the rest of the session.
    if (processType != SrProcess::spectralRecording)
    {
        // In loop mode the modelled channel stays in circuit whatever the process is set
        // to, because in loop mode the channel is the thing being demonstrated and the
        // process is the thing being toggled.
        //
        // It did not, and that made the panel's central comparison say nothing: with the
        // process off the engine was a plain delay, so `off` was a clean signal and `b` was
        // a clean signal, and switching between them was inaudible. What the user is meant
        // to hear is the channel with the noise reduction against the channel without it.
        // The trims come too, so the level does not jump with the switch.
        if (mode == SrMode::loop)
        {
            applyLoopTrim (sendOutTrim, channelData, usable, numSamples);

            for (int ch = 0; ch < usable; ++ch)
                tapePointers[(size_t) ch] = channelData[ch];

            tape.process (tapePointers.data(), usable, numSamples);

            applyLoopTrim (returnInTrim, channelData, usable, numSamples);
            applyLoopTrim (returnOutTrim, channelData, usable, numSamples);
        }

        const int delay = getLatencySamples();

        for (int ch = 0; ch < usable; ++ch)
        {
            auto& c = *channels[(size_t) ch];

            c.modeDelay.setDelay (delay);
            c.modeDelay.process (channelData[ch], numSamples);
        }

        return;
    }

    lastDominants.store (0, std::memory_order_relaxed);

    // The display probes follow whichever direction is running: the encode boost
    // in encode and loop modes, the decode cut in decode mode — a decode-only
    // instance used to publish nothing and the plot froze on its last frame.
    // Re-wired per block, because setMode is a flag flip on the audio thread, not
    // a re-prepare; these are four pointer stores.
    {
        auto& c0 = *channels.front();
        const bool showDecode = mode == SrMode::decode;

        c0.hfEncode.setProbe (showDecode ? nullptr : &hfDisplayProbe);
        c0.lfEncode.setProbe (showDecode ? nullptr : &lfDisplayProbe);
        c0.hfDecode.setProbe (showDecode ? &hfDisplayProbe : nullptr);
        c0.lfDecode.setProbe (showDecode ? &lfDisplayProbe : nullptr);
    }

    for (int ch = 0; ch < usable; ++ch)
    {
        auto& c = *channels[(size_t) ch];

        if (mode == SrMode::encode || mode == SrMode::loop)
            encodePass (c, channelData[ch], numSamples);
    }

    // The impairment the noise reduction exists to defeat sits between the two
    // passes, which is the whole point of loop mode. The send out and return in trims sit
    // either side of it, exactly where they do on the hardware: send out scales what is
    // printed, return in scales what comes back.
    if (mode == SrMode::loop)
    {
        applyLoopTrim (sendOutTrim, channelData, usable, numSamples);

        for (int ch = 0; ch < usable; ++ch)
            tapePointers[(size_t) ch] = channelData[ch];

        tape.process (tapePointers.data(), usable, numSamples);

        applyLoopTrim (returnInTrim, channelData, usable, numSamples);
    }

    for (int ch = 0; ch < usable; ++ch)
    {
        auto& c = *channels[(size_t) ch];

        const bool decoding = (mode == SrMode::decode || mode == SrMode::loop);

        if (decoding)
            decodePass (c, channelData[ch], numSamples);

        // Every pass not taken has to be paid for in delay instead, or the plugin's
        // latency would change with a switch.
        const int passesTaken = ((mode == SrMode::encode || mode == SrMode::loop) ? 1 : 0)
                              + (decoding ? 1 : 0);
        const int padding = (2 - passesTaken) * geometry.passLatencySamples();

        c.modeDelay.setDelay (padding);
        c.modeDelay.process (channelData[ch], numSamples);
    }

    // Last in the loop, as it is on the hardware: the make-up after decoding.
    if (mode == SrMode::loop)
        applyLoopTrim (returnOutTrim, channelData, usable, numSamples);
}

} // namespace lime
