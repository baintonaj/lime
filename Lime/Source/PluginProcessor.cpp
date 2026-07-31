/*
  ==============================================================================

    Lime — plugin wrapper. All DSP lives under Source/dsp and knows nothing about
    JUCE, so it can be exercised by the standalone tests in Source/tests.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

namespace
{
    // Parameter identifiers. The control set mirrors the reference hardware front
    // panel; every internal constant is compiled in rather than exposed.
    constexpr auto processId = "process";
    constexpr auto modeId    = "mode";
    constexpr auto bypassId  = "bypass";
    constexpr auto recInId   = "recIn";
    constexpr auto recOutId  = "recOut";
    constexpr auto playInId  = "playIn";
    constexpr auto playOutId = "playOut";
    constexpr auto setUpId   = "setUp";
    constexpr auto inputGainId = "inputGain";
    constexpr auto outputGainId = "outputGain";
    constexpr auto mixId = "mix";
    constexpr auto polarityId = "polarity";
    constexpr auto autoGainId = "autoGain";
    constexpr auto scHpfId = "scHpf";
    constexpr auto scLpfId = "scLpf";
    constexpr auto scHpfOnId = "scHpfOn";
    constexpr auto scLpfOnId = "scLpfOn";
    constexpr auto upAttackId = "upAttack";
    constexpr auto upReleaseId = "upRelease";
    constexpr auto downAttackId = "downAttack";
    constexpr auto downReleaseId = "downRelease";

    /** The band split's resolution, compiled in rather than exposed.

        The control this replaces was measured to matter only at the bottom of its
        range, where it made the round trip audibly worse; from the default upward
        every setting sounded and nulled the same, which is the definition of a
        knob not worth the panel space. 2048 sections sits on the measured null
        plateau — the best figure the range produces — and stays comfortably below
        one section per bin at every rate. See SrSideChain::setSections.
    */
    constexpr int fixedSectionCount = 2048;

    /** The hardware's calibration trimmers cover roughly -10 to +10 dBr around reference
        level, which is the range quoted for matching studio lines. */
    juce::NormalisableRange<float> trimRange()
    {
        return { -10.0f, 10.0f, 0.01f };
    }

    /** The masters are not calibration trimmers and are not bound by the reference
        hardware's range: they exist to match the plugin to whatever is either side of it
        in a session, which can be a long way from reference level. */
    juce::NormalisableRange<float> masterGainRange()
    {
        return { -24.0f, 24.0f, 0.01f };
    }

    /** Zero is the dry signal, a hundred is the processed one, and the range carries on to
        two hundred rather than stopping there.

        Past a hundred the blend extrapolates instead of interpolating: the output becomes
        the processed signal plus the difference between it and the dry one again, so
        whatever the process did is applied twice over. It is the same arithmetic either
        side of a hundred — only the sign of what is being added changes — and it is the
        cheapest way to hear what a subtle setting is actually doing. */
    juce::NormalisableRange<float> mixRange()
    {
        return { 0.0f, 200.0f, 0.01f };
    }

    /** The side-chain high pass sweeps the whole audio band with its skew centred
        on the 80 Hz default — the corner a side-chain high pass classically rests
        at, keeping bass out of the detectors — so the resting position points
        straight up and every turn away from noon is a departure from stock. Far
        from a log taper, but a default a user can read at a glance across a
        session is worth more than symmetry. */
    juce::NormalisableRange<float> scHpfRange()
    {
        juce::NormalisableRange<float> range { 20.0f, 20000.0f, 1.0f };
        range.setSkewForCentre (80.0f);
        return range;
    }

    /** The matching low pass, same span, its skew centred on its own 6 kHz
        default for the same reason: noon is stock. */
    juce::NormalisableRange<float> scLpfRange()
    {
        juce::NormalisableRange<float> range { 20.0f, 20000.0f, 1.0f };
        range.setSkewForCentre (6000.0f);
        return range;
    }

    /** An attack runs 1 to 100 ms — the span settled by ear. The published
        fast trackers ahead of the governed smoother respond at 5 to 15 ms
        regardless, so the low end approaches them rather than diving a decade
        under; the floor also coincides with the a-type's 1 ms fast-attack
        cap, so no knob position asks for what the process cannot run. Skew
        centred on the 5 ms default: noon is stock. */
    juce::NormalisableRange<float> attackRange()
    {
        juce::NormalisableRange<float> range { 1.0f, 100.0f, 0.01f };
        range.setSkewForCentre (5.0f);
        return range;
    }

    /** A release spans 10 ms to 2.5 s around its 100 ms default — the
        neighbourhood of the published 80/150 ms finals — with the same rule:
        noon is stock. */
    juce::NormalisableRange<float> releaseRange()
    {
        juce::NormalisableRange<float> range { 10.0f, 2500.0f, 0.1f };
        range.setSkewForCentre (100.0f);
        return range;
    }

    /** Two decimals, because the interesting settings of a parallel blend are not the
        round numbers and a knob that reads "50 %" over a range of half a percent is
        telling the user less than it appears to. */
    juce::AudioParameterFloatAttributes mixAttributes()
    {
        return juce::AudioParameterFloatAttributes()
            .withLabel ("%")
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (value, 2) + " %";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("-+.0123456789").getFloatValue();
            });
    }

    /** Trims read out as a signed decibel offset to one place.

        The default float-parameter formatting gave "-0.00" at the centre of the range:
        two decimals nobody can set by hand, and a minus sign on a value that is zero,
        because the normalised round trip lands a hair below it. */
    juce::AudioParameterFloatAttributes trimAttributes()
    {
        return juce::AudioParameterFloatAttributes()
            .withLabel ("dB")
            .withStringFromValueFunction ([] (float value, int)
            {
                auto shown = std::round (value * 10.0f) / 10.0f;

                if (shown == 0.0f)
                    shown = 0.0f;      // clears a negative zero's sign bit

                // Carries its unit, as the footer's level readout does. A bare "0.0"
                // beside a knob could be a ratio, a percentage or a level.
                return (shown > 0.0f ? "+" : "") + juce::String (shown, 1) + " dB";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("-+.0123456789").getFloatValue();
            });
    }

    /** Frequencies read in the unit a musician says them: "6.00 kHz" above a
        thousand, whole hertz below. Typed values accept "6000", "6k" and "6 kHz"
        alike — the only 'k' a frequency box ever sees is a kilo prefix, so its
        presence anywhere in the text means thousands; the range clamps afterwards. */
    juce::AudioParameterFloatAttributes hzAttributes()
    {
        return juce::AudioParameterFloatAttributes()
            .withLabel ("Hz")
            .withStringFromValueFunction ([] (float value, int)
            {
                if (value >= 1000.0f)
                    return juce::String (value / 1000.0f, 2) + " kHz";

                return juce::String ((int) std::round (value)) + " Hz";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                const auto number = text.retainCharacters ("-+.0123456789").getFloatValue();
                return text.containsIgnoreCase ("k") ? number * 1000.0f : number;
            });
    }

    /** Times read in the unit a musician says them: "5.0 ms" below a second,
        "1.20 s" from there up. Typed values: a bare number is milliseconds, and
        text carrying an 's' with no 'm' — "1.2 s", "2s" — is seconds; the range
        clamps afterwards. */
    juce::AudioParameterFloatAttributes msAttributes()
    {
        return juce::AudioParameterFloatAttributes()
            .withLabel ("ms")
            .withStringFromValueFunction ([] (float value, int)
            {
                if (value >= 1000.0f)
                    return juce::String (value / 1000.0f, 2) + " s";

                return juce::String (value, 1) + " ms";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                const auto number = text.retainCharacters ("-+.0123456789").getFloatValue();
                const bool seconds = text.containsIgnoreCase ("s")
                                  && ! text.containsIgnoreCase ("m");

                return seconds ? number * 1000.0f : number;
            });
    }
}

//==============================================================================
juce::AudioProcessorValueTreeState::ParameterLayout LimeAudioProcessor::createParameterLayout()
{
    using namespace juce;

    AudioProcessorValueTreeState::ParameterLayout layout;

    // The three-position process switch: Off, the four-band A process, and the staggered
    // spectral B process.
    layout.add (std::make_unique<AudioParameterChoice> (
        ParameterID { processId, 1 }, "Process",
        StringArray { "Off", "A", "B" }, 2));

    // Up boosts quiet detail for the channel; Down reverses it; Both does the round trip
    // through the channel in one instance, so the noise reduction is audible.
    layout.add (std::make_unique<AudioParameterChoice> (
        ParameterID { modeId, 1 }, "Mode",
        StringArray { "Up", "Down", "Both" }, 2));

    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { bypassId, 1 }, "Bypass", false));

    // Named for what they do to the user, in the modes' own vocabulary: the in-trims
    // set how hard the programme drives each pass — how much of it falls into the
    // action region, which is the threshold control of this architecture — and the
    // out-trims are that pass's level compensation. The identifiers keep their
    // original hardware-derived spelling: they are invisible to the user, and
    // renaming one would silently drop its value out of every session already saved.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { recInId, 1 }, "Up In", trimRange(), 0.0f, trimAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { recOutId, 1 }, "Up Out", trimRange(), 0.0f, trimAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { playInId, 1 }, "Down In", trimRange(), 0.0f, trimAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { playOutId, 1 }, "Down Out", trimRange(), 0.0f, trimAttributes()));

    // Set Up sends the calibration signal and, on playback, runs Auto Compare.
    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { setUpId, 1 }, "Set Up", false));

    // The masters, which the reference hardware has no equivalent of because it did not
    // need one: it sat in a console at a known level. A plugin does not, so there is a
    // gain either side of the whole thing and a blend against the dry path.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { inputGainId, 1 }, "Input", masterGainRange(), 0.0f, trimAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { outputGainId, 1 }, "Output", masterGainRange(), 0.0f, trimAttributes()));

    // Loudness-matching makeup for the wet path, so bypass and mix compare
    // character rather than volume. A dynamics unit without makeup gain makes
    // every comparison a loudness test; see AutoGain.
    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { autoGainId, 1 }, "Auto Gain", false));

    // The side-chain high pass. It filters what the compressors' level detection
    // reads, never the audio: below the corner the programme reads quieter than
    // it is, so the process treats it as low-level material, and the applied
    // gain curve on the plot moves accordingly. Because encode and decode weight
    // their detection identically, any setting leaves the round trip exact —
    // agree on it across two instances, like the trims. 18 dB an octave,
    // Butterworth, so the corner reads true. Off the panel since the time
    // constants took the knobs, but still full parameters, the arrangement mix
    // and Set Up have — the filters themselves are unchanged.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { scHpfId, 1 }, "Sidechain HPF", scHpfRange(), 80.0f, hzAttributes()));

    // Its low-pass partner: the two corners band-limit what the detection hears.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { scLpfId, 1 }, "Sidechain LPF", scLpfRange(), 6000.0f, hzAttributes()));

    // Each filter has its own switch, and both rest off: a fresh instance detects
    // with the full band, and the knobs' corners are positions held in readiness
    // rather than filters already in circuit. Engaging one is deliberate — the
    // switch, not the knob leaving its default, is what puts it in the path.
    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { scHpfOnId, 1 }, "Sidechain HPF On", false));
    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { scLpfOnId, 1 }, "Sidechain LPF On", false));

    // The per-half time constants, always in circuit: each pair overrides the
    // dominant final smoothers of its own direction's detection — attack while
    // the detected level rises, release while it falls — leaving the fast
    // published trackers in place. The defaults are therefore the stock
    // ballistics of a fresh instance. The halves are deliberately independent,
    // and the trade is stated plainly: the round trip is exact only while up
    // and down run the same times, so matching the pairs is what preserves
    // the null, and splitting them is a creative choice that forfeits it.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { upAttackId, 1 }, "Up Attack", attackRange(), 5.0f, msAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { upReleaseId, 1 }, "Up Release", releaseRange(), 100.0f, msAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { downAttackId, 1 }, "Down Attack", attackRange(), 5.0f, msAttributes()));
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { downReleaseId, 1 }, "Down Release", releaseRange(), 100.0f, msAttributes()));

    // Off the panel since the side-chain low pass took its knob, but still a full
    // parameter — reachable from the host's list and automatable, the same
    // arrangement Set Up has. The blend arithmetic, and the guarantee that zero
    // is bit-exact dry, are unchanged.
    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { mixId, 1 }, "Mix", mixRange(), 100.0f, mixAttributes()));

    layout.add (std::make_unique<AudioParameterBool> (
        ParameterID { polarityId, 1 }, "Polarity", false));

    return layout;
}

//==============================================================================
LimeAudioProcessor::LimeAudioProcessor()
     : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
       parameters (*this, nullptr, "Lime", createParameterLayout())
{
    processParam = parameters.getRawParameterValue (processId);
    modeParam    = parameters.getRawParameterValue (modeId);
    bypassParam  = parameters.getRawParameterValue (bypassId);
    recInParam   = parameters.getRawParameterValue (recInId);
    recOutParam  = parameters.getRawParameterValue (recOutId);
    playInParam  = parameters.getRawParameterValue (playInId);
    playOutParam = parameters.getRawParameterValue (playOutId);
    setUpParam   = parameters.getRawParameterValue (setUpId);
    inputGainParam = parameters.getRawParameterValue (inputGainId);
    outputGainParam = parameters.getRawParameterValue (outputGainId);
    mixParam = parameters.getRawParameterValue (mixId);
    polarityParam = parameters.getRawParameterValue (polarityId);
    autoGainParam = parameters.getRawParameterValue (autoGainId);
    scHpfParam = parameters.getRawParameterValue (scHpfId);
    scLpfParam = parameters.getRawParameterValue (scLpfId);
    scHpfOnParam = parameters.getRawParameterValue (scHpfOnId);
    scLpfOnParam = parameters.getRawParameterValue (scLpfOnId);
    upAttackParam = parameters.getRawParameterValue (upAttackId);
    upReleaseParam = parameters.getRawParameterValue (upReleaseId);
    downAttackParam = parameters.getRawParameterValue (downAttackId);
    downReleaseParam = parameters.getRawParameterValue (downReleaseId);
}

LimeAudioProcessor::~LimeAudioProcessor() = default;

//==============================================================================
const juce::String LimeAudioProcessor::getName() const              { return JucePlugin_Name; }
bool LimeAudioProcessor::acceptsMidi() const                        { return false; }
bool LimeAudioProcessor::producesMidi() const                       { return false; }
bool LimeAudioProcessor::isMidiEffect() const                       { return false; }

double LimeAudioProcessor::getTailLengthSeconds() const
{
    // The engine's own delay is reported as latency, not tail. Beyond it the longest
    // ringing element is the decoder's de-skewing network, whose slowest pole sits at
    // 9.5 Hz — the zero of the 40 Hz spectral skewing shelf. With Butterworth damping
    // that decays with a time constant near 24 ms, so a quarter of a second covers it
    // to well below audibility.
    return 0.25;
}

int LimeAudioProcessor::getNumPrograms()                            { return 1; }
int LimeAudioProcessor::getCurrentProgram()                         { return 0; }
void LimeAudioProcessor::setCurrentProgram (int)                    {}
const juce::String LimeAudioProcessor::getProgramName (int)         { return {}; }
void LimeAudioProcessor::changeProgramName (int, const juce::String&) {}

//==============================================================================
void LimeAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    const int numChannels = juce::jmax (getTotalNumInputChannels(), getTotalNumOutputChannels());

    // Everything downstream is sized to this; processDouble splits anything longer
    // into pieces of it, so a host that breaks its promise costs a few extra passes
    // rather than an allocation on the audio thread.
    preparedMaxBlock = juce::jmax (1, samplesPerBlock);

    engine.prepare (sampleRate, preparedMaxBlock, numChannels);

    // The band split's resolution is a constant of the design now, pinned to the
    // measured null plateau; see fixedSectionCount.
    engine.setSections (fixedSectionCount);
    atype.prepare (sampleRate, preparedMaxBlock, numChannels);

    calibrationNoise.prepare (sampleRate);
    calibrationTone.prepare (sampleRate);
    referenceNoise.prepare (sampleRate);
    noiseDetector.prepare (sampleRate);
    autoCompare.prepare (sampleRate);

    calibrationBuffer.setSize (2, preparedMaxBlock, false, true, true);

    conversionBuffer.setSize (numChannels, preparedMaxBlock, false, true, true);
    channelPointers.assign ((size_t) juce::jmax (1, numChannels), nullptr);

    setLatencySamples (engine.getLatencySamples());

    // The meter's decay assumes a continuous stream; across a stop it would
    // otherwise open showing the last thing the previous run played.
    meterDb.store (-100.0f, std::memory_order_relaxed);

    // Engaged from the first sample rather than faded in: a host restoring a session is not
    // a control movement, and 75 ms of the wrong settings at the top of a bounce would be a
    // bug rather than a nicety.
    engagedProcessIndex = -1;
    engagedModeIndex = -1;
    configurationChanging = false;
    applyParameters (true);

    const auto trims = currentTrims();

    inTrim.prepare (sampleRate, controlSettlingSeconds, trims.in);
    outTrim.prepare (sampleRate, controlSettlingSeconds, trims.out);
    engine.setLoopTrims (trims.loopSendOut, trims.loopReturnIn, trims.loopReturnOut, true);

    inputGain.prepare (sampleRate, controlSettlingSeconds,
                       juce::Decibels::decibelsToGain (inputGainParam->load()));
    outputGain.prepare (sampleRate, controlSettlingSeconds, currentOutputGain());

    bypassAmount.prepare (sampleRate, controlSettlingSeconds, *bypassParam > 0.5f ? 0.0 : 1.0);
    mixAmount.prepare (sampleRate, controlSettlingSeconds, mixParam->load() / 100.0);
    transition.prepare (sampleRate, transitionSeconds, 1.0);

    calibrationAmount.prepare (sampleRate, controlSettlingSeconds, 0.0);

    autoGain.prepare (sampleRate);
    autoGain.setEnabled (*autoGainParam > 0.5f);

    // The dry path is delayed to sit exactly under the process's output, so the crossfade
    // compares like with like.
    dryDelays.resize ((size_t) juce::jmax (1, numChannels));

    for (auto& delay : dryDelays)
    {
        delay.prepare (engine.getLatencySamples());
        delay.setDelay (engine.getLatencySamples());
        delay.reset();
    }

    dryBuffer.setSize (numChannels, preparedMaxBlock, false, true, true);
}

void LimeAudioProcessor::releaseResources()
{
    // State is cleared; allocations are kept. Freeing the buffers here would turn a
    // stray block from a host that calls once more after releasing — which happens —
    // into an allocation on the audio thread. The memory is a few blocks' worth and
    // the next prepareToPlay resizes it anyway.
    engine.reset();
    atype.reset();
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool LimeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto& out = layouts.getMainOutputChannelSet();

    if (out != juce::AudioChannelSet::mono() && out != juce::AudioChannelSet::stereo())
        return false;

    return out == layouts.getMainInputChannelSet();
}
#endif

//==============================================================================
LimeAudioProcessor::Trims LimeAudioProcessor::currentTrims() const
{
    const auto gain = [] (const std::atomic<float>* p)
    {
        return juce::Decibels::decibelsToGain ((double) p->load());
    };

    Trims trims;

    switch (engine.getMode())
    {
        case lime::SrMode::decode:
            trims.in = gain (playInParam);
            trims.out = gain (playOutParam);
            break;

        case lime::SrMode::loop:
            // The whole chain is inside one instance, so every trim has somewhere real to
            // be. Only send in is outside the engine; the other three sit around the
            // modelled channel and the decoder, which is where the hardware puts them and
            // where they were doing nothing at all before.
            trims.in = gain (recInParam);
            trims.loopSendOut = gain (recOutParam);
            trims.loopReturnIn = gain (playInParam);
            trims.loopReturnOut = gain (playOutParam);
            break;

        case lime::SrMode::encode:
        default:
            trims.in = gain (recInParam);
            trims.out = gain (recOutParam);
            break;
    }

    return trims;
}

double LimeAudioProcessor::currentOutputGain() const
{
    const auto gain = juce::Decibels::decibelsToGain ((double) outputGainParam->load());

    return (*polarityParam > 0.5f ? -1.0 : 1.0) * gain;
}

//==============================================================================
void LimeAudioProcessor::applyParameters (bool immediately)
{
    const auto processIndex = (int) processParam->load();
    const auto modeIndex = (int) modeParam->load();

    if (processIndex == engagedProcessIndex && modeIndex == engagedModeIndex)
    {
        configurationChanging = false;
        return;
    }

    // Neither can be interpolated: Off, A and B are different processes, not three
    // values of one thing, and the same goes for encode against decode. So the
    // switch is thrown while the wet path is silent — faded out to the delayed dry
    // signal, changed, then faded back. Fading through dry rather than through
    // silence keeps the programme present the whole way across.
    //
    // Bypass is its own reason for the wet path being silent, and a good one: if nothing is
    // being heard there is nothing to fade, so the switch is thrown at once rather than
    // making the user wait out a transition they cannot hear.
    const bool inaudible = bypassAmount.getCurrent() <= 1.0e-4;

    if (! immediately && ! inaudible && transition.getCurrent() > 1.0e-6)
    {
        configurationChanging = true;
        return;
    }

    switch (processIndex)
    {
        case 0:  engine.setProcess (lime::SrProcess::off);   break;
        case 1:  engine.setProcess (lime::SrProcess::atype); break;
        default: engine.setProcess (lime::SrProcess::spectralRecording); break;
    }

    switch (modeIndex)
    {
        case 0:  engine.setMode (lime::SrMode::encode); break;
        case 1:  engine.setMode (lime::SrMode::decode); break;
        default: engine.setMode (lime::SrMode::loop);   break;
    }

    // Everything the old configuration was holding goes with it.
    //
    // The engine's state belongs to the process it was built for: the overlap-add buffers
    // carry a window of the old process's output, the side chains carry its control
    // history, and the padding delay carries whatever it last had a use for. Left alone,
    // switching from B to A handed back a quarter of a second of audio the delay line had
    // been sitting on since the last time it was in circuit — B's, arriving under A. There
    // is nothing in any of it worth keeping across a change of process, and this happens
    // while the wet path is silent, so clearing it costs nothing audible.
    engine.reset();
    atype.reset();

    engagedProcessIndex = processIndex;
    engagedModeIndex = modeIndex;
    engagedProcess.store (processIndex, std::memory_order_relaxed);
    engagedMode.store (modeIndex, std::memory_order_relaxed);
    configurationChanging = false;
}

void LimeAudioProcessor::applySmoothedGain (lime::ControlSmoother& smoother,
                                            juce::AudioBuffer<double>& buffer,
                                            int numChannels)
{
    const int numSamples = buffer.getNumSamples();

    // Arrived: one multiply per sample by a constant, or nothing at all at unity.
    if (! smoother.isMoving())
    {
        const auto gain = smoother.getCurrent();

        if (gain != 1.0)
            for (int ch = 0; ch < numChannels; ++ch)
                buffer.applyGain (ch, 0, numSamples, gain);

        return;
    }

    // Moving: the smoother advances once per sample, and every channel is given the same
    // value, so a stereo pair cannot drift apart mid-ramp.
    for (int n = 0; n < numSamples; ++n)
    {
        const auto gain = smoother.next();

        for (int ch = 0; ch < numChannels; ++ch)
            buffer.getWritePointer (ch)[n] *= gain;
    }
}

//==============================================================================
void LimeAudioProcessor::processDouble (juce::AudioBuffer<double>& buffer)
{
    const int totalSamples = buffer.getNumSamples();

    if (totalSamples <= 0 || ! engine.isPrepared() || preparedMaxBlock <= 0)
        return;

    if (totalSamples <= preparedMaxBlock)
    {
        processChunk (buffer);
        return;
    }

    // The host broke its block-size promise — offline bounces in several hosts do.
    // Everything downstream was sized to the promise in prepareToPlay, so the block
    // is taken in promised-size pieces rather than grown for. The smoothers advance
    // per piece, which is at least as smooth as per block.
    const int numChannels = buffer.getNumChannels();
    double* pointers[32];

    if (numChannels > 32)
    {
        // Unreachable — the buses are mono or stereo — but if it ever fires,
        // silence is a safer thing to hand back than unprocessed input.
        buffer.clear();
        return;
    }

    for (int offset = 0; offset < totalSamples; offset += preparedMaxBlock)
    {
        const int length = juce::jmin (preparedMaxBlock, totalSamples - offset);

        for (int ch = 0; ch < numChannels; ++ch)
            pointers[ch] = buffer.getWritePointer (ch) + offset;

        // A view over the piece; with this few channels JUCE keeps the pointer
        // table in preallocated space, so constructing it does not allocate.
        juce::AudioBuffer<double> view (pointers, numChannels, length);
        processChunk (view);
    }
}

void LimeAudioProcessor::processChunk (juce::AudioBuffer<double>& buffer)
{
    const int numSamples = buffer.getNumSamples();
    const int numInputs = juce::jmin (buffer.getNumChannels(), getTotalNumInputChannels());
    const int numOutputs = juce::jmin (buffer.getNumChannels(), getTotalNumOutputChannels());

    // Clear any output channels the input did not supply.
    for (int ch = numInputs; ch < numOutputs; ++ch)
        buffer.clear (ch, 0, numSamples);

    const int numToProcess = juce::jmin (numInputs, numOutputs);

    if (numToProcess <= 0 || (int) dryDelays.size() < numToProcess)
        return;

    // A restored session's settings must be in force on the first sample, not faded
    // in over 1.5 s of whatever the previous patch was doing; see setStateInformation.
    // Acquire pairs with the release there, so every parameter write is visible.
    if (stateJustRestored.exchange (false, std::memory_order_acquire))
    {
        applyParameters (true);

        const auto restoredTrims = currentTrims();
        inTrim.snapTo (restoredTrims.in);
        outTrim.snapTo (restoredTrims.out);
        engine.setLoopTrims (restoredTrims.loopSendOut, restoredTrims.loopReturnIn,
                             restoredTrims.loopReturnOut, true);

        inputGain.snapTo (juce::Decibels::decibelsToGain ((double) inputGainParam->load()));
        outputGain.snapTo (currentOutputGain());
        bypassAmount.snapTo (isBypassed() ? 0.0 : 1.0);
        mixAmount.snapTo ((double) mixParam->load() / 100.0);
        transition.snapTo (1.0);
        autoGain.reset();
        autoGain.setEnabled (*autoGainParam > 0.5f);

        const bool restoredSendingCalibration = *setUpParam > 0.5f
            && engine.getMode() != lime::SrMode::decode
            && engine.getProcess() != lime::SrProcess::off;
        calibrationAmount.snapTo (restoredSendingCalibration ? 1.0 : 0.0);
    }

    // A change of process or mode is held back until the wet path has faded out; see
    // applyParameters.
    applyParameters (false);

    // The side-chain filters, to both engines every chunk: a switched-off filter
    // reads as no corner at all, which the engines take as unity weighting. Each
    // setter is a no-op unless its value moved, and needs no smoothing of its
    // own: it feeds level detection, whose own time constants — and the slew
    // limit on the composite transfer — are what pace the resulting gain
    // movement, switch throws included.
    const double hpfCorner = *scHpfOnParam > 0.5f ? (double) scHpfParam->load() : 0.0;
    const double lpfCorner = *scLpfOnParam > 0.5f ? (double) scLpfParam->load() : 0.0;

    engine.setSidechainHighPass (hpfCorner);
    engine.setSidechainLowPass (lpfCorner);
    atype.setSidechainHighPass (hpfCorner);
    atype.setSidechainLowPass (lpfCorner);

    // The per-half time constants, by the same arrangement — always in
    // circuit, the knobs' values are the ballistics. (The engines' zero
    // sentinel, meaning the published constants, remains a library capability
    // the wrapper no longer reaches for.) Each setter is a no-op unless a
    // value moved, and both engines read the same four loads.
    const double upAttack = (double) upAttackParam->load();
    const double upRelease = (double) upReleaseParam->load();
    const double downAttack = (double) downAttackParam->load();
    const double downRelease = (double) downReleaseParam->load();

    engine.setTimeConstants (upAttack, upRelease, downAttack, downRelease);
    atype.setTimeConstants (upAttack, upRelease, downAttack, downRelease);

    bypassAmount.setTarget (isBypassed() ? 0.0 : 1.0);
    transition.setTarget (configurationChanging ? 0.0 : 1.0);
    mixAmount.setTarget ((double) mixParam->load() / 100.0);

    // The master input gain, ahead of everything including the dry copy — so it is what
    // the whole plugin sees, bypass included, rather than a drive control on the process.
    inputGain.setTarget (juce::Decibels::decibelsToGain ((double) inputGainParam->load()));
    applySmoothedGain (inputGain, buffer, numToProcess);

    // The dry copy, taken before the input trim so a bypass is exactly the input, and
    // delayed to sit under the process's output. The engine keeps running while bypassed:
    // it costs the same as not bypassing, and it means coming out of bypass finds the side
    // chain already tracking the programme instead of starting cold.
    //
    // Sized in prepareToPlay; processDouble never hands this more than the promise.
    jassert (dryBuffer.getNumChannels() >= numToProcess && dryBuffer.getNumSamples() >= numSamples);

    if (dryBuffer.getNumChannels() < numToProcess || dryBuffer.getNumSamples() < numSamples)
        return;

    for (int ch = 0; ch < numToProcess; ++ch)
    {
        auto* dry = dryBuffer.getWritePointer (ch);
        juce::FloatVectorOperations::copy (dry, buffer.getReadPointer (ch), numSamples);
        dryDelays[(size_t) ch].process (dry, numSamples);
    }

    // The calibration trims sit either side of the process, exactly as the hardware's in and
    // out trimmers do. In loop mode three of the four fall inside the engine, between the
    // encoder, the modelled channel and the decoder; see currentTrims. All are smoothed, so
    // any of them can be swept while audio is running.
    const auto trims = currentTrims();

    inTrim.setTarget (trims.in);
    outTrim.setTarget (trims.out);
    engine.setLoopTrims (trims.loopSendOut, trims.loopReturnIn, trims.loopReturnOut);

    applySmoothedGain (inTrim, buffer, numToProcess);

    jassert ((int) channelPointers.size() >= numToProcess);   // sized in prepareToPlay

    if ((int) channelPointers.size() < numToProcess)
        return;

    for (int ch = 0; ch < numToProcess; ++ch)
        channelPointers[(size_t) ch] = buffer.getWritePointer (ch);

    const auto processType = engine.getProcess();
    const bool setUp = *setUpParam > 0.5f;
    const bool encoding = engine.getMode() != lime::SrMode::decode;

    // Set Up, following section 4.3. While recording, the calibration signal replaces the
    // programme: shaped noise at 15 dB below reference for B, a tone at reference for A. The
    // 15 dB is deliberate — at reference level the noise would risk saturating the extremes,
    // so what came back would not reflect the recorder's true response.
    //
    // Faded in over the programme rather than cut to, so engaging Set Up down a live line
    // does not put a step into it.
    const bool sendingCalibration = setUp && encoding && processType != lime::SrProcess::off;

    calibrationAmount.setTarget (sendingCalibration ? 1.0 : 0.0);

    if (sendingCalibration || calibrationAmount.getCurrent() > 1.0e-4)
    {
        jassert (calibrationBuffer.getNumSamples() >= numSamples);   // sized in prepareToPlay

        auto* generated = calibrationBuffer.getWritePointer (0);

        if (processType == lime::SrProcess::atype)
            calibrationTone.process (generated, numSamples);
        else
            calibrationNoise.process (generated, numSamples);

        auto* const* io = buffer.getArrayOfWritePointers();

        for (int n = 0; n < numSamples; ++n)
        {
            const auto amount = calibrationAmount.next();

            for (int ch = 0; ch < numToProcess; ++ch)
                io[ch][n] += amount * (generated[n] - io[ch][n]);
        }
    }

    if (processType == lime::SrProcess::atype)
    {
        atype.setMode (encoding ? lime::AtypeMode::encode : lime::AtypeMode::decode);
        atype.process (channelPointers.data(), numToProcess, numSamples);

        // A-type is a time-domain process with no latency of its own, so it has to be
        // padded to the latency the plugin reports for every mode.
        engine.setProcess (lime::SrProcess::off);
        engine.process (channelPointers.data(), numToProcess, numSamples);
        engine.setProcess (processType);
    }
    else
    {
        engine.process (channelPointers.data(), numToProcess, numSamples);
    }

    // Auto Compare, section 4.3: on playback, once the calibration noise is detected,
    // alternate the tape return with an internally generated uninterrupted reference in four
    // second segments. The reference has no nicks, which is how they are told apart.
    if (setUp && ! encoding && processType == lime::SrProcess::spectralRecording)
    {
        noiseDetector.push (buffer.getReadPointer (0), numSamples);
        autoCompare.setRunning (noiseDetector.isCalibrationNoisePresent());

        jassert (calibrationBuffer.getNumSamples() >= numSamples);   // sized in prepareToPlay

        auto* reference = calibrationBuffer.getWritePointer (0);
        referenceNoise.process (reference, numSamples);

        // One call across all channels: the segment clock advances per sample, not
        // per channel, so a stereo pair switches source on the same sample.
        autoCompare.process (buffer.getArrayOfWritePointers(), reference,
                             buffer.getArrayOfWritePointers(), numToProcess, numSamples);
    }

    applySmoothedGain (outTrim, buffer, numToProcess);

    // Loudness-matching makeup, measured against the delayed dry copy and applied
    // to the wet path before the blend — so bypass, mix and the process switches
    // compare character at matched loudness rather than louder-wins. The
    // measurement deliberately runs while the switch is off: engaging then finds
    // a warm 3-second estimate instead of ramping in from a cold one.
    autoGain.setEnabled (*autoGainParam > 0.5f);
    autoGain.measure (dryBuffer.getArrayOfReadPointers(),
                      buffer.getArrayOfReadPointers(), numToProcess, numSamples);
    autoGain.apply (buffer.getArrayOfWritePointers(), numToProcess, numSamples);

    // Blend against the delayed dry path. This is what bypass is, what covers a change of
    // process or mode, and what the mix control sets: at zero the output is the input
    // delayed, exactly. The three multiply, so bypass still wins outright whatever mix is
    // set to, and a switch still fades cleanly through dry at any blend.
    //
    // Past unity the same line extrapolates rather than interpolates, which is how mix
    // reaches 200%. Nothing here has to know that: dry + a(wet - dry) walks past wet on its
    // own when a exceeds one.
    //
    // All three have to be advanced together or not at all. Advancing one in a block that
    // skips the others would let them slide out of step with each other, so the test asks
    // whether *any* of them has something left to do.
    const double settled = bypassAmount.getCurrent() * transition.getCurrent()
                         * mixAmount.getCurrent();

    const bool blending = bypassAmount.isMoving() || transition.isMoving() || mixAmount.isMoving()
                       || std::abs (settled - 1.0) > 1.0e-9;

    if (blending)
    {
        auto* const* wet = buffer.getArrayOfWritePointers();
        auto* const* dry = dryBuffer.getArrayOfWritePointers();

        for (int n = 0; n < numSamples; ++n)
        {
            const auto amount = bypassAmount.next() * transition.next() * mixAmount.next();

            for (int ch = 0; ch < numToProcess; ++ch)
                wet[ch][n] = dry[ch][n] + amount * (wet[ch][n] - dry[ch][n]);
        }
    }

    // The master output gain, carrying the polarity switch in its sign. Last, so the meter
    // and the level readout below report what actually leaves.
    outputGain.setTarget (currentOutputGain());
    applySmoothedGain (outputGain, buffer, numToProcess);

    // Peak across the block for the meter, decaying between blocks so it reads rather
    // than flickers.
    double peak = 0.0;

    for (int ch = 0; ch < numToProcess; ++ch)
        peak = juce::jmax (peak, (double) buffer.getMagnitude (ch, 0, numSamples));

    const float blockDb = peak > 0.0 ? (float) juce::Decibels::gainToDecibels (peak) : -100.0f;
    const float previous = meterDb.load (std::memory_order_relaxed);
    const float decayPerBlock = 12.0f * (float) numSamples / (float) juce::jmax (1.0, getSampleRate());

    meterDb.store (juce::jmax (blockDb, previous - decayPerBlock), std::memory_order_relaxed);
}

void LimeAudioProcessor::processBlock (juce::AudioBuffer<double>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    hostBypassRequested.store (false, std::memory_order_relaxed);
    processDouble (buffer);
}

void LimeAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    hostBypassRequested.store (false, std::memory_order_relaxed);
    processFloat (buffer);
}

//==============================================================================
juce::AudioProcessorParameter* LimeAudioProcessor::getBypassParameter() const
{
    return parameters.getParameter (bypassId);
}

void LimeAudioProcessor::processBlockBypassed (juce::AudioBuffer<double>& buffer, juce::MidiBuffer&)
{
    // In the standard AU/VST3 wrappers this never runs: getBypassParameter() is
    // non-null, so the host drives the bypass parameter and the normal path's
    // crossfade handles it. The override exists for embedded hosts — an
    // AudioProcessorGraph, a JUCE host — that call processBlockBypassed directly;
    // it runs the same path with the crossfade forced to the delayed dry, so the
    // reported latency stays honoured and the engine stays warm. The flag is
    // scoped to the call so it cannot latch if the host stops processing.
    juce::ScopedNoDenormals noDenormals;
    hostBypassRequested.store (true, std::memory_order_relaxed);
    processDouble (buffer);
    hostBypassRequested.store (false, std::memory_order_relaxed);
}

void LimeAudioProcessor::processBlockBypassed (juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    hostBypassRequested.store (true, std::memory_order_relaxed);
    processFloat (buffer);
    hostBypassRequested.store (false, std::memory_order_relaxed);
}

void LimeAudioProcessor::processFloat (juce::AudioBuffer<float>& buffer)
{
    const int numChannels = juce::jmin (buffer.getNumChannels(),
                                        conversionBuffer.getNumChannels());
    const int totalSamples = buffer.getNumSamples();

    if (totalSamples <= 0 || numChannels <= 0 || ! engine.isPrepared() || preparedMaxBlock <= 0)
        return;

    // Widen into 64-bit, process, narrow back — in promised-size pieces, so the
    // persistent conversion buffer sized in prepareToPlay never has to grow on
    // the audio thread however long a block the host sends.
    for (int offset = 0; offset < totalSamples; offset += preparedMaxBlock)
    {
        const int length = juce::jmin (preparedMaxBlock, totalSamples - offset);

        juce::AudioBuffer<double> view (conversionBuffer.getArrayOfWritePointers(),
                                        numChannels, length);

        for (int ch = 0; ch < numChannels; ++ch)
        {
            const auto* src = buffer.getReadPointer (ch) + offset;
            auto* dst = view.getWritePointer (ch);

            for (int n = 0; n < length; ++n)
                dst[n] = (double) src[n];
        }

        processChunk (view);

        for (int ch = 0; ch < numChannels; ++ch)
        {
            const auto* src = view.getReadPointer (ch);
            auto* dst = buffer.getWritePointer (ch) + offset;

            for (int n = 0; n < length; ++n)
                dst[n] = (float) src[n];
        }
    }
}

//==============================================================================
bool LimeAudioProcessor::hasEditor() const { return true; }

juce::AudioProcessorEditor* LimeAudioProcessor::createEditor()
{
    return new LimeAudioProcessorEditor (*this);
}

//==============================================================================
void LimeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    if (auto state = parameters.copyState(); auto xml = state.createXml())
    {
        // Stamped so a future change of parameter layout can migrate old sessions
        // deliberately instead of guessing from their shape.
        xml->setAttribute ("stateVersion", 1);
        copyXmlToBinary (*xml, destData);
    }
}

void LimeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    const auto xml = getXmlFromBinary (data, sizeInBytes);

    if (xml == nullptr || ! xml->hasTagName (parameters.state.getType()))
        return;

    // Version 1 is the current layout; states saved before the stamp existed carry
    // the same parameter set and load identically, so nothing branches on it yet.
    //
    // Defaults first, then the saved values over them. replaceState leaves any
    // parameter the incoming tree does not name at its *current* value, so loading
    // a state saved before a control existed would silently keep whatever the
    // previous patch had set that control to. The defaults are grafted into the
    // tree before the single replaceState rather than notified in afterwards:
    // setValueNotifyingHost fires host listeners, and a host writing automation
    // would record a restore as if it were a hand on the controls.
    auto restored = juce::ValueTree::fromXml (*xml);

    for (auto* p : getParameters())
    {
        if (auto* ranged = dynamic_cast<juce::RangedAudioParameter*> (p))
        {
            if (! restored.getChildWithProperty ("id", ranged->paramID).isValid())
            {
                juce::ValueTree param ("PARAM");
                param.setProperty ("id", ranged->paramID, nullptr);
                param.setProperty ("value",
                                   (double) ranged->convertFrom0to1 (ranged->getDefaultValue()),
                                   nullptr);
                restored.appendChild (param, nullptr);
            }
        }
    }

    parameters.replaceState (restored);

    // The next audio block throws the switches and snaps the smoothers: a restored
    // session is not a control movement. Release, paired with the acquire in
    // processChunk, so the flag publishes the parameter writes above. Until a block
    // runs the panel may briefly show the previous patch's engaged process on the
    // plot — with the transport stopped there is no safe way to reach into the
    // engine from here.
    stateJustRestored.store (true, std::memory_order_release);
}

//==============================================================================
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new LimeAudioProcessor();
}
