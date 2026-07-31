/*
  ==============================================================================

    Lime — a companding noise reduction system reinterpreted as a per-bin STFT process.

    The whole signal path is 64-bit. processBlock (AudioBuffer<double>&) is the
    real implementation; the float overload converts up into a persistent double
    buffer, runs the same engine, and converts back, so the process behaves
    identically whatever precision the host asks for.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "dsp/AtypeEngine.h"
#include "dsp/AutoGain.h"
#include "dsp/Calibration.h"
#include "dsp/ControlSmoother.h"
#include "dsp/Crossfade.h"
#include "dsp/DelayLine64.h"
#include "dsp/SrEngine.h"

//==============================================================================
class LimeAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    LimeAudioProcessor();
    ~LimeAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    //==============================================================================
    // Double precision is the native path; float is supported by conversion.
    bool supportsDoublePrecisionProcessing() const override    { return true; }

    void processBlock (juce::AudioBuffer<double>&, juce::MidiBuffer&) override;
    void processBlock (juce::AudioBuffer<float>&,  juce::MidiBuffer&) override;

    /** The host's bypass drives the plugin's own bypass parameter rather than a
        synthesised zero-latency one. Left to its default, the JUCE wrapper invents a
        bypass that passes audio with no delay while the host still applies this
        plugin's 200-odd milliseconds of compensation to everything else — the
        bypassed signal jumps forward in time. The panel's switch and the host's
        control are also the same parameter this way, so they cannot disagree. */
    juce::AudioProcessorParameter* getBypassParameter() const override;

    /** For embedded hosts only — the standard wrappers never call these once
        getBypassParameter() is non-null. Bypassed processing runs the normal path
        with the crossfade forced to dry, so the reported latency stays honoured
        and the engine keeps tracking the programme. */
    void processBlockBypassed (juce::AudioBuffer<double>&, juce::MidiBuffer&) override;
    void processBlockBypassed (juce::AudioBuffer<float>&,  juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    //==============================================================================
    juce::AudioProcessorValueTreeState& getState() noexcept { return parameters; }

    /** Reported latency in samples, or 0 before prepareToPlay. */
    int getEngineLatencySamples() const noexcept { return engine.getLatencySamples(); }

    lime::SrEngine& getEngine() noexcept { return engine; }

    /** Which source Auto Compare is currently passing, for the front-panel
        indicator: red for the internal reference, green for the tape return. */
    bool isAutoCompareShowingReference() const noexcept
    {
        return autoCompare.getCurrentSource() == lime::AutoCompareSource::reference;
    }

    bool isSetUpEngaged() const noexcept { return setUpParam != nullptr && *setUpParam > 0.5f; }

    bool isBypassed() const noexcept
    {
        return (bypassParam != nullptr && *bypassParam > 0.5f)
            || hostBypassRequested.load (std::memory_order_relaxed);
    }

    /** Which process is actually in circuit, as against which the parameter asks for.

        The plot follows this rather than the parameter, so it changes when the audio does
        rather than most of a second before it. Written on the audio thread, read by the editor. */
    int getEngagedProcessIndex() const noexcept
    {
        return engagedProcess.load (std::memory_order_relaxed);
    }

    /** Which mode is actually in circuit, on the same terms as the process above:
        the plot signs the fixed main path's contribution by it. */
    int getEngagedModeIndex() const noexcept
    {
        return engagedMode.load (std::memory_order_relaxed);
    }

    /** Peak output level in dB for the editor's meter. Written on the audio thread and
        read on the message thread, which for a single value a display samples is what
        a relaxed atomic is for. */
    float getMeterDb() const noexcept { return meterDb.load (std::memory_order_relaxed); }

    /** Every control reaches the signal path through a one-pole, so nothing a hand or a
        host automation lane does puts a step into the audio. 95% of any change lands
        inside this time — see ControlSmoother for why it is quoted that way. */
    static constexpr double controlSettlingSeconds = lime::defaultSettlingSeconds;

    /** How long a change of process or mode takes, one direction.

        Those two switches cannot be interpolated, so the wet path is faded out to the
        delayed dry, the switch is thrown while it is silent, and it is faded back — 1500 ms
        across the whole move. That is a great deal slower than a control smoother, and
        deliberately: a switch is not a knob, and at the speed a knob wants the change
        arrives as an event rather than as a transition. See lime::Crossfade for why it is
        a shaped ramp rather than a one-pole. */
    static constexpr double transitionSeconds = 0.75;

private:
    //==============================================================================
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

    /** Pulls the current parameter values into the engine.

        @param immediately  true from prepareToPlay, where a saved session's settings must
                            be in force on the first sample rather than faded in. Otherwise
                            a change of process or mode waits for the wet path to fade out,
                            since those two settings are different processes rather than two
                            values of one and cannot be interpolated. */
    void applyParameters (bool immediately);

    /** The shared 64-bit implementation both processBlock overloads funnel into.
        Splits any block longer than the size prepareToPlay was promised into
        promised-size pieces, so nothing downstream ever has to grow a buffer on
        the audio thread — hosts do break the promise, notably in offline bounces. */
    void processDouble (juce::AudioBuffer<double>& buffer);

    /** One promised-size piece of the block. This is the whole signal path; every
        buffer it touches was sized in prepareToPlay. */
    void processChunk (juce::AudioBuffer<double>& buffer);

    /** The float path: widen, process, narrow, chunk by chunk. */
    void processFloat (juce::AudioBuffer<float>& buffer);

    /** Applies a smoothed gain across a block, advancing the smoother once per sample so
        every channel follows the same trajectory. */
    void applySmoothedGain (lime::ControlSmoother&, juce::AudioBuffer<double>&, int numChannels);

    /** Which of the four calibration trims is in circuit where, for the current mode.

        The four come in pairs on the hardware and only the pair belonging to the direction
        is live — except in loop mode, where the channel is inside the instance and all four
        have a real position: send in and send out either side of the encoder, return in and
        return out either side of the decoder. The three that fall inside the loop are
        applied by the engine; see SrEngine::setLoopTrims. */
    struct Trims
    {
        double in = 1.0, out = 1.0;
        double loopSendOut = 1.0, loopReturnIn = 1.0, loopReturnOut = 1.0;
    };

    Trims currentTrims() const;

    /** The output master, signed by the polarity switch.

        Folded into one number deliberately: ControlSmoother is linear in gain, so a flip
        ramps from +g through zero to -g rather than stepping, and a polarity invert that
        does not step is the only kind worth having. */
    double currentOutputGain() const;

    juce::AudioProcessorValueTreeState parameters;

    std::atomic<float>* processParam = nullptr;
    std::atomic<float>* modeParam = nullptr;
    std::atomic<float>* bypassParam = nullptr;
    std::atomic<float>* recInParam = nullptr;
    std::atomic<float>* recOutParam = nullptr;
    std::atomic<float>* playInParam = nullptr;
    std::atomic<float>* playOutParam = nullptr;
    std::atomic<float>* setUpParam = nullptr;
    std::atomic<float>* inputGainParam = nullptr;
    std::atomic<float>* outputGainParam = nullptr;
    std::atomic<float>* mixParam = nullptr;
    std::atomic<float>* polarityParam = nullptr;
    std::atomic<float>* autoGainParam = nullptr;
    std::atomic<float>* scHpfParam = nullptr;
    std::atomic<float>* scLpfParam = nullptr;
    std::atomic<float>* scHpfOnParam = nullptr;
    std::atomic<float>* scLpfOnParam = nullptr;
    std::atomic<float>* upAttackParam = nullptr;
    std::atomic<float>* upReleaseParam = nullptr;
    std::atomic<float>* downAttackParam = nullptr;
    std::atomic<float>* downReleaseParam = nullptr;
    std::atomic<float>* upTimeOnParam = nullptr;
    std::atomic<float>* downTimeOnParam = nullptr;

    lime::SrEngine engine;
    lime::AtypeEngine atype;

    // Loudness-matching makeup on the processed path, measured and applied before
    // the blend, so every comparison the panel offers is level-matched.
    lime::AutoGain autoGain;

    // Calibration, per the reference hardware: calibration noise for SR, calibration tone for A-type, and
    // Auto Compare to hear the tape return against an internal reference.
    lime::CalibrationNoiseGenerator calibrationNoise;
    lime::CalibrationToneGenerator calibrationTone;
    lime::PinkNoiseGenerator referenceNoise;
    lime::CalibrationNoiseDetector noiseDetector;
    lime::AutoCompare autoCompare;

    juce::AudioBuffer<double> calibrationBuffer;

    // The smoothed controls, all linear gains. inTrim and outTrim are the calibration
    // pair in circuit for the current mode; inputGain and outputGain are the masters
    // outside everything, with polarity folded into outputGain's sign so a flip ramps
    // through zero rather than stepping; calibrationAmount fades the Set Up signal in over
    // the programme rather than cutting to it.
    lime::ControlSmoother inTrim, outTrim, calibrationAmount;
    lime::ControlSmoother inputGain, outputGain;

    // The three factors that decide how much of the wet path is heard. They multiply, and
    // they are separate because they answer at different speeds and for different reasons:
    // bypass is a monitoring switch and has to be quick, the transition covers a change of
    // process or mode and is deliberately slow, and mix is a control the user sets.
    lime::ControlSmoother bypassAmount, mixAmount;
    lime::Crossfade transition;

    // The dry path exists only for the crossfade, delayed to sit exactly under the
    // process's output. It also makes a full bypass exact rather than merely equivalent:
    // with the blend at zero the output is the input, delayed, to the last bit.
    std::vector<lime::DelayLine64> dryDelays;
    juce::AudioBuffer<double> dryBuffer;

    // What the engine is currently set to, as against what the parameters ask for. They
    // differ only while a change is fading through.
    int engagedProcessIndex = -1, engagedModeIndex = -1;
    bool configurationChanging = false;

    std::atomic<int> engagedProcess { -1 };
    std::atomic<int> engagedMode { -1 };

    std::atomic<float> meterDb { -100.0f };

    /** True while the host is driving processBlockBypassed rather than the bypass
        parameter. Folded into the same crossfade; atomic only so the editor can
        grey the plot from the message thread. */
    std::atomic<bool> hostBypassRequested { false };

    /** Set by setStateInformation, consumed by the next audio block: a restored
        state is not a control movement, so the switches are thrown and the
        smoothers snapped rather than faded over 1.5 s of the old settings. */
    std::atomic<bool> stateJustRestored { false };

    /** The block size prepareToPlay was promised, which is the most any chunk
        handed downstream may carry. */
    int preparedMaxBlock = 0;

    // Scratch used only by the float overload, allocated in prepareToPlay so the
    // audio thread never has to.
    juce::AudioBuffer<double> conversionBuffer;

    std::vector<double*> channelPointers;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LimeAudioProcessor)
};
