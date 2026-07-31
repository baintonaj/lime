/*
  ==============================================================================

    Lime — editor, in the Aptitude Audio house style established by Apti-Q.

    Same knobs, same toggles, same branding at the bottom left and level readout at the
    bottom right, same vertical meter down the right edge, laid out on the same 140 px
    grid — but in lime rather than periwinkle, since that is the plugin's name.

    Where Lime differs from Apti-Q it is because the process differs. Apti-Q is sixteen
    knobs; Lime is eight — two aligned rows of four, each filling the grid — and the
    freed space goes to the plot of the applied gain curve over the input spectrum,
    which is the most informative thing the plugin can show: boost rides above the zero
    line wherever the signal is quiet, cut mirrors it below in down mode.

    The top row is the four pass trims, named in the modes' own vocabulary — the in-trims
    set how hard the programme drives each pass, which is this architecture's threshold
    control, and the out-trims are each pass's level compensation. The bottom row is the
    plugin's own: the side-chain high and low passes, then the input and output
    gains. Every switch lives in one column to the right of both rows, so the eye finds
    them in one place rather than three.

    The selectors are rows of toggles rather than combo boxes, because Off/A/B is a
    three-position switch on the hardware this follows, not a menu, and because that is
    the vocabulary the range already uses. They carry words rather than abbreviations —
    the earlier "chk" and "byp" saved space the panel did not need and cost the reader the
    meaning.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "LimeLookAndFeel.h"
#include "MeterBar.h"
#include "PluginProcessor.h"
#include "SpectrumDisplay.h"

//==============================================================================
class LimeAudioProcessorEditor  : public juce::AudioProcessorEditor,
                                  private juce::Timer
{
public:
    explicit LimeAudioProcessorEditor (LimeAudioProcessor&);
    ~LimeAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    void timerCallback() override;

    /** A knob with its value box and caption, as the range lays them out. */
    struct Trim
    {
        juce::Slider slider;
        juce::Label caption;
        std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> attachment;
    };

    /** A row of square toggles standing in for one choice parameter, which is how a
        three-position switch reads on hardware. */
    struct Selector
    {
        std::vector<std::unique_ptr<juce::TextButton>> buttons;
        juce::Label caption;
        juce::RangedAudioParameter* parameter = nullptr;
    };

    /** @param textBoxWidthOverride  wider than the shared width where the value needs it;
                                     "100.00 %" does not fit what "+3.5 dB" was sized for. */
    void setUpTrim (Trim&, const juce::String& parameterId, const juce::String& caption,
                    int decimalPlaces = 1, int textBoxWidthOverride = 0);
    void setUpSelector (Selector&, const juce::String& parameterId,
                        const juce::String& caption);
    void setUpToggle (juce::TextButton&, const juce::String& parameterId,
                      const juce::String& caption,
                      std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment>&);

    void refreshSelectors();

    static juce::String formatLevel (float db);

    /** How wide a switch has to be to carry its word, measured in the same font the
        look-and-feel draws it with. The switches used to be fixed squares, which forced
        three-letter abbreviations — "chk", "byp" — that nobody can read cold. */
    static int widthFor (const juce::String& text);

    /** How wide a row of a selector's switches actually needs to be. */
    static int selectorWidth (const Selector&);

    /** How much room the whole switch block on the right needs — the widest of the three
        rows in it. Measured rather than reserved, so the panel is exactly as wide as its
        contents and the plot gets everything left over. */
    int buttonColumnWidth() const;

    LimeAudioProcessor& audioProcessor;
    LimeLookAndFeel lookAndFeel;

    SpectrumDisplay display;
    MeterBar meter;

    // Two aligned rows of four. The top row is the pass trims; the bottom is the
    // plugin's own — the side-chain high and low passes, then the masters that
    // exist because a plugin, unlike a module in a console, does not sit at a known
    // level. Declared in panel order, which is not signal order.
    Trim recIn, recOut, playIn, playOut;
    Trim scHpf, scLpf, inputGain, outputGain;
    Selector processSelector, modeSelector;

    juce::TextButton polarityButton, bypassButton, autoGainButton;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment>
        polarityAttachment, bypassAttachment, autoGainAttachment;

    // The side-chain filters' switches, in their own captioned row of the switch
    // column: the knobs hold the corners, these put them in the path.
    juce::TextButton scHpfButton, scLpfButton;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment>
        scHpfOnAttachment, scLpfOnAttachment;
    juce::Label scFilterCaption;

    // The panel carries no explanatory text: no running description of what the process is
    // doing, and no tooltips. Both were tried. The description took a whole row and changed
    // under the reader; tooltips only appear for someone already hovering the thing they
    // wanted explained. What a control does belongs in MANUAL.md, which can be read once
    // instead of rediscovered on every hover.

    juce::Label brandingLabel, levelLabel, autoCompareLabel;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LimeAudioProcessorEditor)
};
