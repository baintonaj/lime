/*
  ==============================================================================

    Lime — editor.

  ==============================================================================
*/

#include "PluginEditor.h"

#include "LimeSurfaces.h"

using namespace lime::style;

namespace
{
    /** Velocity sensitivity used across the range. */
    const double velocitySensitivity = std::sqrt (2.0);

    /** Caption size: the range's base type scaled by the same log10 idiom. */
    const float captionFontSize = fontSize
        * (float) std::pow (2.0, std::pow (10.0, -15.0 * std::log10 (std::sqrt (2.0)) / 3.0));

    /** Switch legends are set in capitals.

        Not the knob captions, which stay lower case: those name a control you are looking
        at, while a switch names a state the panel is asserting, and small capitals are how
        the rest of the range has always drawn a state. It also fixes the ragged x-height a
        row of equal slabs made obvious — "Off" beside "Both" set the two words at visibly
        different weights when only one of them had an ascender. */
    juce::String switchLegend (const juce::String& text)
    {
        return text.toUpperCase();
    }

    // One set of spacings for the whole panel, so nothing is aligned by coincidence.
    constexpr int margin = 23;
    constexpr int gap = 8;
    constexpr int captionHeight = 26;
    constexpr int footerHeight = 46;
    constexpr int buttonGap = 21;
}

//==============================================================================
LimeAudioProcessorEditor::LimeAudioProcessorEditor (LimeAudioProcessor& p)
    : AudioProcessorEditor (&p),
      audioProcessor (p),
      meter ([proc = &p] { return proc->getMeterDb(); })
{
    setLookAndFeel (&lookAndFeel);

    addAndMakeVisible (meter);
    addAndMakeVisible (display);

    // The engine only publishes curve snapshots while something is looking at them; see
    // SrEngine::setProbesActive. Switched off again in the destructor.
    //
    // The probes themselves are fetched on every timer tick rather than once here. They
    // belong to the engine's per-channel side chains, and SrEngine::prepare rebuilds those
    // from scratch — so a pointer taken at construction is null if the host builds the
    // editor before it prepares the processor, and dangling if it prepares again afterwards.
    // Either way the plot read "no signal" against a perfectly live engine.
    p.getEngine().setProbesActive (true);

    // The pass trims, named in the modes' own vocabulary: the in-trims set how hard
    // the programme drives each pass — the threshold control of this architecture —
    // and the out-trims are that pass's level compensation.
    setUpTrim (recIn, "recIn", "up in");
    setUpTrim (recOut, "recOut", "up out");
    setUpTrim (playIn, "playIn", "down in");
    setUpTrim (playOut, "playOut", "down out");

    // The bottom row, in the order the signal meets it: in, then the blend against
    // dry, then out. Mix carries two decimals and needs a wider box for them —
    // "100.00 %" is three characters longer than anything the trims print.
    setUpTrim (inputGain, "inputGain", "input");
    setUpTrim (mix, "mix", "mix", 2, 96);
    setUpTrim (outputGain, "outputGain", "output");

    setUpSelector (processSelector, "process", "process");
    setUpSelector (modeSelector, "mode", "mode");

    // Set Up is deliberately not on the panel. It sends a calibration signal so an unknown
    // channel can be aligned, which is a hardware-insert job: across a bit-transparent path
    // between two instances there is nothing to align, so the switch would be one a DAW user
    // never reaches for. The parameter and the whole calibration path remain — reachable from
    // the host's parameter list, and still measured by CalibrationTests — so the capability is
    // intact without spending panel space on it.
    //
    // Check Source is gone rather than hidden. It suppressed the decode so the return could
    // be heard as recorded, which is a real thing to want down a hardware insert — but in
    // down mode, which is where a DAW user would reach for it, suppressing the only pass in
    // circuit leaves delay and nothing else. It read as a second bypass switch, and a panel
    // with two switches that do the same thing in the common case is worse than a panel
    // without one of them.
    setUpToggle (polarityButton, "polarity", switchLegend ("polarity"), polarityAttachment);
    setUpToggle (bypassButton, "bypass", switchLegend ("bypass"), bypassAttachment);

    // Auto makeup: a dynamics unit without it makes every comparison a loudness
    // test. See AutoGain for what it matches and how slowly.
    setUpToggle (autoGainButton, "autoGain", switchLegend ("auto"), autoGainAttachment);

    addAndMakeVisible (brandingLabel);
    brandingLabel.setText ("Aptitude Audio | Lime | v" + juce::String (JucePlugin_VersionString),
                           juce::dontSendNotification);
    brandingLabel.setFont (juce::Font (juce::FontOptions (fontSize * (float) std::sqrt (2.0),
                                                          juce::Font::plain)));
    brandingLabel.setColour (juce::Label::textColourId, branding());
    brandingLabel.setJustificationType (juce::Justification::bottomLeft);
    brandingLabel.setBufferedToImage (true);

    addAndMakeVisible (levelLabel);
    levelLabel.setFont (juce::Font (juce::FontOptions (fontSize, juce::Font::plain)));
    levelLabel.setColour (juce::Label::textColourId, text());
    levelLabel.setJustificationType (juce::Justification::right);

    addAndMakeVisible (autoCompareLabel);
    autoCompareLabel.setFont (juce::Font (juce::FontOptions (fontSize * 0.7f, juce::Font::plain)));
    autoCompareLabel.setJustificationType (juce::Justification::centredLeft);

    // Fixed size, and honestly so.
    //
    // The panel was resizable, but resized() lays out in pixels — it removes slices of a
    // fixed height and steps the knobs along a fixed grid pitch — so dragging the corner
    // never scaled anything. It cropped: made larger the panel grew empty margin, and made
    // smaller the knobs silently dropped off the right-hand end one at a time, because each
    // one is placed only `if (row.getWidth() >= gridPitch)`. A handle that quietly deletes
    // controls is worse than no handle. Scaling properly would mean every metric in
    // LimeStyle becoming a ratio and the artwork re-rasterising per size, which is real work
    // for a panel whose one job is to be read; until that is done the window is the size it
    // was designed at.
    setResizable (false, false);

    // Four knobs abreast, twice, with every switch in one column beside them and the meter
    // strip down the edge. The width is measured rather than reserved — the switch column is
    // exactly as wide as the widest row of switches in it — so the plot gets whatever is
    // genuinely left over instead of a round number somebody guessed.
    setSize (4 * gridPitch + gap + buttonColumnWidth() + meterStrip + 2 * margin,
             547 + knobSize + textBoxHeight + captionHeight + gap);

    // Slow, and it only sets text and toggle states. Nothing here repaints the panel:
    // the meter and the plot animate themselves, and everything else is static.
    startTimerHz (8);
    timerCallback();

    // Development aid, off unless asked for: setting LIME_SNAPSHOT to a path makes the
    // panel photograph itself once, shortly after it opens. Worth keeping, because a
    // panel that cannot be looked at cannot be checked — this renders through exactly the
    // same paint calls the screen uses, so what it shows is what is on screen.
    //
    // Debug builds only. A release binary that deletes and writes a path taken from the
    // environment is a liability in a signed plugin, and the feature is for development.
   #if JUCE_DEBUG || defined (LIME_ENABLE_SNAPSHOT)
    if (const auto* path = std::getenv ("LIME_SNAPSHOT"))
    {
        juce::Timer::callAfterDelay (1500,
            [safe = juce::Component::SafePointer<LimeAudioProcessorEditor> (this),
             file = juce::File (juce::String (path))]
            {
                if (safe == nullptr)
                    return;

                const auto shot = safe->createComponentSnapshot (safe->getLocalBounds(), false);

                file.deleteFile();

                if (auto stream = file.createOutputStream())
                    juce::PNGImageFormat().writeImageToStream (shot, *stream);
            });
    }
   #endif
}

LimeAudioProcessorEditor::~LimeAudioProcessorEditor()
{
    audioProcessor.getEngine().setProbesActive (false);
    setLookAndFeel (nullptr);
}

//==============================================================================
void LimeAudioProcessorEditor::setUpTrim (Trim& trim, const juce::String& parameterId,
                                          const juce::String& caption,
                                          int decimalPlaces, int textBoxWidthOverride)
{
    const int boxWidth = textBoxWidthOverride > 0 ? textBoxWidthOverride : textBoxWidth;

    trim.slider.setSliderStyle (juce::Slider::RotaryHorizontalVerticalDrag);
    trim.slider.setTextBoxStyle (juce::Slider::TextBoxBelow, false, boxWidth, textBoxHeight);
    trim.slider.setVelocityModeParameters (velocitySensitivity, 1, 0.0, true,
                                           juce::ModifierKeys::ctrlAltCommandModifiers);
    trim.slider.setScrollWheelEnabled (true);
    trim.slider.setNumDecimalPlacesToDisplay (decimalPlaces);
    trim.slider.setPopupMenuEnabled (false);
    trim.slider.setColour (juce::Slider::textBoxHighlightColourId,
                           base().brighter().withSaturation (satMod()));
    addAndMakeVisible (trim.slider);

    trim.caption.setText (caption, juce::dontSendNotification);
    trim.caption.setFont (juce::Font (juce::FontOptions (captionFontSize, juce::Font::plain)));
    trim.caption.setColour (juce::Label::textColourId, text());
    trim.caption.setJustificationType (juce::Justification::centred);
    trim.caption.setBufferedToImage (true);
    addAndMakeVisible (trim.caption);

    trim.attachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment> (
        audioProcessor.getState(), parameterId, trim.slider);
}

void LimeAudioProcessorEditor::setUpSelector (Selector& selector, const juce::String& parameterId,
                                              const juce::String& caption)
{
    selector.parameter = audioProcessor.getState().getParameter (parameterId);

    juce::StringArray choices;

    if (auto* choice = dynamic_cast<juce::AudioParameterChoice*> (selector.parameter))
        choices = choice->choices;

    for (int i = 0; i < choices.size(); ++i)
    {
        auto button = std::make_unique<juce::TextButton>();

        // The parameter's own words, in capitals; see switchLegend.
        button->setButtonText (switchLegend (choices[i]));
        button->setColour (juce::TextButton::buttonColourId, buttonOff());
        button->setColour (juce::TextButton::buttonOnColourId, buttonOn());
        button->setColour (juce::TextButton::textColourOffId, text());

        // Dark text on the engaged face: the "on" colour is the brightest lime on the
        // panel, and near-white on it is barely legible.
        button->setColour (juce::TextButton::textColourOnId, background().darker (0.4f));
        button->setClickingTogglesState (false);
        button->setWantsKeyboardFocus (false);

        // Radio behaviour driven straight at the parameter, so host automation and the
        // panel never disagree about which position the switch is in. Wrapped in a
        // change gesture: hosts recording touch or latch automation key on the
        // begin/end pair, and a bare setValueNotifyingHost is invisible to some of
        // them — the attachment classes do exactly this for the knobs.
        auto* parameter = selector.parameter;
        const int index = i;
        const int count = choices.size();

        button->onClick = [parameter, index, count]
        {
            if (parameter != nullptr && count > 1)
            {
                parameter->beginChangeGesture();
                parameter->setValueNotifyingHost ((float) index / (float) (count - 1));
                parameter->endChangeGesture();
            }
        };

        addAndMakeVisible (*button);
        selector.buttons.push_back (std::move (button));
    }

    selector.caption.setText (caption, juce::dontSendNotification);
    selector.caption.setFont (juce::Font (juce::FontOptions (captionFontSize, juce::Font::plain)));
    selector.caption.setColour (juce::Label::textColourId, text());

    // Left aligned, so the caption starts over its own first switch. Centring it in the
    // cell floated it away from the row it names, since the switches are only as wide as
    // their words.
    selector.caption.setJustificationType (juce::Justification::centredLeft);
    selector.caption.setBufferedToImage (true);
    addAndMakeVisible (selector.caption);
}

void LimeAudioProcessorEditor::setUpToggle (
    juce::TextButton& button, const juce::String& parameterId, const juce::String& caption,
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment>& attachment)
{
    button.setButtonText (caption);
    button.setColour (juce::TextButton::buttonColourId, buttonOff());
    button.setColour (juce::TextButton::buttonOnColourId, buttonOn());
    button.setColour (juce::TextButton::textColourOffId, text());
    button.setColour (juce::TextButton::textColourOnId, background().darker (0.4f));
    button.setClickingTogglesState (true);
    button.setWantsKeyboardFocus (false);
    addAndMakeVisible (button);

    attachment = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment> (
        audioProcessor.getState(), parameterId, button);
}

//==============================================================================
juce::String LimeAudioProcessorEditor::formatLevel (float db)
{
    // The range shows one decimal down to -10 and whole numbers below, so the readout
    // never changes width mid-performance. Carries its unit, because a bare number beside
    // a meter could as easily be a percentage.
    if (db > -10.0f)
        return juce::String (std::round (db * 10.0f) / 10.0f, 1) + " dB";

    return juce::String ((int) std::round (db)) + " dB";
}

int LimeAudioProcessorEditor::widthFor (const juce::String& text)
{
    return juce::jmax (buttonSize,
                       juce::GlyphArrangement::getStringWidthInt (buttonFont(), text) + 30);
}

int LimeAudioProcessorEditor::selectorWidth (const Selector& selector)
{
    int width = 0;

    for (const auto& button : selector.buttons)
        width += widthFor (button->getButtonText()) + buttonGap;

    return juce::jmax (0, width - buttonGap);
}

int LimeAudioProcessorEditor::buttonColumnWidth() const
{
    const int toggles = widthFor (polarityButton.getButtonText()) + buttonGap
                      + widthFor (bypassButton.getButtonText()) + buttonGap
                      + widthFor (autoGainButton.getButtonText());

    // Two of the knob grid's own columns, so the switch block lines up with the rows beside
    // it rather than ending wherever its longest word happened to. The measured widths are
    // a floor rather than the answer: they only matter if a word ever outgrows the grid.
    return juce::jmax (2 * gridPitch,
                       selectorWidth (processSelector),
                       selectorWidth (modeSelector),
                       toggles);
}

//==============================================================================
void LimeAudioProcessorEditor::paint (juce::Graphics& g)
{
    // Static. The meter and the plot are components and animate themselves.
    lime::surfaces::paintPanel (g, getLocalBounds());
}

void LimeAudioProcessorEditor::resized()
{
    // Laid out entirely by removing slices, so no two things can overlap and the
    // vertical baselines are shared by construction rather than by arithmetic. Reading
    // order top to bottom: what the process is doing to the signal, said as a curve;
    // then said in words; then the switches that decide it; then the trims; then
    // identity and level.
    auto bounds = getLocalBounds();

    meter.setBounds (bounds.removeFromRight (meterStrip));

    auto inner = bounds.reduced (margin, margin);

    auto footer = inner.removeFromBottom (footerHeight);
    brandingLabel.setBounds (footer.removeFromLeft (juce::jmin (460, footer.getWidth() / 2)));
    levelLabel.setBounds (footer.removeFromRight (110));

    // Whatever is left in the middle. Empty unless Set Up has been engaged from the host, so
    // it costs nothing to leave room for it.
    autoCompareLabel.setBounds (footer);

    inner.removeFromBottom (gap);

    const int rowHeight = knobSize + textBoxHeight + captionHeight;

    // Two rows, bottom up. The masters go under the calibration trims rather than beside
    // them: eight knobs abreast would have made the panel half as wide again as the plot
    // needs, and the two groups are not the same kind of control. The top row is what is set
    // once against a channel — the trims, and how finely the band split is resolved. The
    // bottom row is the plugin's own, reached for while listening.
    auto masterRow = inner.removeFromBottom (rowHeight);
    inner.removeFromBottom (gap);
    auto trimRow = inner.removeFromBottom (rowHeight);

    // Every switch on the panel lives in one column on the right, against the two knob rows.
    // Scattering them — a pair beside the selectors, another pair under them — meant the eye
    // had to find them twice, and left the two rows different widths.
    auto switchColumn = trimRow.getUnion (masterRow)
                               .removeFromRight (buttonColumnWidth() + gap)
                               .withTrimmedLeft (gap);

    trimRow = trimRow.withTrimmedRight (buttonColumnWidth() + gap);
    masterRow = masterRow.withTrimmedRight (buttonColumnWidth() + gap);

    const auto placeTrim = [&] (Trim& trim, juce::Rectangle<int> cell)
    {
        const int boxWidth = juce::jmax (textBoxWidth, trim.slider.getTextBoxWidth());

        trim.slider.setBounds (cell.removeFromTop (knobSize + textBoxHeight)
                                   .withSizeKeepingCentre (knobSize + boxWidth,
                                                           knobSize + textBoxHeight));
        trim.caption.setBounds (cell);
    };

    // Four over three, the three centred under the four — the knob block reads as
    // an inverted trapezoid, wide row of pass trims above, the masters tucked
    // beneath. The half-pitch inset is what centres a row of three on a grid of
    // four.
    for (auto* trim : { &recIn, &recOut, &playIn, &playOut })
        if (trimRow.getWidth() >= gridPitch)
            placeTrim (*trim, trimRow.removeFromLeft (gridPitch));

    masterRow.removeFromLeft (gridPitch / 2);

    for (auto* trim : { &inputGain, &mix, &outputGain })
        if (masterRow.getWidth() >= gridPitch)
            placeTrim (*trim, masterRow.removeFromLeft (gridPitch));

    // Justified, like a line of type: whatever a row holds, it spans the column exactly.
    // Every switch therefore shares one of three widths rather than being cut to its own
    // word, which is what makes the block read as a grid instead of as a ragged list — and
    // it lets the switches be half again as tall as they were, since the rows now divide the
    // full height of the two knob rows beside them rather than sitting at a fixed size with
    // the remainder left empty underneath.
    const auto placeRow = [] (juce::Rectangle<int> row,
                              const std::vector<juce::Component*>& items)
    {
        if (items.empty())
            return;

        const int count = (int) items.size();
        const int spanned = row.getWidth() - (count - 1) * buttonGap;
        const int left = row.getX();

        for (int i = 0; i < count; ++i)
        {
            // Edges computed from the row's own width rather than accumulated, so rounding
            // cannot leave the last switch a pixel short of the margin.
            const int x0 = left + (spanned * i) / count + i * buttonGap;
            const int x1 = left + (spanned * (i + 1)) / count + i * buttonGap;

            items[(size_t) i]->setBounds (x0, row.getY(), x1 - x0, row.getHeight());
        }
    };

    const auto placeSelector = [&] (Selector& selector, juce::Rectangle<int> cell)
    {
        selector.caption.setBounds (cell.removeFromTop (captionHeight));

        std::vector<juce::Component*> items;

        for (auto& button : selector.buttons)
            items.push_back (button.get());

        placeRow (cell, items);
    };

    // Three rows in the switch column: the two selectors, then the two plain toggles side by
    // side. Each row is a caption slot over a row of switches, so every switch on the panel
    // sits on one of three shared baselines.
    {
        auto column = switchColumn;

        // Three rows sharing the full height, each a caption slot over its switches. The
        // toggles keep an empty caption slot rather than growing into it, so all nine
        // switch positions sit on the same three baselines.
        const int rows = 3;
        const int available = column.getHeight() - (rows - 1) * gap;
        const int rowStep = available / rows;

        placeSelector (processSelector, column.removeFromTop (rowStep));
        column.removeFromTop (gap);
        placeSelector (modeSelector, column.removeFromTop (rowStep));
        column.removeFromTop (gap);

        auto toggles = column.removeFromTop (rowStep);
        toggles.removeFromTop (captionHeight);

        placeRow (toggles, { &polarityButton, &autoGainButton, &bypassButton });
    }

    inner.removeFromBottom (gap);
    display.setBounds (inner);
}

//==============================================================================
void LimeAudioProcessorEditor::refreshSelectors()
{
    const auto sync = [] (Selector& selector)
    {
        if (auto* choice = dynamic_cast<juce::AudioParameterChoice*> (selector.parameter))
        {
            const int current = choice->getIndex();

            for (size_t i = 0; i < selector.buttons.size(); ++i)
                selector.buttons[i]->setToggleState ((int) i == current,
                                                     juce::dontSendNotification);
        }
    };

    sync (processSelector);
    sync (modeSelector);
}

void LimeAudioProcessorEditor::timerCallback()
{
    refreshSelectors();

    levelLabel.setText (formatLevel (audioProcessor.getMeterDb()), juce::dontSendNotification);

    const auto sampleRate = audioProcessor.getSampleRate();

    // The plot's fixed main-path response is designed per sample rate; the call is
    // idempotent, so telling it every tick costs nothing and it can never be stale.
    // Re-fetched every tick, so the plot follows the engine across a re-prepare instead of
    // holding a pointer into side chains that no longer exist. Idempotent and cheap.
    display.setProbes (audioProcessor.getEngine().getHfProbe (0),
                       audioProcessor.getEngine().getLfProbe (0));

    display.setSampleRate (sampleRate);
    display.setSections (audioProcessor.getEngine().getSections());
    display.setBypassed (audioProcessor.isBypassed());

    // The engaged process and mode rather than the parameters, so the plot changes
    // when the audio does rather than 400 ms ahead of it; see
    // LimeAudioProcessor::getEngagedProcessIndex.
    display.setProcessIndex (audioProcessor.getEngagedProcessIndex());
    display.setModeIndex (audioProcessor.getEngagedModeIndex());

    // The footer used to carry "64-bit", the latency in milliseconds and the sample rate.
    // None of the three was something a user could act on: the precision never changes, the
    // rate is the session's and is already on screen in the host, and the latency is
    // reported to the host and compensated automatically. What is left is the plot, the
    // meter and the controls, which is what the panel is for.

    // Auto Compare only means anything while Set Up is engaged. Red for the internal
    // reference, green for the tape return, as the module's lamps are.
    if (! audioProcessor.isSetUpEngaged())
    {
        autoCompareLabel.setText ({}, juce::dontSendNotification);
    }
    else
    {
        const bool reference = audioProcessor.isAutoCompareShowingReference();

        autoCompareLabel.setText (reference ? "auto compare: reference" : "auto compare: tape",
                                  juce::dontSendNotification);
        autoCompareLabel.setColour (juce::Label::textColourId,
                                   reference ? juce::Colours::orangered : accent());
    }
}
