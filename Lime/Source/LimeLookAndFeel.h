/*
  ==============================================================================

    LimeLookAndFeel — the house style: knobs and switches, drawn from the SVG artwork in
    LimeGraphics.

    The knob is a round object lit from the upper left, with a separate pointer rotated
    into place. The switches are rounded slabs, raised when off and pressed in when on,
    sized to the word they carry rather than to a fixed square — three-letter
    abbreviations saved space and cost the reader the meaning.

    The modifier constants in LimeStyle are powers of the square root of two written
    through log10, which is how the rest of the range expresses them — they are ratios in
    the same sense a decibel is, not arbitrary numbers.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "LimeGraphics.h"
#include "LimeStyle.h"

//==============================================================================
class LimeLookAndFeel  : public juce::LookAndFeel_V4
{
public:
    LimeLookAndFeel()
    {
        setColour (juce::Slider::textBoxTextColourId, lime::style::text());
        setColour (juce::Slider::textBoxBackgroundColourId, juce::Colours::transparentBlack);
        setColour (juce::Slider::textBoxOutlineColourId, juce::Colours::transparentBlack);
        setColour (juce::Label::textColourId, lime::style::text());
        setColour (juce::ComboBox::backgroundColourId, lime::style::buttonOff());
        setColour (juce::ComboBox::textColourId, lime::style::text());
        setColour (juce::ComboBox::outlineColourId, juce::Colours::transparentBlack);
        setColour (juce::ComboBox::arrowColourId, lime::style::text());
        setColour (juce::PopupMenu::backgroundColourId, lime::style::buttonOff());
        setColour (juce::PopupMenu::textColourId, lime::style::text());
        setColour (juce::PopupMenu::highlightedBackgroundColourId, lime::style::buttonOn());
    }

    void drawRotarySlider (juce::Graphics& g, int x, int y, int width, int height,
                           float sliderPos, float rotaryStartAngle, float rotaryEndAngle,
                           juce::Slider&) override
    {
        const float diameter = (float) juce::jmin (width, height);
        const auto centre = juce::Point<float> ((float) x + (float) width * 0.5f,
                                                (float) y + (float) height * 0.5f);

        // Square, so the artwork's circles stay circular: the cast shadow lives inside the
        // artwork's own box, so the knob is drawn at full size rather than inset and the
        // shadow falls where the light says it should.
        const auto area = juce::Rectangle<float> (diameter, diameter).withCentre (centre);

        knobBody.draw (g, area);

        // The pointer is authored pointing up and rotated into place, so its own shading
        // stays consistent with the light source at every angle. Rendered at the
        // display's physical resolution and scaled back down in the transform, so it
        // stays crisp on high-density displays; see CachedSvg::imageFor.
        const float angle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
        const int side = (int) std::ceil (diameter);
        const float pixelScale = lime::graphics::CachedSvg::physicalScaleFor (g);

        if (const auto& pointer = knobPointer.imageFor (side, side, pixelScale); pointer.isValid())
        {
            const auto transform = juce::AffineTransform::scale (1.0f / pixelScale)
                                       .rotated (angle, diameter * 0.5f, diameter * 0.5f)
                                       .translated (area.getX(), area.getY());

            g.drawImageTransformed (pointer, transform);
        }
    }

    void drawButtonBackground (juce::Graphics& g, juce::Button& button,
                               const juce::Colour&, bool hovered, bool down) override
    {
        const bool on = button.getToggleState() || down;

        (on ? buttonOnFace : buttonOffFace).draw (g, button.getLocalBounds().toFloat());

        // A hover lift, so the panel answers the mouse without changing state.
        if (hovered && ! down)
        {
            g.setColour (juce::Colours::white.withAlpha (0.07f));
            g.fillRoundedRectangle (button.getLocalBounds().toFloat().reduced (1.5f), 6.0f);
        }
    }

    juce::Font getTextButtonFont (juce::TextButton&, int) override
    {
        // The same font the editor measures the buttons with, so a label can never be
        // wider than the slab it sits on.
        return lime::style::buttonFont();
    }

    juce::Font getLabelFont (juce::Label& label) override
    {
        return label.getFont();
    }

private:
    // Rasterised once per size and reused, since knobs and buttons redraw on every mouse
    // move while the artwork itself never changes. The button faces are authored per size
    // rather than stretched, so their corners stay circular at any width.
    lime::graphics::CachedSvg knobBody { lime::graphics::knobBodySvg() };
    lime::graphics::CachedSvg knobPointer { lime::graphics::knobPointerSvg() };

    lime::graphics::CachedSvg buttonOffFace {
        lime::graphics::CachedSvg::Factory ([] (int w, int h)
            { return lime::graphics::buttonSvg (false, w, h); }) };

    lime::graphics::CachedSvg buttonOnFace {
        lime::graphics::CachedSvg::Factory ([] (int w, int h)
            { return lime::graphics::buttonSvg (true, w, h); }) };
};
