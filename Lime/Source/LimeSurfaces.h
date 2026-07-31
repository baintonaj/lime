/*
  ==============================================================================

    LimeSurfaces — the large surfaces: the panel, and the wells cut into it.

    Drawn with JUCE gradients rather than SVG. A stretched SVG rectangle gains nothing
    over a gradient and costs a rasterisation on every resize, so the vector artwork is
    reserved for the intricate fixed-aspect pieces — knobs and buttons — and the flat
    expanses are shaded directly.

    Same light source as the artwork: upper left. A well is therefore darkest along its
    top and left inner edges, where the lip casts into it, and carries a thin highlight
    along the bottom and right, where the far wall catches the light. Getting that
    backwards is what makes a recess look like a raised slab.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "LimeStyle.h"

namespace lime::surfaces
{
    /** The panel: a shallow vertical gradient plus a vignette, so a large area does not
        read as one flat wash of colour. */
    inline void paintPanel (juce::Graphics& g, juce::Rectangle<int> bounds)
    {
        auto area = bounds.toFloat();
        const auto base = style::background();

        g.setGradientFill (juce::ColourGradient (base.brighter (0.10f), area.getCentreX(), area.getY(),
                                                 base.darker (0.13f), area.getCentreX(), area.getBottom(),
                                                 false));
        g.fillRect (area);

        // Corners fall away slightly, which reads as a panel rather than a rectangle.
        juce::ColourGradient vignette (juce::Colours::transparentBlack,
                                       area.getCentreX(), area.getCentreY(),
                                       juce::Colours::black.withAlpha (0.16f),
                                       area.getX(), area.getBottom(), true);
        g.setGradientFill (vignette);
        g.fillRect (area);

        // A lit top edge and a shaded bottom one give the panel a thickness.
        g.setColour (juce::Colours::white.withAlpha (0.16f));
        g.fillRect (area.withHeight (1.0f));
        g.setColour (juce::Colours::black.withAlpha (0.18f));
        g.fillRect (area.removeFromBottom (1.0f));
    }

    /** A recess cut into the panel, for the plot and the meter to sit in. Returns the
        interior, so callers draw their contents inside the lip rather than under it. */
    inline juce::Rectangle<float> paintWell (juce::Graphics& g, juce::Rectangle<float> area,
                                             juce::Colour interior, float corner = 4.0f)
    {
        // The lip: light below and right, dark above and left, which is the inverse of a
        // raised object and is what makes the eye read a hole.
        g.setColour (juce::Colours::white.withAlpha (0.14f));
        g.fillRoundedRectangle (area.translated (0.0f, 1.0f), corner);

        g.setGradientFill (juce::ColourGradient (interior.darker (0.22f), area.getCentreX(), area.getY(),
                                                 interior.brighter (0.05f), area.getCentreX(), area.getBottom(),
                                                 false));
        g.fillRoundedRectangle (area, corner);

        // Shadow cast inward from the top lip. Stacked strips rather than a blur, to stay
        // cheap on a component that repaints.
        for (int i = 0; i < 5; ++i)
        {
            const float alpha = 0.16f * (1.0f - (float) i / 5.0f);
            g.setColour (juce::Colours::black.withAlpha (alpha));
            g.drawRoundedRectangle (area.reduced ((float) i * 0.9f + 0.5f), corner, 1.0f);
        }

        g.setColour (juce::Colours::black.withAlpha (0.35f));
        g.drawRoundedRectangle (area, corner, 1.0f);

        return area.reduced (3.0f);
    }
}
