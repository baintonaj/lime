/*
  ==============================================================================

    MeterBar — the vertical level meter down the right edge.

    Its own component with its own timer, deliberately. Painting the meter inside the
    editor's paint() means the editor has to repaint to animate it, which redraws every
    knob, caption and button underneath thirty times a second and reads as flicker.
    Confining the animation to the one component that actually changes leaves the rest
    of the panel static.

    Scaled and captioned, because a bare coloured column tells the reader how much of
    something without saying what or how much.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "LimeStyle.h"
#include "LimeSurfaces.h"

//==============================================================================
class MeterBar  : public juce::Component,
                  private juce::Timer
{
public:
    /** @param source  returns the current level in dB; polled, never blocking. */
    explicit MeterBar (std::function<float()> source)
        : levelSource (std::move (source))
    {
        setOpaque (true);
        startTimerHz (30);
    }

    void paint (juce::Graphics& g) override
    {
        using namespace lime::style;

        // Everything that never moves — the caption, the well, the scale — is
        // painted once into a cached image at the display's physical resolution.
        // Repainting it live thirty times a second to animate a one-column meter
        // was most of this component's cost.
        const float pixelScale = juce::jmax (1.0f,
                                     g.getInternalContext().getPhysicalPixelScaleFactor());

        if (! backdrop.isValid() || std::abs (pixelScale - backdropScale) > 0.01f)
        {
            // If this paint is a partial one — the timer only invalidates the well —
            // the freshly re-rendered caption and scale would be clipped away and the
            // old raster would linger. Queue a full pass; this only happens when the
            // component first paints or moves between displays of different density.
            renderBackdrop (pixelScale);
            repaint();
        }

        g.drawImage (backdrop, getLocalBounds().toFloat());

        // The column itself, lit from the left like everything else on the panel.
        auto bar = well;
        bar.removeFromTop (bar.getHeight() * -displayedDb / lowerLimitDb);

        const auto fill = clipped ? juce::Colours::red : accent();

        g.setGradientFill (juce::ColourGradient (fill.brighter (0.35f), bar.getX(), bar.getCentreY(),
                                                 fill.darker (0.30f), bar.getRight(), bar.getCentreY(),
                                                 false));
        g.fillRect (bar);

        // A bright lip on top of the column so the level reads at a glance.
        if (bar.getHeight() > 2.0f)
        {
            g.setColour (fill.brighter (0.8f));
            g.fillRect (bar.withHeight (1.5f));
        }
    }

    void resized() override
    {
        backdrop = juce::Image();   // sizes changed; rebuild on the next paint
    }

private:
    void renderBackdrop (float pixelScale)
    {
        using namespace lime::style;

        backdropScale = pixelScale;
        backdrop = juce::Image (juce::Image::ARGB,
                                juce::jmax (1, juce::roundToInt ((float) getWidth() * pixelScale)),
                                juce::jmax (1, juce::roundToInt ((float) getHeight() * pixelScale)),
                                true);

        juce::Graphics g (backdrop);
        g.addTransform (juce::AffineTransform::scale (pixelScale));

        g.fillAll (background());

        // Inset so the caption sits on the same line as the top of the plot rather than
        // against the window frame.
        auto area = getLocalBounds().reduced (0, 20);

        g.setColour (muted());
        g.setFont (juce::FontOptions (fontSize * 0.6f));
        g.drawText ("out", area.removeFromTop (18), juce::Justification::centred);

        area.reduce (6, 4);
        const auto scale = area.removeFromLeft (30);

        well = lime::surfaces::paintWell (g, area.toFloat(),
                                          plotBackground().darker (0.35f), 3.0f);

        // The scale, so the column means a number of decibels rather than an amount.
        g.setFont (juce::FontOptions (fontSize * 0.52f));

        for (int db = 0; db >= (int) -lowerLimitDb; db -= 10)
        {
            const float y = well.getY() + well.getHeight() * (float) -db / lowerLimitDb;
            const bool labelled = db % 20 == 0;

            g.setColour (muted().withAlpha (labelled ? 0.75f : 0.35f));
            g.fillRect ((float) scale.getRight() - (labelled ? 6.0f : 3.0f), y - 0.5f,
                        labelled ? 6.0f : 3.0f, 1.0f);

            if (labelled)
            {
                g.setColour (muted());
                g.drawText (juce::String (db),
                            juce::Rectangle<float> ((float) scale.getX(), y - 7.0f,
                                                    (float) scale.getWidth() - 8.0f, 14.0f),
                            juce::Justification::centredRight);
            }
        }
    }

    void timerCallback() override
    {
        const float db = levelSource != nullptr ? levelSource() : -100.0f;
        const bool nowClipped = db >= 0.0f;
        const float limited = juce::jlimit (-lowerLimitDb, 0.0f, db);

        // Only repaint when the bar would actually move, so a quiet passage costs
        // nothing — and then only the well, which is the one region that animates.
        if (std::abs (limited - displayedDb) < 0.05f && nowClipped == clipped)
            return;

        displayedDb = limited;
        clipped = nowClipped;

        if (well.isEmpty())
            repaint();
        else
            repaint (well.getSmallestIntegerContainer().expanded (1));
    }

    static constexpr float lowerLimitDb = 60.0f;

    std::function<float()> levelSource;
    float displayedDb = -lowerLimitDb;
    bool clipped = false;

    juce::Image backdrop;
    float backdropScale = 0.0f;
    juce::Rectangle<float> well;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MeterBar)
};
