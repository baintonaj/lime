/*
  ==============================================================================

    LimeStyle — the palette and panel metrics.

    Its own header because both the artwork and the look-and-feel need it, and having
    either include the other would be a cycle.

    Every surface is a shade of the same thing: one hue family from 76 to 100 degrees,
    which is yellow-green through to green. The plugin is called Lime, so the panel is
    lime — but a panel painted in one flat bright green would be both garish and
    unreadable, so the hue is held roughly constant while saturation and brightness carry
    the difference. Dark surfaces, light text, and the brightest, most saturated lime
    reserved for the few things that need to be found at a glance: an engaged switch, the
    level column, the encode curve.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

namespace lime::style
{
    // Functions, not namespace-scope objects, and deliberately so.
    //
    // The obvious spelling is `inline const juce::Colour base = ...;` but that is a static
    // initialisation order bug waiting to happen when the value depends on a global in
    // another translation unit, and the order between them is undefined. When it loses the
    // race every colour here comes out garbage or transparent — which showed up as a
    // default-grey panel with the desktop visible through the plot, dark knobs and
    // invisible captions. Returning by value from a function has no initialisation order
    // at all.

    /** One place where a lime is mixed, so no shade in the panel can drift off the family.
        @param hueDegrees  76 to 100; below that it turns yellow, above it turns grass. */
    inline juce::Colour lime (float hueDegrees, float saturation, float brightness)
    {
        return juce::Colour::fromHSV (hueDegrees / 360.0f, saturation, brightness, 1.0f);
    }

    /** The knob body, and the seed for the artwork's shading. */
    inline juce::Colour base()          { return lime (88.0f, 0.60f, 0.46f); }

    /** Captions, values and the knob pointers: near-white, with just enough lime in it
        that it belongs to the panel rather than sitting on top of it. */
    inline juce::Colour text()          { return lime (80.0f, 0.05f, 0.97f); }

    /** Secondary text — units, status, section labels. Dim enough to recede, light enough
        to read, which is the whole job. */
    inline juce::Colour muted()         { return lime (84.0f, 0.12f, 0.72f); }

    /** The panel: deep and only lightly saturated, so everything on it can be brighter
        than it. */
    inline juce::Colour background()    { return lime (96.0f, 0.34f, 0.21f); }

    /** Branding sits in bright lime. It was black, which reads as engraved on a pale
        panel and as absent on a dark one. */
    inline juce::Colour branding()      { return lime (84.0f, 0.62f, 0.74f); }

    /** The one loud colour, for the few things that must be found immediately. */
    inline juce::Colour accent()        { return lime (78.0f, 0.85f, 0.86f); }

    inline juce::Colour buttonOff()     { return lime (94.0f, 0.22f, 0.34f); }
    inline juce::Colour buttonOn()      { return lime (82.0f, 0.80f, 0.62f); }

    /** The plot sits darker than the panel so the curve has something to read against,
        while staying obviously part of the same object. */
    inline juce::Colour plotBackground() { return lime (100.0f, 0.34f, 0.17f); }

    /** Ratios in the same sense a decibel is, which is how the rest of the range writes
        them: powers of the square root of two through log10.

            knobMod  ~ 1.0595   the rim's size ratio
            alphaMod = 0.5      the rim's transparency
            satMod   = 2.0      the rim's saturation
    */
    inline float knobMod()  { return std::pow (10.0f,  0.5f * std::log10 (std::sqrt (2.0f)) / 3.0f); }
    inline float alphaMod() { return std::pow (10.0f, -6.0f * std::log10 (std::sqrt (2.0f)) / 3.0f); }
    inline float satMod()   { return std::pow (10.0f,  6.0f * std::log10 (std::sqrt (2.0f)) / 3.0f); }

    /** Base type size for the range; captions and titles are ratios of it. */
    inline constexpr float fontSize = 20.0f;

    /** The switches carry words rather than three-letter abbreviations, so their font is
        defined once here and used both to draw them and to measure how wide they need to
        be. */
    inline juce::Font buttonFont()
    {
        // A root of two up from where it started, less 7%. The switches were sized to their
        // words when the words were the only thing setting their size; now that they fill a
        // grid the slabs are much larger than the type on them, and 0.72 of the base size
        // read as small lettering on a big button rather than as a label. The 7% back off
        // is what set in capitals took to stop crowding the slab it sits on — capitals have
        // no descenders, so they fill more of a line than the same size in lower case.
        return juce::Font (juce::FontOptions (fontSize * 0.72f * (float) std::sqrt (2.0) * 0.93f,
                                              juce::Font::plain));
    }

    /** Grid pitch between knob centres, and the text box geometry beneath them. */
    inline constexpr int gridPitch = 140;
    inline constexpr int textBoxWidth = 68;
    inline constexpr int textBoxHeight = 22;
    inline constexpr int knobSize = 101;
    inline constexpr int buttonSize = 45;
    inline constexpr int meterStrip = 74;
}
