/*
  ==============================================================================

    LimeGraphics — the panel's SVG artwork, authored here rather than imported.

    The markup is built at runtime from the palette in LimeStyle, so the colours have one
    source of truth and restyling the range does not mean hand-editing hex in a dozen path
    fills.

    Everything is gradients and paths, with no SVG filter elements. JUCE's SVG parser
    covers linear and radial gradients well but its filter support is partial, so a
    feGaussianBlur drop shadow would either be ignored or render inconsistently. Soft
    edges are therefore built from stacked gradient stops, which is what gives the
    shading its depth without depending on anything the parser might skip.

    One light source throughout: upper left, as a panel lit by a room usually is. Every
    highlight sits on the upper-left arc and every contact shadow on the lower-right,
    because depth reads as wrong the moment two objects disagree about where the light
    is.

    Large resizable areas — the panel itself, the plot well, the meter — are drawn with
    JUCE gradients rather than SVG. A stretched SVG rectangle gains nothing over a
    gradient and costs a rasterisation, so SVG is used where the shape is intricate: the
    knobs and the buttons.

    ---------------------------------------------------------------------------
    Rendering notes, all of them found by rendering the artwork and looking at it rather
    than by reading documentation.

    **JUCE 9 replaced JUCE's own SVG parser with lunasvg**, so two of the three notes below
    are now history. They are kept because they explain why the markup is written the way it
    is, and deleting them would invite someone to undo the shape of it.

    1.  **Never scale SVG with Drawable::drawWithin().** It fits `getDrawableBounds()` to
        the target, and for a parsed SVG that is the union of the *drawn shapes*, not the
        viewBox. Artwork that deliberately occupies a corner of its box therefore gets
        stretched across the whole target. The knob pointer — a 3.6 x 22 strip in a
        100 x 100 box — came out as a 101 x 101 white slab covering the knob entirely,
        which is why the knobs looked square. `CachedSvg` scales by the declared viewBox
        instead, so where a shape sits in its box is preserved.

        Measured under the old parser and not re-tested under lunasvg, because it does not
        need to be: scaling by the declared box is correct whatever the parser thinks the
        drawn bounds are, so `CachedSvg` is right either way.

    2.  **Radial gradient radii had to be percentages — under the old parser only.** For the
        default `gradientUnits="objectBoundingBox"`, JUCE scaled `cx`/`cy` by the bounding
        box but computed `r` with `getCoordLength (r, width)`, which only scales a value
        written as a percentage. `r="0.78"` was therefore taken as 0.78 *pixels* — a
        sub-pixel dot, leaving the whole shape flooded with the gradient's last stop.

        lunasvg is a real SVG implementation and does not have this bug. Percentages are
        still used for every gradient coordinate here, because they are equally correct and
        the markup should not be rewritten on the strength of a parser change alone.

    3.  **The skirt gradient changed appearance at JUCE 9, and the new render is the right
        one.** Measured by photographing the panel before and after and differencing it: the
        plot well, meter and panel background are bit-identical, the buttons differ by at
        most 2/255, and the knobs differ by up to 79/255 — all of it in the skirt annulus
        between the dome edge and the outer rim, and all of it on the upper-left side:

            skirt colour, by angle       JUCE 8          JUCE 9
            315 degrees (upper left)     99,114,81       129,148,107
            0 degrees (top)              49,68,27        109,137,77
            45 / 135 / 180 / 225         unchanged       unchanged

        The old parser was compressing the lit end of the `skirt` linear gradient into a
        small patch near the top left; lunasvg spreads it across the whole upper-left half,
        which is what the markup asks for and what the one-light-source rule above says it
        should do. So this is a bug being fixed rather than a regression, and the artwork is
        left alone. `docs/panel.png` is the JUCE 9 render.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include <map>

#include "LimeStyle.h"

namespace lime::graphics
{
    /** Hex without the alpha, which is what SVG fill attributes want. */
    inline juce::String hex (juce::Colour c)
    {
        return "#" + c.toDisplayString (false);
    }

    //==========================================================================
    /** The knob body: skirt, dome, contact shadow and specular highlight.

        Drawn in a 100 x 100 box so the caller can scale it to any size. The stack from
        back to front is a cast shadow, the skirt ring catching light on its upper left,
        the domed cap, a rim shadow where the cap meets the skirt, and a soft specular.
    */
    inline juce::String knobBodySvg()
    {
        const auto body = style::base();
        const auto lit = body.brighter (0.9f);
        const auto dim = body.darker (0.55f);
        const auto skirtLit = body.brighter (0.35f);
        const auto skirtDim = body.darker (0.7f);

        juce::String svg;
        svg << R"SVG(<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">)SVG"
            << "<defs>"

            // Skirt: lit from upper left, falling away to the lower right.
            << R"SVG(<linearGradient id="skirt" x1="15%" y1="5%" x2="85%" y2="95%">)SVG"
            << R"SVG(<stop offset="0" stop-color=")SVG" << hex (skirtLit) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="0.5" stop-color=")SVG" << hex (body) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="1" stop-color=")SVG" << hex (skirtDim) << R"SVG("/>)SVG"
            << "</linearGradient>"

            // Cap: a dome, so the gradient is radial and offset toward the light.
            << R"SVG(<radialGradient id="dome" cx="36%" cy="30%" r="78%">)SVG"
            << R"SVG(<stop offset="0" stop-color=")SVG" << hex (lit) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="0.45" stop-color=")SVG" << hex (body) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="1" stop-color=")SVG" << hex (dim) << R"SVG("/>)SVG"
            << "</radialGradient>"

            // Contact shadow under the cap's lower-right edge. Several stops rather than
            // a blur, since the parser may ignore filters.
            << R"SVG(<radialGradient id="rim" cx="50%" cy="50%" r="50%">)SVG"
            << R"SVG(<stop offset="0.74" stop-color="#000000" stop-opacity="0"/>)SVG"
            << R"SVG(<stop offset="0.92" stop-color="#000000" stop-opacity="0.16"/>)SVG"
            << R"SVG(<stop offset="1" stop-color="#000000" stop-opacity="0.42"/>)SVG"
            << "</radialGradient>"

            // Specular: strongest at the top left, gone by the middle.
            << R"SVG(<radialGradient id="gloss" cx="34%" cy="24%" r="70%">)SVG"
            << R"SVG(<stop offset="0" stop-color="#ffffff" stop-opacity="0.34"/>)SVG"
            << R"SVG(<stop offset="0.55" stop-color="#ffffff" stop-opacity="0.07"/>)SVG"
            << R"SVG(<stop offset="1" stop-color="#ffffff" stop-opacity="0"/>)SVG"
            << "</radialGradient>"

            // The shadow the whole knob casts on the panel, below and to the right.
            << R"SVG(<radialGradient id="cast" cx="50%" cy="50%" r="50%">)SVG"
            << R"SVG(<stop offset="0.62" stop-color="#000000" stop-opacity="0.30"/>)SVG"
            << R"SVG(<stop offset="0.86" stop-color="#000000" stop-opacity="0.10"/>)SVG"
            << R"SVG(<stop offset="1" stop-color="#000000" stop-opacity="0"/>)SVG"
            << "</radialGradient>"
            << "</defs>"

            << R"SVG(<circle cx="52" cy="53" r="48" fill="url(#cast)"/>)SVG"
            << R"SVG(<circle cx="50" cy="50" r="46" fill="url(#skirt)"/>)SVG"
            << R"SVG(<circle cx="50" cy="50" r="46" fill="none" stroke="#ffffff")SVG"
            << R"SVG( stroke-opacity="0.10" stroke-width="1"/>)SVG"
            << R"SVG(<circle cx="50" cy="50" r="37" fill="url(#dome)"/>)SVG"
            << R"SVG(<circle cx="50" cy="50" r="37" fill="url(#rim)"/>)SVG"
            << R"SVG(<ellipse cx="41" cy="34" rx="24" ry="19" fill="url(#gloss)"/>)SVG"
            << "</svg>";

        return svg;
    }

    //==========================================================================
    /** The pointer, in the same 100 x 100 box, pointing up from the centre.

        A dark copy sits a fraction below the light one so the pointer reads as sitting
        on the cap rather than being painted onto it.
    */
    inline juce::String knobPointerSvg()
    {
        juce::String svg;
        svg << R"SVG(<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">)SVG"
            << R"SVG(<rect x="47.6" y="15.5" width="4.8" height="26" rx="2.4")SVG"
            << R"SVG( fill="#000000" fill-opacity="0.32"/>)SVG"
            << R"SVG(<rect x="47.4" y="14.5" width="4.8" height="26" rx="2.4" fill=")SVG"
            << hex (style::text()) << R"SVG("/>)SVG"
            << "</svg>";

        return svg;
    }

    //==========================================================================
    /** A button face, authored at the size it will be drawn.

        The viewBox is the pixel box, so the corner radius stays circular whatever the
        button's aspect ratio — which matters because the labels are words of different
        lengths and the buttons are sized to their text.

        Raised when off — lit along the top edge, shaded along the bottom, casting a
        shadow. Inset when on, with the shading inverted and a shadow along the *top*
        inner edge, so an engaged switch reads as pressed into the panel rather than
        merely recoloured.
    */
    inline juce::String buttonSvg (bool on, int width, int height)
    {
        const auto face = on ? style::buttonOn() : style::buttonOff();
        const auto top = on ? face.darker (0.30f) : face.brighter (0.32f);
        const auto bottom = on ? face.brighter (0.16f) : face.darker (0.34f);

        const auto w = juce::String (juce::jmax (1, width));
        const auto h = juce::String (juce::jmax (1, height));

        // The face is inset by a pixel and a half all round, leaving room for the cast
        // shadow below it without the button growing.
        const auto faceW = juce::String (juce::jmax (1.0f, (float) width - 3.0f), 1);
        const auto faceH = juce::String (juce::jmax (1.0f, (float) height - 4.0f), 1);

        juce::String svg;
        svg << R"SVG(<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 )SVG"
            << w << " " << h << R"SVG(">)SVG"
            << "<defs>"
            << R"SVG(<linearGradient id="face" x1="0%" y1="0%" x2="0%" y2="100%">)SVG"
            << R"SVG(<stop offset="0" stop-color=")SVG" << hex (top) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="0.55" stop-color=")SVG" << hex (face) << R"SVG("/>)SVG"
            << R"SVG(<stop offset="1" stop-color=")SVG" << hex (bottom) << R"SVG("/>)SVG"
            << "</linearGradient>";

        if (on)
        {
            // Pressed: shadow cast inward from the upper edge.
            svg << R"SVG(<linearGradient id="inner" x1="0%" y1="0%" x2="0%" y2="100%">)SVG"
                << R"SVG(<stop offset="0" stop-color="#000000" stop-opacity="0.34"/>)SVG"
                << R"SVG(<stop offset="0.32" stop-color="#000000" stop-opacity="0.05"/>)SVG"
                << R"SVG(<stop offset="1" stop-color="#000000" stop-opacity="0"/>)SVG"
                << "</linearGradient>";
        }

        svg << "</defs>";

        const auto faceRect = [&] (const juce::String& fill, float dy)
        {
            return R"SVG(<rect x="1.5" y=")SVG" + juce::String (1.5f + dy, 1)
                 + R"SVG(" width=")SVG" + faceW + R"SVG(" height=")SVG" + faceH
                 + R"SVG(" rx="6" )SVG" + fill + "/>";
        };

        if (! on)
            svg << faceRect (R"SVG(fill="#000000" fill-opacity="0.22")SVG", 1.5f);

        svg << faceRect (R"SVG(fill="url(#face)")SVG", 0.0f);

        if (on)
            svg << faceRect (R"SVG(fill="url(#inner)")SVG", 0.0f);

        svg << faceRect (R"SVG(fill="none" stroke="#ffffff" stroke-opacity=")SVG"
                         + juce::String (on ? "0.10" : "0.22")
                         + R"SVG(" stroke-width="1")SVG", 0.0f)
            << "</svg>";

        return svg;
    }

    //==========================================================================
    /** Parses one of the strings above. Null only if the markup is malformed, which
        would be a programming error here rather than a runtime condition.

        This went through juce::parseXML and Drawable::createFromSVG until JUCE 9, which
        removed that overload: the SVG support is lunasvg now, and lunasvg does its own XML
        parsing rather than accepting a juce::XmlElement. Handing it the string directly is
        what the replacement wants, and it is one step shorter than what it replaced. */
    inline std::unique_ptr<juce::Drawable> parse (const juce::String& svg)
    {
        return juce::Drawable::createFromSVGString (svg);
    }

    /** The declared viewBox extent, which is the artwork's real box — see the note at the
        top of this file about why the drawable's own bounds must not be used. */
    inline juce::Point<float> viewBoxOf (const juce::String& svg)
    {
        const auto values = juce::StringArray::fromTokens (
            svg.fromFirstOccurrenceOf ("viewBox=\"", false, false)
               .upToFirstOccurrenceOf ("\"", false, false), " ,", "");

        if (values.size() == 4)
        {
            const auto w = values[2].getFloatValue();
            const auto h = values[3].getFloatValue();

            if (w > 0.0f && h > 0.0f)
                return { w, h };
        }

        return { 100.0f, 100.0f };
    }

    //==========================================================================
    /** Rasterises artwork once per size and reuses it.

        Knobs and buttons redraw on every mouse move and every parameter change, and
        re-rasterising vector artwork each time is wasteful when the artwork itself never
        changes. Sizes are kept in a map rather than a single slot because one look-and-feel
        serves buttons of several widths, and a single slot would have them evicting each
        other on every repaint.

        Constructed either with fixed markup, or with a factory that authors the markup for
        a given pixel size — which is how the buttons keep their corner radius circular at
        any width.
    */
    class CachedSvg
    {
    public:
        using Factory = std::function<juce::String (int, int)>;

        explicit CachedSvg (juce::String markup)
            : fixed (std::move (markup)) {}

        explicit CachedSvg (Factory f)
            : factory (std::move (f)) {}

        void draw (juce::Graphics& g, juce::Rectangle<float> area)
        {
            // Rasterised at the display's physical resolution, not the logical one:
            // a 101-point knob drawn from a 101-pixel bitmap on a 2x display is
            // upscaled by the compositor and reads soft. The markup is still
            // authored at the logical size — insets and corner radii are designed
            // in points — and only the raster is denser.
            const auto pixelScale = physicalScaleFor (g);
            const auto& rendered = imageFor ((int) std::ceil (area.getWidth()),
                                             (int) std::ceil (area.getHeight()), pixelScale);

            if (rendered.isValid())
                g.drawImage (rendered, area, juce::RectanglePlacement::stretchToFit);
        }

        /** The display scale a graphics context is heading for, quantised so cache
            keys stay discrete. */
        static float physicalScaleFor (juce::Graphics& g)
        {
            return (float) juce::roundToInt (4.0f * juce::jmax (1.0f,
                       g.getInternalContext().getPhysicalPixelScaleFactor())) * 0.25f;
        }

        /** For callers that need to draw it through a transform, such as a rotating
            pointer. `width` and `height` are logical; the returned image is
            `pixelScale` times denser, so such callers must fold 1/pixelScale into
            their transform. */
        const juce::Image& imageFor (int width, int height, float pixelScale = 1.0f)
        {
            width = juce::jmax (1, width);
            height = juce::jmax (1, height);
            pixelScale = juce::jmax (1.0f, pixelScale);

            auto& slot = images[{ width, height, juce::roundToInt (pixelScale * 100.0f) }];

            if (! slot.isValid())
                slot = render (width, height, pixelScale);

            return slot;
        }

    private:
        juce::Image render (int width, int height, float pixelScale) const
        {
            const int pixelWidth = juce::jmax (1, juce::roundToInt ((float) width * pixelScale));
            const int pixelHeight = juce::jmax (1, juce::roundToInt ((float) height * pixelScale));

            juce::Image image (juce::Image::ARGB, pixelWidth, pixelHeight, true);

            // The factory is asked for the *logical* size, so markup that speaks in
            // pixels — a 1.5 point inset, a 6 point corner — keeps its designed
            // proportions on every display; only the raster density changes.
            const auto markup = factory != nullptr ? factory (width, height) : fixed;

            if (auto drawable = parse (markup))
            {
                const auto box = viewBoxOf (markup);

                // Scaled by the declared viewBox, never by drawWithin: see the note at the
                // top of this file.
                juce::Graphics g (image);
                drawable->draw (g, 1.0f,
                                juce::AffineTransform::scale ((float) pixelWidth / box.x,
                                                              (float) pixelHeight / box.y));
            }

            return image;
        }

        juce::String fixed;
        Factory factory;
        std::map<std::tuple<int, int, int>, juce::Image> images;
    };
}
