/*
  ==============================================================================

    SpectrumDisplay — the input spectrum with the applied gain curve over it.

    Watching the gain curve move against the signal is what makes the least treatment
    principle legible. The scale is symmetric about an emphasised 0 dB centre line:
    in up and both modes the encode boost rides above it, sitting high wherever the
    signal is quiet and dipping only around whatever is loud; in down mode the
    complementary cut mirrors that below the line.

    It also has to be honest about when it is not looking at anything. The curve comes from
    the spectral process's own side chains, so with the process off or on A there is nothing
    publishing frames, and under bypass what it would show is not in circuit. In all three
    the trace says so rather than leaving the last frame up as though it were current.

    Annotated deliberately. An unlabelled pair of traces tells a reader nothing, so
    there is a labelled frequency axis, a boost scale, and a legend naming both traces.

    The two halves are analysed at different resolutions and are *combined* rather than
    spliced where they meet at 800 Hz — see gather(), which is where an apparent
    discontinuity in the process turned out to be one in the plot.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "LimeStyle.h"
#include "LimeSurfaces.h"
#include "dsp/Biquad64.h"
#include "dsp/FilterDesign.h"
#include "dsp/SpectrumProbe.h"
#include "dsp/SrSideChain.h"

//==============================================================================
class SpectrumDisplay  : public juce::Component,
                         private juce::Timer
{
public:
    SpectrumDisplay()
    {
        setOpaque (true);

        // 15 Hz is plenty for reading a curve, and the callback only repaints when the
        // engine has published a new frame, so a silent session costs nothing and the
        // panel stays still.
        startTimerHz (15);
    }

    /** Probes are owned by the engine and outlive the editor. Either may be null. */
    void setProbes (lime::SpectrumProbe* hfSource, lime::SpectrumProbe* lfSource)
    {
        hfProbe = hfSource;
        lfProbe = lfSource;
    }

    /** The fixed main-path networks are designed per sample rate, so the curve has to be
        told when it changes. Cheap and idempotent, so the editor can call it every tick. */
    void setSampleRate (double rate)
    {
        if (rate <= 0.0 || juce::exactlyEqual (rate, fixedResponseRate))
            return;

        fixedResponseRate = rate;
        fixedResponseDb.assign ((size_t) gridPoints, 0.0);

        const auto cascade = lime::srNetworks::mainPathBiquadsForRate (rate);

        for (int i = 0; i < gridPoints; ++i)
        {
            const double omega = twoPi * gridHz (i) / rate;
            std::complex<double> h { 1.0, 0.0 };

            for (const auto& section : cascade)
                h *= section.responseAt (omega);

            const auto magnitude = std::abs (h);
            fixedResponseDb[(size_t) i] = 20.0 * std::log10 (std::max (1.0e-12, magnitude));
        }

        repaint();
    }

    /** How finely the band-defining split between the two halves is resolved, which the
        panel draws as section boundaries. The editor hands this over on every tick, so a
        count that has not moved has to cost nothing rather than provoke a repaint. */
    void setSections (int count)
    {
        if (count <= 0 || count == sections)
            return;

        sections = count;
        repaint();
    }

    /** Bypassed: the plot holds its last frame and greys.

        A live trace under a bypassed plugin says the process is doing something it is not.
        Freezing rather than blanking keeps the last thing it *was* doing on screen to
        compare against, which is usually why bypass was pressed. */
    void setBypassed (bool shouldBeBypassed)
    {
        if (shouldBeBypassed == bypassed)
            return;

        bypassed = shouldBeBypassed;
        refreshAndRepaint();
    }

    /** Which process is engaged: 0 off, 1 A, 2 B. */
    void setProcessIndex (int index)
    {
        if (index == processIndex)
            return;

        processIndex = index;
        refreshAndRepaint();
    }

    /** Which mode is engaged: 0 up, 1 down, 2 both. Down shows the applied cut,
        so the fixed main path joins the curve inverted, exactly as it runs. */
    void setModeIndex (int index)
    {
        if (index == modeIndex)
            return;

        modeIndex = index;
        refreshAndRepaint();
    }

    void paint (juce::Graphics& g) override
    {
        using namespace lime::style;

        // The panel shows through the well's lip, so paint it first.
        g.fillAll (background());

        auto full = lime::surfaces::paintWell (g, getLocalBounds().toFloat().reduced (1.0f),
                                               plotBackground()).getSmallestIntegerContainer()
                        .reduced (3);

        // Room for the scales, so no trace ever runs underneath its own labels.
        const auto scaleColumn = full.removeFromRight (34);
        const auto axisRow = full.removeFromBottom (16);
        const auto plot = full.toFloat();

        drawGrid (g, plot, axisRow, scaleColumn);

        // Switched off, the process is applying nothing, and nothing is what the plot should
        // say. The probes are not running either — SrEngine::process pays the latency as
        // pure delay and never reaches a side chain — so the last input spectrum is stale
        // and is left out rather than drawn as though it were current. A flat line at 0 dB
        // is the whole truth here.
        if (processIndex == 0)
        {
            // Flat at 0 dB, which on the symmetric scale is the centre of the plot:
            // "the process is applying nothing" drawn as exactly that.
            const float y = yFor (plot, 0.0, minGainDb, maxGainDb);

            g.setColour (dim (text()).withMultipliedAlpha (0.85f));
            g.drawLine (plot.getX(), y, plot.getRight(), y, 2.0f);

            g.setColour (muted().withAlpha (0.85f));
            g.setFont (juce::FontOptions (fontSize * 0.7f));

            // Plain ASCII, and it has to be. juce::String's const char* constructor reads
            // its input as CharPointer_ASCII — one character per byte — so a UTF-8 em dash
            // arrives as the three bytes it is made of and draws as "a-with-circumflex, euro
            // sign, box". Anything non-ASCII on this panel has to go through
            // String::fromUTF8, as the degree sign on the phase readout does. Punctuation
            // that needs no escaping is the safer answer for a caption.
            g.drawText (bypassed ? "process off, bypassed" : "process off",
                        plot, juce::Justification::centred);

            drawLegend (g, plot);
            return;
        }

        const bool ready = (hfProbe != nullptr && hfProbe->isReady())
                        || (lfProbe != nullptr && lfProbe->isReady());

        if (! ready)
        {
            g.setColour (muted().withAlpha (0.7f));
            g.setFont (juce::FontOptions (fontSize * 0.7f));
            g.drawText ("no signal", plot, juce::Justification::centred);
            return;
        }

        // The traces are gathered on the timer rather than here: paint() can run
        // several times between frames — an exposed window, an overlapping repaint —
        // and re-walking 480 grid points of pow and log10 for identical data buys
        // nothing. Frozen states simply keep whatever the timer last produced.

        // Input spectrum: filled, behind, in the knob colour so it reads as signal.
        if (auto path = buildPath (plot, levelPoints, minLevelDb, maxLevelDb); ! path.isEmpty())
        {
            auto filled = path;
            filled.lineTo (plot.getRight(), plot.getBottom());
            filled.lineTo (plot.getX(), plot.getBottom());
            filled.closeSubPath();

            g.setColour (dim (base().brighter (0.2f)).withMultipliedAlpha (0.35f));
            g.fillPath (filled);

            g.setColour (dim (base().brighter (0.5f)));
            g.strokePath (path, juce::PathStrokeType (1.0f));
        }

        // The boost curve, in the same whitesmoke as the knob pointers, so the eye reads
        // it as the control surface rather than as signal.
        if (auto path = buildPath (plot, gainPoints, minGainDb, maxGainDb); ! path.isEmpty())
        {
            g.setColour (dim (text()));
            g.strokePath (path, juce::PathStrokeType (2.0f));
        }

        drawLegend (g, plot);

        if (isFrozen())
        {
            g.setColour (muted().withAlpha (0.85f));
            g.setFont (juce::FontOptions (fontSize * 0.7f));
            g.drawText (bypassed ? "bypassed" : "process a: no spectral analysis",
                        plot, juce::Justification::centred);
        }
    }

private:
    /** True where the trace is holding a frame rather than following the signal.

        Two cases, and they look the same on purpose because they mean the same thing to a
        reader: what you are seeing is not what is happening now. Bypass, where the process
        is out of circuit; and process A, whose four fixed bands are a time-domain network
        with no spectral analysis behind it, so there is nothing publishing frames. */
    bool isFrozen() const noexcept { return bypassed || processIndex == 1; }

    /** A state change that un-freezes the trace gathers before it paints, so the
        first frame after releasing bypass or re-engaging the process is current
        rather than whatever the timer last produced before the freeze — which may
        be minutes old, and this panel's one rule is never to present a held frame
        as a live one. */
    void refreshAndRepaint()
    {
        const bool ready = (hfProbe != nullptr && hfProbe->isReady())
                        || (lfProbe != nullptr && lfProbe->isReady());

        if (! isFrozen() && processIndex != 0 && ready)
            gather();

        repaint();
    }

    /** Drains the colour out of a frozen trace, so it reads as a held frame at a glance
        rather than as a live one that happens not to be moving. */
    juce::Colour dim (juce::Colour colour) const
    {
        if (! isFrozen())
            return colour;

        return colour.withSaturation (0.0f).withMultipliedAlpha (0.45f);
    }

    void timerCallback() override
    {
        // Nothing to animate: the frame on screen is the one that should stay there, and
        // repainting fifteen times a second to draw it again would only cost battery.
        if (isFrozen() || processIndex == 0)
            return;

        const bool ready = (hfProbe != nullptr && hfProbe->isReady())
                        || (lfProbe != nullptr && lfProbe->isReady());

        if (! ready)
        {
            // Paint the empty state once rather than fifteen times a second.
            if (! paintedEmptyState)
            {
                paintedEmptyState = true;
                repaint();
            }

            return;
        }

        paintedEmptyState = false;

        // Gather here, once per frame, so paint() only ever draws. See paint().
        gather();
        repaint();
    }

    static constexpr double minHz = 20.0;
    static constexpr double maxHz = 20000.0;
    static constexpr double minLevelDb = -110.0;
    static constexpr double maxLevelDb = 6.0;

    // Symmetric, because the curve is a transfer rather than a boost: encode sits
    // above the zero line, decode's cut sits below it, and 0 dB is the centre.
    static constexpr double minGainDb = -28.0;
    static constexpr double maxGainDb = 28.0;
    int sections = lime::SrSideChain::defaultSections;

    /** The curve is built on a log grid rather than on either half's bins, so the two
        analysis resolutions can be combined onto one continuous trace. 480 points is more
        than a pixel apart at any window width the panel allows. */
    static constexpr int gridPoints = 480;
    static constexpr double twoPi = 6.283185307179586476925286766559;

    /** The grid the curve is built on, shared by the probe sampling and the fixed-response
        table so the two cannot land on different frequencies. */
    static double gridHz (int i)
    {
        const double t = (double) i / (double) (gridPoints - 1);
        return minHz * std::pow (maxHz / minHz, t);
    }

    struct Tick { double hz; const char* caption; };

    void drawGrid (juce::Graphics& g, juce::Rectangle<float> plot,
                   juce::Rectangle<int> axisRow, juce::Rectangle<int> scaleColumn)
    {
        using namespace lime::style;

        g.setFont (juce::FontOptions (fontSize * 0.55f));

        // Frequency axis, labelled: a log sweep is unreadable without it.
        for (auto tick : { Tick { 50.0, "50" },   Tick { 100.0, "100" },
                           Tick { 500.0, "500" }, Tick { 1000.0, "1k" },
                           Tick { 2000.0, "2k" }, Tick { 5000.0, "5k" },
                           Tick { 10000.0, "10k" } })
        {
            const float x = xFor (plot, tick.hz);

            g.setColour (text().withAlpha (0.12f));
            g.drawVerticalLine ((int) x, plot.getY(), plot.getBottom());

            g.setColour (muted());
            g.drawText (tick.caption,
                        juce::Rectangle<int> ((int) x - 20, axisRow.getY(), 40, axisRow.getHeight()),
                        juce::Justification::centred);
        }

        // The section boundaries: the log-spaced grid the fixed 800 Hz split between the two
        // staggered halves is sampled onto, so they show how coarsely or finely the process
        // is allocating gain across the band. Fainter than the single strong line this panel
        // used to draw at 800 Hz, because these are structure rather than a value to read off
        // the axis, but stronger than the decade gridlines, because they are the one thing on
        // this plot the Sections control changes.
        //
        // Only up to 64 of them. Past that the boundaries fall closer together than the
        // pixels between them, and instead of structure the eye gets an even accent wash
        // over the whole plot, which tells a reader less than drawing nothing would.
        if (sections <= 64)
        {
            // Taken from the DSP's own constants rather than this display's minHz/maxHz.
            // The two coincide today — both are 20 Hz to 20 kHz, which is why the
            // boundaries span the full plot width — but if one of them ever moved the ticks
            // have to follow where the sections actually are, not where the plot begins.
            const double lo = lime::SrSideChain::sectionLoHz;
            const double span = lime::SrSideChain::sectionHiHz / lime::SrSideChain::sectionLoHz;

            g.setColour (accent().withAlpha (0.22f));

            // Interior boundaries only. The outermost two land exactly on the plot's left
            // and right edges, where a line adds nothing to the frame already drawn there.
            for (int s = 1; s < sections; ++s)
            {
                const double hz = lo * std::pow (span, (double) s / (double) sections);
                g.drawVerticalLine ((int) xFor (plot, hz), plot.getY(), plot.getBottom());
            }
        }

        // Gain scale, so the curve's height means something. Symmetric about an
        // emphasised zero line: boost reads above it, the decode cut below.
        for (double db = -24.0; db <= 24.0; db += 8.0)
        {
            const float y = yFor (plot, db, minGainDb, maxGainDb);
            const bool zero = std::abs (db) < 0.5;

            if (zero)
            {
                g.setColour (text().withAlpha (0.55f));
                g.fillRect (plot.getX(), y - 0.75f, plot.getWidth(), 1.5f);
            }
            else
            {
                g.setColour (text().withAlpha (0.12f));
                g.drawHorizontalLine ((int) y, plot.getX(), plot.getRight());
            }

            g.setColour (muted());
            g.drawText (juce::String ((int) db),
                        juce::Rectangle<int> (scaleColumn.getX() + 3, (int) y - 8,
                                              scaleColumn.getWidth() - 3, 16),
                        juce::Justification::centredLeft);
        }

        // Units at the start of each scale rather than the end. They used to sit in the
        // bottom-right corner within a few pixels of each other and of the "0", which read
        // as one crowded smudge.
        g.setColour (muted());
        g.drawText ("dB",
                    juce::Rectangle<int> (scaleColumn.getX() + 3, (int) plot.getY() - 2,
                                          scaleColumn.getWidth() - 3, 14),
                    juce::Justification::centredLeft);
        g.drawText ("Hz",
                    juce::Rectangle<int> ((int) plot.getX(), axisRow.getY(),
                                          30, axisRow.getHeight()),
                    juce::Justification::centredLeft);
    }

    void drawLegend (juce::Graphics& g, juce::Rectangle<float> plot)
    {
        using namespace lime::style;

        auto row = plot.reduced (8.0f).removeFromTop (16.0f);
        g.setFont (juce::FontOptions (fontSize * 0.6f));

        const auto entry = [&] (juce::Colour colour, bool filled, const juce::String& caption)
        {
            auto swatch = row.removeFromLeft (18.0f).reduced (0.0f, 5.0f);

            g.setColour (colour);

            if (filled)
                g.fillRect (swatch);
            else
                g.fillRect (swatch.withSizeKeepingCentre (swatch.getWidth(), 2.0f));

            row.removeFromLeft (5.0f);

            g.setColour (text().withAlpha (0.85f));
            g.drawText (caption, row.removeFromLeft (100.0f), juce::Justification::centredLeft);

            row.removeFromLeft (10.0f);
        };

        entry (dim (base().brighter (0.5f)), true, "input");
        entry (dim (text()), false, "gain");
    }

    static float xFor (juce::Rectangle<float> area, double hz)
    {
        const double t = std::log (juce::jlimit (minHz, maxHz, hz) / minHz) / std::log (maxHz / minHz);
        return area.getX() + (float) t * area.getWidth();
    }

    static float yFor (juce::Rectangle<float> area, double value, double lo, double hi)
    {
        const double t = juce::jlimit (0.0, 1.0, (value - lo) / (hi - lo));
        return area.getBottom() - (float) t * area.getHeight();
    }

    /** One frame from one half, plus what is needed to sample it at an arbitrary frequency.

        The frame is claimed whole from the probe's triple buffer, and its geometry —
        bin count, FFT size, sample rate — travels inside the snapshot rather than
        being fetched from probe getters, so a sample-rate change can never mismatch
        a buffer against a bin count. */
    struct Frame
    {
        lime::SpectrumProbe::Snapshot snapshot;
        bool valid = false;

        void read (lime::SpectrumProbe* probe)
        {
            valid = probe != nullptr && probe->read (snapshot) && snapshot.bins > 1;
        }

        const std::vector<double>& levels() const noexcept { return snapshot.levelDb; }
        const std::vector<double>& gains() const noexcept { return snapshot.gainDb; }

        /** Linear interpolation between bins, in dB. Returns the fallback where this half
            has nothing to say. */
        double at (const std::vector<double>& values, double hz, double fallback) const
        {
            if (! valid || snapshot.fftSize <= 0 || snapshot.sampleRate <= 0.0 || values.size() < 2)
                return fallback;

            const double position = hz * (double) snapshot.fftSize / snapshot.sampleRate;

            if (position < 1.0 || position > (double) (values.size() - 1))
                return fallback;

            const auto lower = (size_t) position;
            const double fraction = position - (double) lower;

            return values[lower] + fraction * (values[lower + 1] - values[lower]);
        }

        /** Hz per bin, needed to put two different analysis resolutions on one scale. */
        double binWidthHz() const
        {
            return (valid && snapshot.fftSize > 0) ? snapshot.sampleRate / (double) snapshot.fftSize : 1.0;
        }
    };

    /** Combines both halves into one continuous curve on a log frequency grid.

        This used to splice — the long window below 800 Hz, the short one above — and that
        put a step of about 12 dB at the join into a plot of a process that has no step in
        it. Measure before smoothing: through the assembled encoder at 48 kHz with a tone at
        -85 dB, the real boost runs 19.88, 20.14, 20.38, 20.61, 20.82 dB through 700, 750,
        800, 850, 900 Hz, which is smooth and monotone. The step was the display's, because
        where the halves meet *each half on its own* is only doing part of the work, and
        showing one of them there shows about half the boost.

        So the halves are combined rather than chosen between. They are cascaded in the
        engine, so their boosts multiply and the decibels add. The input levels are power
        summed, since each half sees the programme through its own band split.

        Interpolating or splining across the join would have been the wrong fix: it would
        have drawn a smooth curve through the wrong values, and the reader would have had no
        way to tell.
    */
    void gather()
    {
        levelPoints.clearQuick();
        gainPoints.clearQuick();

        lf.read (lfProbe);
        hf.read (hfProbe);

        // Both halves are put on the short window's bin width, so broadband material does
        // not step where they meet purely because one analysis has finer bins than the
        // other and therefore less energy in each. It makes the input trace a spectral
        // density, which shifts a *tone* in the low half down by 10*log10(N_lf/N_hf) — a
        // fixed offset on an axis that carries no numbers, against a visible artefact at a
        // frequency the eye is already drawn to.
        const double reference = hf.valid ? hf.binWidthHz() : lf.binWidthHz();
        const double lfLevelOffset = lf.valid ? 10.0 * std::log10 (reference / lf.binWidthHz()) : 0.0;
        const double hfLevelOffset = hf.valid ? 10.0 * std::log10 (reference / hf.binWidthHz()) : 0.0;

        for (int i = 0; i < gridPoints; ++i)
        {
            const double hz = gridHz (i);

            // Cascaded, so the transfers multiply: decibels add. The fixed main path — the
            // spectral skewing and antisaturation networks — is part of the response
            // and is what makes it taper at both extremes, so it belongs in the curve. The
            // side chain alone would show a flat 24 dB out to 20 kHz, where the assembled
            // encoder actually delivers about 13 dB at 15 kHz and 9 dB at 30 Hz. In down
            // mode the fixed networks run inverted, so their contribution flips sign with
            // the side chains' cut.
            const double fixedSign = modeIndex == 1 ? -1.0 : 1.0;
            const double fixedDb = (i < (int) fixedResponseDb.size())
                                     ? fixedSign * fixedResponseDb[(size_t) i] : 0.0;
            const double gainDb = lf.at (lf.gains(), hz, 0.0) + hf.at (hf.gains(), hz, 0.0) + fixedDb;
            gainPoints.add ({ (float) hz, (float) gainDb });

            // Band split, so the levels add in power.
            const double power = std::pow (10.0, (lf.at (lf.levels(), hz, -300.0) + lfLevelOffset) / 10.0)
                               + std::pow (10.0, (hf.at (hf.levels(), hz, -300.0) + hfLevelOffset) / 10.0);

            levelPoints.add ({ (float) hz,
                               (float) (10.0 * std::log10 (std::max (1.0e-16, power))) });
        }
    }

    juce::Path buildPath (juce::Rectangle<float> area,
                          const juce::Array<juce::Point<float>>& points,
                          double lo, double hi) const
    {
        juce::Path path;

        for (int i = 0; i < points.size(); ++i)
        {
            const float x = xFor (area, points[i].x);
            const float y = yFor (area, points[i].y, lo, hi);

            if (i == 0)
                path.startNewSubPath (x, y);
            else
                path.lineTo (x, y);
        }

        return path;
    }

    lime::SpectrumProbe* hfProbe = nullptr;
    lime::SpectrumProbe* lfProbe = nullptr;

    bool paintedEmptyState = false;
    bool bypassed = false;

    // 2 is B, the process the plot was written for; the panel corrects this on its first tick.
    int processIndex = 2;

    // 2 is both; the panel corrects this on its first tick. 1 (down) flips the
    // fixed main path's contribution to match the inverted networks.
    int modeIndex = 2;

    Frame lf, hf;
    std::vector<double> fixedResponseDb;
    double fixedResponseRate = 0.0;
    juce::Array<juce::Point<float>> levelPoints, gainPoints;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SpectrumDisplay)
};
