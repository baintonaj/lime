/*
  ==============================================================================

    FractionalDelay64 — a delay of any length, not just a whole number of samples.

    DelayLine64 is deliberately integer-only: everything it aligns is an integer number of
    samples by construction, and interpolating would only add error to a delay that has an
    exact answer. This is the other case. A delay a user sets in hundredths of a millisecond
    is 0.48 samples at 48 kHz, and rounding it to whole samples would make the bottom two
    digits of the control do nothing at all.

    Third-order Lagrange interpolation, over four taps. Linear interpolation was the obvious
    cheaper choice and is not good enough here: its response falls away with frequency and
    with the fractional part, so a half-sample delay loses several dB at the top of the band
    and none at all at zero — which, on a control whose whole purpose is to be blended
    against an undelayed copy, reads as the delay dulling the sound rather than moving it.
    Lagrange at third order holds the magnitude flat to well under a tenth of a dB to 15 kHz.

    The delay is smoothed, so it can be swept. That makes it a varying delay line while it
    moves, which is to say it pitches slightly, exactly as moving a physical distance would.

    An engaged delay is floored at one sample — 0.021 ms at 48 kHz. The kernel's earliest
    tap sits one sample ahead of the read position, and ahead of the sample just written
    the ring holds nothing but the oldest audio it has, so a read position inside the first
    sample would blend in whatever went through tens of milliseconds ago. A delay set to
    exactly zero never reaches the interpolator at all: process() bypasses the line
    entirely, so the floor costs nothing where zero is actually wanted.

  ==============================================================================
*/

#pragma once

#include "ControlSmoother.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace lime
{

class FractionalDelay64
{
public:
    void prepare (double sampleRate, double maxDelaySeconds)
    {
        // The interpolator's deepest tap reaches two samples past the whole delay, and the
        // slot the write index is about to reuse must stay clear of it, so three extra
        // samples are strictly required. Eight is that requirement with margin; a few
        // spare doubles cost nothing against the price of an off-by-one here.
        maxSamples = std::max (1.0, std::ceil (maxDelaySeconds * sampleRate));
        capacity = (int) maxSamples + 8;

        buffer.assign ((size_t) capacity, 0.0);
        writeIndex = 0;

        smoothed.prepare (sampleRate, defaultSettlingSeconds, 0.0);
    }

    void reset()
    {
        std::fill (buffer.begin(), buffer.end(), 0.0);
        writeIndex = 0;
        smoothed.snapTo (smoothed.getTarget());
    }

    /** @param samples  delay in samples, fractional. Clamped to what was prepared for. */
    void setDelaySamples (double samples) noexcept
    {
        smoothed.setTarget (std::clamp (samples, 0.0, maxSamples));
    }

    void snapDelaySamples (double samples) noexcept
    {
        smoothed.snapTo (std::clamp (samples, 0.0, maxSamples));
    }

    /** Delays a block in place, or does nothing at all if it is set to zero.

        The skip is the reason this is a block call rather than a bare per-sample one. A
        line that is not being written holds whatever was last put through it, so simply
        resuming would replay however many milliseconds of stale audio the buffer happened
        to be sitting on — the same trap DelayLine64 fell into when the engine's padding
        changed. Waking up therefore clears the buffer first: the first samples out are
        silence rather than something from minutes ago, and at the moment of waking the
        delay is still ramping up from zero, so there is nothing there to hear anyway. */
    void process (double* samples, int numSamples) noexcept
    {
        // Dormant means "the control is heading for zero and the line has nothing
        // left to say", not "the smoother's value is exactly zero" — a one-pole
        // never actually arrives, so testing the value would leave the line running
        // (and floored at one sample) for the many seconds the residual takes to
        // underflow after the knob returns to rest. Below half a sample the floor
        // has already pinned the read position to one sample, so the single-sample
        // splice to passthrough is the same whenever it is taken; taking it here,
        // as the decay crosses the halfway point, is what keeps a released control
        // from leaving a lingering offset against the dry path.
        if (smoothed.getTarget() <= 0.0 && smoothed.getCurrent() < 0.5)
        {
            dormant = true;
            return;
        }

        if (dormant)
        {
            std::fill (buffer.begin(), buffer.end(), 0.0);
            writeIndex = 0;
            dormant = false;
        }

        for (int n = 0; n < numSamples; ++n)
            samples[n] = processSample (samples[n]);
    }

private:
    double processSample (double x) noexcept
    {
        buffer[(size_t) writeIndex] = x;

        // Floored at one sample, because the kernel's earliest tap sits one sample ahead
        // of the read position: below a whole sample that tap would land ahead of the
        // write index, where the ring keeps its oldest audio, not its newest. Exact zero
        // never gets here — process() bypasses the line for it.
        const double requested = smoothed.next();
        const double readPos = std::max (1.0, requested);

        // The four taps straddle the read position: one before it and two after, which is
        // what the third-order kernel is centred on.
        const int whole = (int) readPos;
        const double fraction = readPos - (double) whole;

        const double y = lagrange (readAt (whole - 1), readAt (whole),
                                   readAt (whole + 1), readAt (whole + 2), fraction);

        if (++writeIndex >= capacity)
            writeIndex = 0;

        return y;
    }

    double readAt (int samplesBack) const noexcept
    {
        int index = writeIndex - samplesBack;

        while (index < 0)
            index += capacity;

        while (index >= capacity)
            index -= capacity;

        return buffer[(size_t) index];
    }

    /** Third-order Lagrange, with `t` measured from the second of the four points. The
        points arrive newest-first, so `t` runs from x1 toward x2 as the delay lengthens. */
    static double lagrange (double xm1, double x0, double x1, double x2, double t) noexcept
    {
        const double c0 = x0;
        const double c1 = x1 - (1.0 / 3.0) * xm1 - 0.5 * x0 - (1.0 / 6.0) * x2;
        const double c2 = 0.5 * (xm1 + x1) - x0;
        const double c3 = (1.0 / 6.0) * (x2 - xm1) + 0.5 * (x0 - x1);

        return ((c3 * t + c2) * t + c1) * t + c0;
    }

    std::vector<double> buffer;
    ControlSmoother smoothed;

    double maxSamples = 1.0;
    int capacity = 1, writeIndex = 0;
    bool dormant = true;
};

} // namespace lime
