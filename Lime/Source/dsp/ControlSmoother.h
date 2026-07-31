/*
  ==============================================================================

    ControlSmoother — one pole between a control and the signal path.

    A control read once per block and applied as a constant puts a step into the signal
    every time it moves. A step in a gain is a click, and with a 64-sample block at 48 kHz
    a knob drag produces one every 1.3 ms. Feeding the value through a single pole makes
    the change continuous instead.

    The settling time is quoted the way a listener experiences it: **the time to cover 95%
    of the distance**, which for a one-pole is three time constants. So 75 ms of settling is
    tau = 25 ms. Quoting tau = 75 ms instead would mean 225 ms to arrive, which no longer
    tracks a hand on a knob — the reading matters, so it is stated rather than implied.

    Deliberately not an exponential-only smoother of a gain in decibels: the smoothing is
    applied to whatever quantity the caller passes, and the callers here smooth linear gain,
    so a trim moving from -10 dB to +10 dB passes through unity rather than through the
    geometric mean. For a calibration trim that is the right choice, because the number on
    the panel is a level and the ear follows level.

  ==============================================================================
*/

#pragma once

#include <cmath>

namespace lime
{

/** The settling time every control in the plugin uses, in seconds to 95%.

    Lives here rather than in the wrapper because the engine has trims of its own now —
    the three inside loop mode — and two constants that have to agree are one constant
    waiting to disagree. */
inline constexpr double defaultSettlingSeconds = 0.075;

class ControlSmoother
{
public:
    /** @param sampleRate       samples per second
        @param settlingSeconds  time to cover 95% of a step; see the note above
        @param initial          starting value, which becomes both current and target, so
                                the first block after prepare does not ramp from zero
    */
    void prepare (double sampleRate, double settlingSeconds, double initial)
    {
        // Three time constants to 95%, so tau = settling / 3.
        const auto tau = settlingSeconds / 3.0;

        coefficient = (sampleRate > 0.0 && tau > 0.0)
                        ? 1.0 - std::exp (-1.0 / (tau * sampleRate))
                        : 1.0;

        snapTo (initial);
    }

    void setTarget (double newTarget) noexcept   { target = newTarget; }
    double getTarget() const noexcept            { return target; }

    /** Jumps, for prepareToPlay and for a host loading a session — a saved state should be
        in force on the first sample, not faded in over 75 ms. */
    void snapTo (double value) noexcept          { target = current = value; }

    double getCurrent() const noexcept           { return current; }

    /** One sample of movement toward the target. */
    double next() noexcept
    {
        current += coefficient * (target - current);
        return current;
    }

    /** True while the value is still moving, so a caller can skip the per-sample loop
        entirely once it has arrived. The tolerance is well below anything audible in a
        gain — 1e-9 is about -180 dB. */
    bool isMoving (double tolerance = 1.0e-9) const noexcept
    {
        return std::abs (target - current) > tolerance;
    }

private:
    double coefficient = 1.0;
    double current = 0.0;
    double target = 0.0;
};

} // namespace lime
