/*
  ==============================================================================

    Crossfade — a finite, shaped ramp between two states.

    ControlSmoother is the right thing for a control that is being moved: a one-pole
    tracks a hand on a knob and never has to arrive. It is the wrong thing for a
    transition, for two reasons the panel made audible.

    First, a one-pole is steepest at the moment it starts. Covering 63% of the distance
    in the first time constant is what makes it feel responsive under a knob, and exactly
    what makes a switch sound like a switch: the ear hears the leading edge, not the long
    tail. A raised cosine has zero slope at both ends, so a transition eases out of one
    state and into the other with no corner at either.

    Second, an exponential never actually arrives, so code waiting to act "once it has
    faded out" has to test against a small number and accept whatever time that takes.
    Here the ramp covers the whole distance in exactly the time it was configured with,
    and the caller can throw its switch knowing the wet path is at zero rather than
    near it.

  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <cmath>

namespace lime
{

class Crossfade
{
public:
    /** @param sampleRate  samples per second
        @param seconds     time to cross the whole range, in one direction
        @param initial     starting position, 0 or 1
    */
    void prepare (double sampleRate, double seconds, double initial)
    {
        increment = (sampleRate > 0.0 && seconds > 0.0) ? 1.0 / (seconds * sampleRate) : 1.0;
        snapTo (initial);
    }

    void setTarget (double newTarget) noexcept { target = newTarget > 0.5 ? 1.0 : 0.0; }
    double getTarget() const noexcept          { return target; }

    void snapTo (double value) noexcept
    {
        phase = target = std::clamp (value, 0.0, 1.0);
        shaped = shape (phase);
    }

    /** The shaped position, which is what callers should multiply by. */
    double getCurrent() const noexcept { return shaped; }

    /** One sample of movement toward the target. */
    double next() noexcept
    {
        // A fade spends nearly all of its life parked at one end or the other, where the
        // cosine has nothing new to say — so it is only evaluated while the phase moves.
        if (! isMoving())
            return shaped;

        if (phase < target)
            phase = std::min (target, phase + increment);
        else
            phase = std::max (target, phase - increment);

        shaped = shape (phase);
        return shaped;
    }

    /** Exact, unlike the one-pole's: the ramp reaches its endpoint and stops. The
        comparison is exact on purpose — the ramp lands on the target with min/max,
        so equality is reachable — spelled as two orderings only to say so without
        tripping float-equality warnings. */
    bool isMoving() const noexcept { return phase < target || target < phase; }

private:
    /** Raised cosine: zero slope at 0 and at 1, so neither end has a corner. */
    static double shape (double t) noexcept
    {
        return 0.5 - 0.5 * std::cos (3.14159265358979323846 * t);
    }

    double increment = 1.0;
    double phase = 0.0;
    double target = 0.0;

    // shape(phase), kept current whenever phase changes, so the settled states cost a
    // load rather than a cosine.
    double shaped = 0.0;
};

} // namespace lime
