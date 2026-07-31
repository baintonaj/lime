/*
  ==============================================================================

    DspMath — the two small helpers half this directory used to carry a private
    copy of.

    Three files each had their own single-pole magnitude and their own per-frame
    smoothing coefficient, spelled differently to dodge name collisions, and a
    set of near-duplicates is exactly how one of them eventually drifts onto a
    different convention. One definition each, one convention, stated here.

  ==============================================================================
*/

#pragma once

#include <cmath>

namespace lime
{

/** Magnitude of a single-pole high- or low-pass at a frequency.

    `cornerHz` must be positive; every caller passes a fixed design constant. */
inline double onePoleMagnitude (double hz, double cornerHz, bool highPass) noexcept
{
    const double ratio = hz / cornerHz;
    const double lowPass = 1.0 / std::sqrt (1.0 + ratio * ratio);

    return highPass ? lowPass * ratio : lowPass;
}

/** Per-frame one-pole coefficient for a time constant quoted in milliseconds.

    The convention is exp(-T/tau) — the fraction *kept* each frame — for use as
    state = a * state + (1 - a) * input. A non-positive tau returns zero, which
    makes the smoother instantaneous rather than stuck. (ControlSmoother quotes
    its coefficient the other way round, as the fraction moved per sample; that
    is a per-sample control-rate object with its own documented reading, not a
    fourth copy of this.) */
inline double framePoleCoefficient (double tauMs, double frameSeconds) noexcept
{
    const double tau = tauMs * 0.001;

    return tau <= 0.0 ? 0.0 : std::exp (-frameSeconds / tau);
}

} // namespace lime
