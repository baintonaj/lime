/*
  ==============================================================================

    CompressionLaw — the knee shared by the fixed band and the sliding layer.

    Both stages compute gain = 1 / (1 + (control/threshold)^n), and they must use the
    same n or the two halves of a stage stop being commensurate. It lives here so
    neither can be changed without the other.

    The square is spelled out rather than left to std::pow. n is 2.0 in every
    constructed path, but it is a member rather than a literal, so std::pow could not
    see the exponent and emitted a full runtime call — per bin, per stage, per frame.
    That is around a thousand times the cost of the multiply it stands in for, and
    worse than the cost alone: a libm call in the middle of the frame loop is opaque
    to the vectoriser, so it pinned loops scalar whose per-bin work is otherwise
    entirely independent. Removing it is what lets those loops vectorise at all.

    pow(x, 2.0) and x * x are both the correctly rounded square, so this is exact
    rather than an approximation.

  ==============================================================================
*/

#pragma once

#include <cmath>

namespace lime
{

/** Raises the control-to-threshold ratio to the knee exponent.

    @param ratio     control signal over threshold, non-negative
    @param exponent  the knee; 1.0 and 2.0 are the shapes actually used
*/
inline double shapeRatio (double ratio, double exponent) noexcept
{
    if (exponent == 1.0)
        return ratio;

    if (exponent == 2.0)
        return ratio * ratio;

    return std::pow (ratio, exponent);
}

} // namespace lime
