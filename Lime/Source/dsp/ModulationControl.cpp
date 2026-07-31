/*
  ==============================================================================

    ModulationControl — MC1 to MC8.

  ==============================================================================
*/

#include "ModulationControl.h"

#include "DspMath.h"

#include <algorithm>
#include <cmath>

namespace lime
{

namespace
{
    /** Power-combined weighted level over the spectrum, normalised to the weight's
        own norm — the weighted RMS. See the note in the header on why the
        normalisation matters. */
    double weightedLevel (const double* magnitudes, const std::vector<double>& weights,
                          double weightNorm)
    {
        double sum = 0.0;

        for (size_t k = 0; k < weights.size(); ++k)
        {
            const double v = weights[k] * magnitudes[k];
            sum += v * v;
        }

        return std::sqrt (sum) / weightNorm;
    }

    /** sqrt (sum of w^2), floored so an all-zero weighting cannot divide by zero. */
    double normOf (const std::vector<double>& weights)
    {
        double sum = 0.0;

        for (const double w : weights)
            sum += w * w;

        return std::max (1.0e-12, std::sqrt (sum));
    }
}

//==============================================================================
void ModulationControl::TwoStageSmoother::setTau (double tauMs, double frameSeconds)
{
    coeff = framePoleCoefficient (tauMs, frameSeconds);
}

double ModulationControl::TwoStageSmoother::process (double input)
{
    a = coeff * a + (1.0 - coeff) * input;
    b = coeff * b + (1.0 - coeff) * a;

    return b;
}

//==============================================================================
void ModulationControl::prepare (int bins, double sampleRate, int fftSize, double frameSeconds)
{
    numBins = bins;

    const auto build = [bins, sampleRate, fftSize] (std::vector<double>& table,
                                                    auto&& shape)
    {
        table.assign ((size_t) bins, 0.0);

        for (int k = 0; k < bins; ++k)
        {
            const double hz = (double) k * sampleRate / (double) fftSize;
            table[(size_t) k] = shape (hz);
        }
    };

    using namespace mcWeighting;

    build (weightHfSliding, [] (double hz)
    {
        return onePoleMagnitude (hz, hfSlidingHighPassHz, true);
    });

    build (weightHfFixed, [] (double hz)
    {
        return onePoleMagnitude (hz, hfFixedLowPassOneHz, false)
             * onePoleMagnitude (hz, hfFixedLowPassTwoHz, false);
    });

    build (weightLfSliding, [] (double hz)
    {
        return onePoleMagnitude (hz, lfSlidingLowPassHz, false);
    });

    build (weightLfFixed, [] (double hz)
    {
        return onePoleMagnitude (hz, lfFixedHighPassOneHz, true)
             * onePoleMagnitude (hz, lfFixedHighPassTwoHz, true);
    });

    build (weightTransient, [] (double hz)
    {
        return onePoleMagnitude (hz, transientHighPassHz, true);
    });

    normHfSliding = normOf (weightHfSliding);
    normHfFixed = normOf (weightHfFixed);
    normLfSliding = normOf (weightLfSliding);
    normLfFixed = normOf (weightLfFixed);
    normTransient = normOf (weightTransient);

    smootherMc2.setTau (hfSlidingOvershootTauMs, frameSeconds);
    smootherMc5.setTau (lfSlidingOvershootTauMs, frameSeconds);
    smootherMc7.setTau (lfFixedOvershootTauMs, frameSeconds);

    transientDecay = framePoleCoefficient (transientPeakHoldMs, frameSeconds);

    reset();
}

void ModulationControl::reset()
{
    smootherMc2.reset();
    smootherMc5.reset();
    smootherMc7.reset();

    previousTransientLevel = 0.0;
    transientHold = 0.0;

    values.fill (0.0);
}

//==============================================================================
void ModulationControl::processFrame (const double* stageFeed, const double* skewFeed)
{
    const double mc1 = weightedLevel (stageFeed, weightHfSliding, normHfSliding);
    const double mc3 = weightedLevel (stageFeed, weightHfFixed, normHfFixed);
    const double mc4 = weightedLevel (stageFeed, weightLfSliding, normLfSliding);
    const double mc6 = weightedLevel (stageFeed, weightLfFixed, normLfFixed);

    values[0] = mc1;
    values[1] = smootherMc2.process (mc1);
    values[2] = mc3;
    values[3] = mc4;
    values[4] = smootherMc5.process (mc4);
    values[5] = mc6;
    values[6] = smootherMc7.process (mc6);

    // MC8: high-frequency transient detection. The analog double-differentiates
    // with 15 us constants, far below the frame interval, so the rising part of the
    // frame-to-frame change stands in for it. The 30 ms peak hold is as specified.
    const double transientLevel = weightedLevel (skewFeed, weightTransient, normTransient);
    const double rise = std::max (0.0, transientLevel - previousTransientLevel);

    previousTransientLevel = transientLevel;

    transientHold = std::max (rise, transientHold * transientDecay);
    values[7] = transientHold;
}

double ModulationControl::processFrameEssential (const double* stageFeed, bool highFrequency)
{
    if (highFrequency)
    {
        const double mc1 = weightedLevel (stageFeed, weightHfSliding, normHfSliding);
        values[0] = mc1;
        values[1] = smootherMc2.process (mc1);
        return mc1;
    }

    const double mc4 = weightedLevel (stageFeed, weightLfSliding, normLfSliding);
    values[3] = mc4;
    values[4] = smootherMc5.process (mc4);
    return mc4;
}

} // namespace lime
