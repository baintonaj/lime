/*
  ==============================================================================

    FixedBand — per-bin fixed band for one staggered stage.

  ==============================================================================
*/

#include "FixedBand.h"

#include "DspMath.h"

#include "CompressionLaw.h"

#include <algorithm>
#include <cmath>

namespace lime
{

namespace
{
    double dbToGain (double db) { return std::pow (10.0, db / 20.0); }
}

//==============================================================================
void FixedBand::setBandWeight (const double* bandWeight)
{
    for (int k = 0; k < numBins; ++k)
    {
        const double weight = bandWeight != nullptr ? bandWeight[k] : 1.0;
        transferScale[(size_t) k] = dbToGain (params.lowLevelGainDb * weight) - 1.0;
    }
}

void FixedBand::prepare (const FixedBandParams& p, int bins,
                         double sampleRate, int fftSize, double frameSeconds,
                         const double* bandWeight)
{
    params = p;
    numBins = bins;

    thresholdLinear = dbToGain (params.thresholdDb + params.precedingGainDb);

    transferScale.assign ((size_t) numBins, 0.0);
    setBandWeight (bandWeight);

    // At low level the band output equals the input, so the unnormalised
    // combination would be (1 + feedforwardWeight) times too large and would move
    // the knee. Normalising keeps the threshold where it was specified.
    controlNormalisation = 1.0 / (1.0 + params.feedforwardWeight);

    mainCoeff = framePoleCoefficient (params.mainTauMs, frameSeconds);
    passCoeff = framePoleCoefficient (params.passBandTauMs, frameSeconds);
    finalCoeff = framePoleCoefficient (params.finalTauMs, frameSeconds);

    passWeight.assign ((size_t) numBins, 0.0);

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * sampleRate / (double) fftSize;
        passWeight[(size_t) k] = onePoleMagnitude (hz, params.passBandCornerHz,
                                                  params.passBandIsHighPass);
    }

    mainSmoothed.assign ((size_t) numBins, 0.0);
    passSmoothed.assign ((size_t) numBins, 0.0);
    control.assign ((size_t) numBins, 0.0);
}

void FixedBand::reset()
{
    std::fill (mainSmoothed.begin(), mainSmoothed.end(), 0.0);
    std::fill (passSmoothed.begin(), passSmoothed.end(), 0.0);
    std::fill (control.begin(), control.end(), 0.0);
}

//==============================================================================
void FixedBand::processFrame (const double* magnitudes, double* transferOut)
{
    for (int k = 0; k < numBins; ++k)
    {
        const size_t i = (size_t) k;

        // Feedback control, using the control signal established by previous
        // frames: linear below the threshold, shedding boost above it at a rate set
        // by the knee exponent.
        const double ratio = control[i] / thresholdLinear;
        const double shaped = shapeRatio (ratio, params.kneeExponent);

        const double gain = 1.0 / (1.0 + shaped);

        const double bandOutput = gain * magnitudes[k];

        transferOut[k] = transferScale[i] * gain;

        // Feedback plus feedforward, normalised so the low-level knee stays put.
        // The feedforward term is what limits the differential component at high
        // level; see the comment on feedforwardWeight.
        const double combined = (bandOutput + params.feedforwardWeight * magnitudes[k])
                              * controlNormalisation;

        // Main control path: rectified control signal opposed by the end-stop.
        const double mainInput = std::max (0.0, combined - params.endStopFraction * magnitudes[k]);

        mainSmoothed[i] = mainCoeff * mainSmoothed[i] + (1.0 - mainCoeff) * mainInput;

        // Pass-band control path: frequency-weighted, and deliberately *not*
        // opposed by the end-stop, so that it can take over through the maximum
        // selector when the end-stop suppresses the main path more than the signal
        // conditions warrant.
        const double passInput = combined * passWeight[i];

        passSmoothed[i] = passCoeff * passSmoothed[i] + (1.0 - passCoeff) * passInput;

        const double selected = std::max (mainSmoothed[i], passSmoothed[i]);

        control[i] = finalCoeff * control[i] + (1.0 - finalCoeff) * selected;
    }
}

} // namespace lime
