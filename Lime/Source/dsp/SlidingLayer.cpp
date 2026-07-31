/*
  ==============================================================================

    SlidingLayer — masking-shaped cross-bin profile for one staggered stage.

  ==============================================================================
*/

#include "SlidingLayer.h"

#include "CompressionLaw.h"
#include "DspMath.h"

#include <algorithm>
#include <cmath>

// The dB conversions below run per bin per frame in every stage of all four side
// chains — around five million log10 and pow calls a second at 48 kHz — so on
// Apple platforms they go through vForce, with the scalar loop as the portable
// reference. The same convention as Fft64.
#if ! defined (LIME_HAS_ACCELERATE) && defined (__APPLE__)
 #define LIME_HAS_ACCELERATE 1
#endif

#if LIME_HAS_ACCELERATE
 #include <Accelerate/Accelerate.h>
#endif

namespace lime
{

namespace
{
    constexpr double slidingSilenceDb = -400.0;
}

//==============================================================================
void SlidingLayer::prepare (const SlidingLayerParams& p, int bins,
                            double sampleRate, int fftSize, double frameSeconds)
{
    params = p;
    numBins = bins;

    thresholdLinear = std::pow (10.0, (params.thresholdDb + params.precedingGainDb) / 20.0);

    firstCoeff = framePoleCoefficient (params.firstTauMs, frameSeconds);
    finalCoeff = framePoleCoefficient (params.finalTauMs, frameSeconds);

    spread.prepare (numBins, sampleRate, fftSize);

    controlWeight.assign ((size_t) numBins, 0.0);

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * sampleRate / (double) fftSize;

        // High-frequency stages weight their control high-pass, low-frequency
        // stages low-pass, so each responds to dominants in its own region.
        controlWeight[(size_t) k] = onePoleMagnitude (hz, params.controlCornerHz,
                                                   params.highFrequency);
    }

    levelDb.assign ((size_t) numBins, slidingSilenceDb);
    thresholdDb.assign ((size_t) numBins, slidingSilenceDb);
    spreadLinear.assign ((size_t) numBins, 0.0);
    firstSmoothed.assign ((size_t) numBins, 0.0);
    control.assign ((size_t) numBins, 0.0);

    reset();
}

void SlidingLayer::reset()
{
    std::fill (firstSmoothed.begin(), firstSmoothed.end(), 0.0);
    std::fill (control.begin(), control.end(), 0.0);
    std::fill (spreadLinear.begin(), spreadLinear.end(), 0.0);
    std::fill (thresholdDb.begin(), thresholdDb.end(), slidingSilenceDb);

    modulationControl = 0.0;
    mcFirstSmoothed = 0.0;
    mcSmoothed = 0.0;
    numDominants = 0;
}

//==============================================================================
void SlidingLayer::processFrame (const double* magnitudes, double* fractionOut,
                                double* spreadLevelOut)
{
    // The masking spread function works in dB, so convert once. The magnitudes
    // are floored at the silence level's own linear value first, which folds the
    // zero test into the clamp: 20*log10(1e-20) is exactly the -400 dB floor.
   #if LIME_HAS_ACCELERATE
    {
        const int n = numBins;
        const double floorLinear = 1.0e-20;
        const double twenty = 20.0;

        vDSP_vthrD (magnitudes, 1, &floorLinear, spreadLinear.data(), 1, (vDSP_Length) n);
        vvlog10 (levelDb.data(), spreadLinear.data(), &n);
        vDSP_vsmulD (levelDb.data(), 1, &twenty, levelDb.data(), 1, (vDSP_Length) n);
    }
   #else
    for (int k = 0; k < numBins; ++k)
    {
        const double m = std::max (1.0e-20, magnitudes[k]);
        levelDb[(size_t) k] = 20.0 * std::log10 (m);
    }
   #endif

    // Per-bin threshold produced by every *other* bin. This is the skirt: steep
    // below a dominant, shallow and level-dependent above it, which is exactly the
    // asymmetry the analog's sliding band creates by sliding upward with level.
    spread.compute (levelDb.data(), thresholdDb.data());

    // The spread threshold back on the linear magnitude scale, in one vectorised
    // pass rather than a pow per bin inside the loop below: 10^x = e^(x ln10).
   #if LIME_HAS_ACCELERATE
    {
        const int n = numBins;
        const double ln10Over20 = 0.11512925464970228420089957273422;

        vDSP_vsmulD (thresholdDb.data(), 1, &ln10Over20, spreadLinear.data(), 1, (vDSP_Length) n);
        vvexp (spreadLinear.data(), spreadLinear.data(), &n);
    }
   #else
    for (int k = 0; k < numBins; ++k)
        spreadLinear[(size_t) k] = std::pow (10.0, thresholdDb[(size_t) k] / 20.0);
   #endif

    numDominants = 0;

    // The MC level passes through the same two smoothers as the per-bin control
    // before it opposes anything. The compressors are controlled by well-smoothed
    // signals, and the opposition has to be as smooth as what it opposes: the raw
    // weighted RMS moves at the frame rate, and dividing every bin's ratio by a
    // scalar that steps 190 times a second is a recipe for broadband pumping.
    // Smoothing it here with the identical coefficients also preserves what the
    // old formulation did implicitly — a subtraction ahead of a linear smoother
    // is a smoothed subtraction.
    mcFirstSmoothed = firstCoeff * mcFirstSmoothed + (1.0 - firstCoeff) * modulationControl;
    mcSmoothed = finalCoeff * mcSmoothed + (1.0 - finalCoeff) * mcFirstSmoothed;

    for (int k = 0; k < numBins; ++k)
    {
        const size_t i = (size_t) k;

        // A bin is a dominant when it stands above what its neighbours can mask.
        // An isolated tone sits far above its own threshold, so this discriminates
        // sharply; a bin buried in a skirt does not qualify.
        const bool isLocalMaximum = (k == 0 || levelDb[i] >= levelDb[i - 1])
                                 && (k == numBins - 1 || levelDb[i] >= levelDb[i + 1]);

        if (isLocalMaximum && levelDb[i] > thresholdDb[i])
            ++numDominants;

        // The spread level, weighted by this stage's control filter so the stage
        // only responds to dominants in its own frequency region.
        const double spreadMagnitude = spreadLinear[i];

        if (spreadLevelOut != nullptr)
            spreadLevelOut[k] = spreadMagnitude;

        const double controlInput = spreadMagnitude * controlWeight[i];

        firstSmoothed[i] = firstCoeff * firstSmoothed[i] + (1.0 - firstCoeff) * controlInput;
        control[i] = finalCoeff * control[i] + (1.0 - finalCoeff) * firstSmoothed[i];

        // Opposed by MC1 (high frequencies) or MC4 (low). Section 3.3: "the ratio
        // between the rectified control signal and MC1 monitors the signal
        // attenuation", and that ratio is what forms the end-stop. So the MC level
        // raises the reference the knee measures against rather than being
        // subtracted per bin: as the overall weighted level rises, every bin's
        // ratio falls together and the withholding releases smoothly — the "no
        // reason for further sliding" condition of section 2.3. A subtraction
        // cannot express that: MC is one scalar for the whole spectrum, and taking
        // it off each bin's own control zeroes the layer outright the moment the
        // programme is broadband. At MC = 0 the reference is the bare threshold
        // and the law is exactly the fixed band's.
        //
        // Carried through the stage's own control weighting, as the subtraction
        // was: the opposition belongs to the frequency region the stage responds
        // to, not uniformly to every bin the transform happens to have.
        const double reference = thresholdLinear
                               + params.mcOppositionGain * mcSmoothed * controlWeight[i];

        // Same compression law as the fixed band, so the two are commensurate.
        const double ratio = control[i] / reference;
        const double shaped = shapeRatio (ratio, params.kneeExponent);

        const double pass = 1.0 / (1.0 + shaped);

        if (params.composition == SlidingComposition::restriction)
        {
            // Surviving boost fraction, with the end-stop limiting how much any
            // dominant may strip.
            const double withheld = std::min (params.maximumWithheldFraction, 1.0 - pass);
            fractionOut[k] = 1.0 - withheld;
        }
        else
        {
            // This layer's own boost fraction, to be combined by Eq. (1).
            fractionOut[k] = pass;
        }
    }
}

} // namespace lime
