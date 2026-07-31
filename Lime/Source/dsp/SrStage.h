/*
  ==============================================================================

    SrStage — one staggered stage: a per-bin fixed band and a sliding layer.

    Section 2.2: "both fixed-band and sliding-band dynamic actions are used in each
    of the five stages, a total of ten compressor circuits. While there is an
    effective interaction of the fixed and sliding bands in any particular stage,
    all of the stages operate independently."

    The two are combined by restriction rather than by Eq. (1). That choice was
    settled by measurement, not preference: applying the augmenting law leaves a
    dominant bin at full boost, because a masking threshold excludes the bin it is
    computed for, so the sliding circuit reports no reason to compress exactly where
    the fixed band has correctly decided to. See DESIGN.md and SlidingComposition.

  ==============================================================================
*/

#pragma once

#include "FixedBand.h"
#include "SlidingLayer.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace lime
{

struct SrStageParams
{
    /** Section 2.5: thresholds of approximately -30, -48 and -62 dB. */
    double thresholdDb = -30.0;
    double precedingGainDb = 0.0;
    double lowLevelGainDb = 8.0;

    bool highFrequency = true;

    SlidingComposition composition = SlidingComposition::restriction;
};

//==============================================================================
class SrStage
{
public:
    /** @param bandWeight  per-bin fraction of this stage's boost, in [0,1]. The two
                             frequency halves use power-complementary single-pole
                             weights so their dB contributions add across 800 Hz. */
    void prepare (const SrStageParams& p, int bins,
                  double sampleRate, int fftSize, double frameSeconds,
                  const double* bandWeight = nullptr)
    {
        params = p;
        numBins_ = bins;

        FixedBandParams fp;
        fp.thresholdDb = params.thresholdDb;
        fp.precedingGainDb = params.precedingGainDb;
        fp.lowLevelGainDb = params.lowLevelGainDb;
        fp.passBandCornerHz = params.highFrequency ? 1600.0 : 400.0;
        fp.passBandIsHighPass = params.highFrequency;
        fp.finalTauMs = params.highFrequency ? 160.0 : 300.0;
        fixed.prepare (fp, numBins_, sampleRate, fftSize, frameSeconds, bandWeight);
        transferScale = fixed.getTransferScale();

        SlidingLayerParams sp;
        sp.highFrequency = params.highFrequency;
        sp.thresholdDb = params.thresholdDb;
        sp.precedingGainDb = params.precedingGainDb;
        sp.lowLevelGainDb = params.lowLevelGainDb;
        sp.controlCornerHz = params.highFrequency ? 10000.0 : 80.0;
        sp.firstTauMs = params.highFrequency ? 5.0 : 7.5;
        sp.finalTauMs = params.highFrequency ? 80.0 : 150.0;
        sp.composition = params.composition;
        sliding.prepare (sp, numBins_, sampleRate, fftSize, frameSeconds);

        fixedTransfer.assign ((size_t) numBins_, 0.0);
        slidingFraction.assign ((size_t) numBins_, 0.0);
        spreadLevel.assign ((size_t) numBins_, 0.0);
        fixedFeed.assign ((size_t) numBins_, 0.0);
    }

    /** A new band weighting, which is what changing the section count comes to. Allocation
        free; the transfer scale pointer is unchanged, since only its contents move. */
    void setBandWeight (const double* bandWeight)
    {
        fixed.setBandWeight (bandWeight);
    }

    void reset()
    {
        fixed.reset();
        sliding.reset();
    }

    /** MC1 for the high-frequency stages, MC4 for the low. */
    void setModulationControl (double level) noexcept
    {
        sliding.setModulationControl (level);
    }

    /** The user's attack/release override, to both circuits' final smoothers at
        once — the stage is one detection and its two halves must keep moving
        together. See FixedBand::setTimeConstants for the contract. */
    void setTimeConstants (double attackMs, double releaseMs)
    {
        fixed.setTimeConstants (attackMs, releaseMs);
        sliding.setTimeConstants (attackMs, releaseMs);
    }

    int getNumBins() const noexcept { return numBins_; }
    int getNumDominants() const noexcept { return sliding.getNumDominants(); }

    /** Advances one frame.

        @param magnitudes   this stage's input magnitude per bin
        @param transferOut  per-bin transfer F, so the stage output is
                            input * (1 + F)
    */
    void processFrame (const double* magnitudes, double* transferOut)
    {
        sliding.processFrame (magnitudes, slidingFraction.data(), spreadLevel.data());

        // Action substitution mode also drives the fixed band from the spread level,
        // so that both circuits see the masked region as occupied. Retained for
        // comparison only; restriction is the default.
        const double* feed = magnitudes;

        if (params.composition == SlidingComposition::actionSubstitution)
        {
            for (int k = 0; k < numBins_; ++k)
                fixedFeed[(size_t) k] = std::max (magnitudes[k], spreadLevel[(size_t) k]);

            feed = fixedFeed.data();
        }

        fixed.processFrame (feed, fixedTransfer.data());

        for (int k = 0; k < numBins_; ++k)
        {
            const double scale = transferScale[(size_t) k];

            if (scale <= 0.0)
            {
                transferOut[k] = 0.0;
                continue;
            }

            const double phiFixed = fixedTransfer[(size_t) k] / scale;
            const double phiSliding = slidingFraction[(size_t) k];

            const double phi = params.composition == SlidingComposition::restriction
                                 ? phiFixed * phiSliding
                                 : 1.0 - (1.0 - phiFixed) * (1.0 - phiSliding);

            transferOut[k] = scale * phi;
        }
    }

private:
    SrStageParams params;

    int numBins_ = 0;
    const double* transferScale = nullptr;   // owned by `fixed`

    FixedBand fixed;
    SlidingLayer sliding;

    std::vector<double> fixedTransfer, slidingFraction, spreadLevel, fixedFeed;
};

} // namespace lime
