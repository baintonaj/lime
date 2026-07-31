/*
  ==============================================================================

    AutoGain — loudness-matching makeup for the processed path.

    A dynamics unit that lifts quiet detail is louder than its input by however
    much lifting it did, and comparing the two by ear then measures loudness, not
    processing. This stage matches the wet path's long-term power back to the dry
    path's, so bypass, mix and the process switches compare character against
    character.

    Deliberately slow and bounded rather than clever:

      - Both sides' powers run through the same multi-second one-pole, so the
        correction tracks programme loudness, not envelopes — it is makeup, not a
        compressor stacked on a compressor. The correction target itself is then
        smoothed again on its way to the signal, so nothing steps.
      - The correction is clamped to +/-limitDb. A pathological ratio — silence on
        one side, a fault on the other — can never run the gain away.
      - Below the gate threshold the correction holds. Silence carries no loudness
        to match, and the previous programme's correction is the best available
        estimate for what comes next.

    Measured on the wet path *before* the dry blend, applied in the same place, so
    at mix extremes the arithmetic still compares like with like.

  ==============================================================================
*/

#pragma once

#include "ControlSmoother.h"

#include <algorithm>
#include <cmath>

namespace lime
{

class AutoGain
{
public:
    static constexpr double limitDb = 12.0;

    void prepare (double sampleRate)
    {
        currentSampleRate = sampleRate;

        // Loudness integration: several seconds to 95%, so bars and phrases read
        // as one level rather than pumping the correction.
        const double tau = integrationSeconds / 3.0;
        blockCoeffBase = tau > 0.0 ? 1.0 / (tau * sampleRate) : 1.0;

        correction.prepare (sampleRate, correctionSettlingSeconds, 1.0);
        reset();
    }

    void reset()
    {
        dryPower = 0.0;
        wetPower = 0.0;
        correction.snapTo (1.0);
        engaged = false;
    }

    void setEnabled (bool shouldBeEnabled) noexcept { engaged = shouldBeEnabled; }
    bool isEnabled() const noexcept                 { return engaged; }

    /** The correction currently applied, for a readout. */
    double getCorrectionDb() const noexcept
    {
        return 20.0 * std::log10 (std::max (1.0e-12, correction.getCurrent()));
    }

    /** Accumulates one block of both sides' loudness. Call before apply(). */
    void measure (const double* const* dry, const double* const* wet,
                  int numChannels, int numSamples) noexcept
    {
        if (numChannels <= 0 || numSamples <= 0)
            return;

        double drySum = 0.0, wetSum = 0.0;

        for (int ch = 0; ch < numChannels; ++ch)
        {
            for (int n = 0; n < numSamples; ++n)
            {
                drySum += dry[ch][n] * dry[ch][n];
                wetSum += wet[ch][n] * wet[ch][n];
            }
        }

        const double norm = 1.0 / (double) (numChannels * numSamples);
        const double a = std::min (1.0, blockCoeffBase * (double) numSamples);

        dryPower += a * (drySum * norm - dryPower);
        wetPower += a * (wetSum * norm - wetPower);
    }

    /** Applies the current correction to the wet path. */
    void apply (double* const* wet, int numChannels, int numSamples) noexcept
    {
        if (engaged && dryPower > gatePower && wetPower > gatePower)
        {
            const double limit = std::pow (10.0, limitDb / 20.0);
            const double wanted = std::sqrt (dryPower / wetPower);

            correction.setTarget (std::clamp (wanted, 1.0 / limit, limit));
        }
        else if (! engaged)
        {
            correction.setTarget (1.0);
        }
        // Engaged but gated: hold the last correction — silence has no loudness
        // to match.

        if (! correction.isMoving())
        {
            const double gain = correction.getCurrent();

            if (std::abs (gain - 1.0) < 1.0e-9)
                return;

            for (int ch = 0; ch < numChannels; ++ch)
                for (int n = 0; n < numSamples; ++n)
                    wet[ch][n] *= gain;

            return;
        }

        for (int n = 0; n < numSamples; ++n)
        {
            const double gain = correction.next();

            for (int ch = 0; ch < numChannels; ++ch)
                wet[ch][n] *= gain;
        }
    }

private:
    static constexpr double integrationSeconds = 3.0;
    static constexpr double correctionSettlingSeconds = 0.4;

    /** -70 dBFS in power: below this neither side carries a loudness worth
        matching, and the correction holds. */
    static constexpr double gatePower = 1.0e-7;

    ControlSmoother correction;

    double currentSampleRate = 48000.0;
    double blockCoeffBase = 1.0;
    double dryPower = 0.0;
    double wetPower = 0.0;
    bool engaged = false;
};

} // namespace lime
