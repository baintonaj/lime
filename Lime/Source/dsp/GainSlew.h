/*
  ==============================================================================

    GainSlew — per-bin limit on how fast the transfer may change.

    Stands in for the analog overshoot suppressors. Those cannot be reproduced
    directly: section 2.4 builds them from unsmoothed rectified control signals with
    15 us double differentiators, three hundred times shorter than the 5.33 ms frame
    interval at 48 kHz. What they exist to prevent also does not arise here — an STFT
    side chain cannot overshoot into the channel the way an analog one can.

    What does arise instead is pre-echo: a gain that jumps between frames smears
    across the whole analysis window, before as well as after the event. Limiting the
    per-frame rate of change is the STFT's equivalent problem, so that is what this
    does: a plain clamp on each bin's frame-to-frame change in dB, sized so that for
    most material it never acts and the compressors remain "controlled by
    well-smoothed signals".

    One feature of the analog arrangement is kept because it carries over directly:
    the limits track the modulation control signal. Section 2.4 derives the overshoot
    thresholds from the same signals as the steady-state ones "whereby a tracking
    action between the transient and steady-state behaviour is obtained", and the
    same coupling is wanted here — the limiter should relax under exactly the
    conditions the steady-state circuits are already handling.

    Asymmetric by design. Boost increasing means gain reduction releasing, which is
    the recovery direction and where the paper wants speed: action substitution and
    modulation control are both credited with making recovery fast. Boost decreasing
    means gain reduction engaging, which is where a jump smears.

  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace lime
{

struct GainSlewParams
{
    /** Maximum change per frame, in dB, when gain reduction is engaging. */
    double engagingDbPerFrame = 3.0;

    /** Maximum change per frame when it is releasing. Looser, because recovery speed
        is something the paper explicitly values. */
    double releasingDbPerFrame = 6.0;

    /** How strongly the modulation control signal raises the limits. Section 2.4
        derives the overshoot thresholds from the same signals as the steady-state
        ones "whereby a tracking action between the transient and steady-state
        behaviour is obtained". */
    double modulationControlScaling = 6.0;
};

//==============================================================================
class GainSlew
{
public:
    void prepare (const GainSlewParams& p, int numBins)
    {
        params = p;
        numBins_ = numBins;
        previousFactor.assign ((size_t) numBins, 1.0);
        primed = false;
    }

    void reset()
    {
        std::fill (previousFactor.begin(), previousFactor.end(), 1.0);
        primed = false;
    }

    /** Raises the limits in proportion to the modulation control level, so the
        limiter relaxes under the conditions the steady-state circuits are already
        handling. */
    void setModulationControl (double level) noexcept
    {
        modulationControl = std::max (0.0, level);
    }

    int getNumBins() const noexcept { return numBins_; }

    /** Clamps `transfer` in place. Values are the per-bin transfer F, so the applied
        factor is 1 + F. */
    void process (double* transfer)
    {
        const double relax = 1.0 + params.modulationControlScaling * modulationControl;
        const double engaging = params.engagingDbPerFrame * relax;
        const double releasing = params.releasingDbPerFrame * relax;

        // The state is kept as a linear factor rather than in decibels, because a per-frame
        // clamp in dB does not need decibels at all: |20*log10(f/p)| <= limit is
        // p*lower <= f <= p*upper, and both bounds are constant across the frame.
        //
        // That matters because the limiter is meant to be inert for ordinary material —
        // the limits sit far above the frame-to-frame movement of a well-smoothed control
        // signal — so the overwhelming majority of bins pass the comparison untouched, and
        // a bin that does trip lands exactly on a bound: a multiply, never a log or a pow.
        const double upper = std::pow (10.0, releasing / 20.0);
        const double lower = std::pow (10.0, -engaging / 20.0);

        for (int k = 0; k < numBins_; ++k)
        {
            const size_t i = (size_t) k;
            const double factor = std::max (1.0e-30, 1.0 + transfer[k]);

            if (! primed)
            {
                previousFactor[i] = factor;
                continue;
            }

            const double previous = previousFactor[i];

            if (factor <= previous * upper && factor >= previous * lower)
            {
                previousFactor[i] = factor;
                transfer[k] = factor - 1.0;
                continue;
            }

            const double limited = previous * (factor > previous ? upper : lower);

            previousFactor[i] = limited;
            transfer[k] = limited - 1.0;
        }

        primed = true;
    }

private:
    GainSlewParams params;

    int numBins_ = 0;
    double modulationControl = 0.0;
    bool primed = false;

    // The applied factor, 1 + F, as of the last frame. Linear rather than dB: see process().
    std::vector<double> previousFactor;
};

} // namespace lime
