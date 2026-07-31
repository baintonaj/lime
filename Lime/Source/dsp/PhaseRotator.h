/*
  ==============================================================================

    PhaseRotator — a rotation control, -180 to +180 degrees, flat in magnitude.

    What is constant across frequency is the *rotation*: moving the control from one
    setting to another shifts the phase at every frequency by the same number of degrees,
    the difference between the two settings. The absolute phase against the input is
    another matter — the network's branches are allpasses, so everything that comes out
    carries their own frequency-dependent shift. That shift is the same at every setting
    of the knob, which is why it cancels out of any comparison between settings, and why
    this is a rotation control rather than a phase-alignment one. This is what "phase
    rotation" means on a broadcast processor or a mastering desk, and it is what the
    tests measure: magnitude untouched everywhere, and the difference between any two
    settings a constant rotation at every frequency.

    It works by building a quadrature pair. Two allpass cascades, fed the same signal, are
    designed so their outputs differ by 90 degrees everywhere in the band: call them I and
    Q. Any phase in between is then a plain sum,

        y = I cos(theta) - Q sin(theta)

    which rotates the signal by theta while leaving |y| alone, because I and Q are already
    the two axes of the same vector. At either end of the range this is exactly -I — a
    polarity inversion of what the network passes at zero rotation, not of the input, since
    the I branch carries its own allpass shift; the control passes through it rather than
    stopping at it.

    The range runs either side of zero rather than from zero upward because a rotation has a
    direction, and which way round a blend is leaning is the thing being listened for. A
    control that put no rotation at one end would spend half its travel getting back to where
    it started, and would make +30 and -30 degrees — which sound different — hard to compare.

    The cascades are Niemitalo's four-section-per-branch design, each section a second-order
    allpass in z^-2 with a single real coefficient, plus a one-sample delay on the branch
    that needs it. Its 90 degree difference holds to a fraction of a degree from about
    0.005 to 0.495 of the sample rate — the whole audio band at any rate this plugin runs
    at, and its error at the extremes is far outside anything the process is working on.

    Bypassed outright at zero.

    That matters more here than it would in most plugins. Lime's encode and decode are meant
    to be exact inverses across two instances, and an allpass network in circuit on one side
    and not the other would break that even at "no rotation" — the I branch is magnitude-flat
    but not phase-flat, so passing through it is not the same as not passing through it. At
    zero the signal therefore goes round the network entirely, and the default is zero.

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"
#include "ControlSmoother.h"

#include <vector>

namespace lime
{

class PhaseRotator
{
public:
    void prepare (double sampleRate, int numChannels);
    void reset();

    /** @param degrees  -180 to +180. Smoothed, so it can be swept while audio is running. */
    void setPhaseDegrees (double degrees) noexcept;

    /** Jumps, for prepareToPlay and for a restored session — saved settings must be in
        force on the first sample rather than faded in. */
    void snapPhaseDegrees (double degrees) noexcept;

    /** True while any rotation is in circuit, including while it is fading out. Callers can
        skip the whole network when this is false, and at zero it costs nothing. */
    bool isEngaged() const noexcept;

    /** Rotates in place. */
    void process (double* const* channelData, int numChannels, int numSamples);

    /** The design's own limits, exposed for the tests that measure them. */
    static constexpr int sectionsPerBranch = 4;

private:
    struct Channel
    {
        BiquadCascade64 inPhase, quadrature;
        double delayed = 0.0;   // the one-sample delay on the in-phase branch
    };

    std::vector<Channel> channels;

    // The rotation itself, and a separate fade for whether the network is in circuit at
    // all. Two smoothers rather than one because they are not the same question: the angle
    // moves continuously, while engaging is a step between two different signal paths.
    ControlSmoother angle, engaged;

    double targetDegrees = 0.0;
};

} // namespace lime
