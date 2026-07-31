/*
  ==============================================================================

    PhaseRotator — see the header for what it is and why it is bypassed at zero.

  ==============================================================================
*/

#include "PhaseRotator.h"

#include <algorithm>
#include <cmath>

namespace lime
{

namespace
{
    /** Niemitalo's 90 degree phase-difference network.

        Each value is the `a` of a second-order allpass in z^-2,

            H(z) = (a^2 - z^-2) / (1 - a^2 z^-2)

        so the section's coefficient is a^2. The two branches share a magnitude response of
        exactly one and differ in phase by 90 degrees across the band; the in-phase branch
        carries an extra sample of delay, which is part of the design rather than an
        accident of implementation. */
    constexpr double inPhasePoles[PhaseRotator::sectionsPerBranch] =
    {
        0.6923877778065, 0.9360654322959, 0.9882295226860, 0.9987488452737
    };

    constexpr double quadraturePoles[PhaseRotator::sectionsPerBranch] =
    {
        0.4021921162426, 0.8561710882420, 0.9722909545651, 0.9952884791278
    };

    std::vector<BiquadCoeffs> branchCoeffs (const double (&poles)[PhaseRotator::sectionsPerBranch])
    {
        std::vector<BiquadCoeffs> out;
        out.reserve (PhaseRotator::sectionsPerBranch);

        for (double pole : poles)
        {
            const double c = pole * pole;

            // (c - z^-2) / (1 - c z^-2), in the b/a convention BiquadCoeffs documents.
            BiquadCoeffs section;
            section.b0 = c;
            section.b1 = 0.0;
            section.b2 = -1.0;
            section.a1 = 0.0;
            section.a2 = -c;

            out.push_back (section);
        }

        return out;
    }

    constexpr double pi = 3.14159265358979323846;
}

//==============================================================================
void PhaseRotator::prepare (double sampleRate, int numChannels)
{
    const auto iCoeffs = branchCoeffs (inPhasePoles);
    const auto qCoeffs = branchCoeffs (quadraturePoles);

    channels.clear();
    channels.resize ((size_t) std::max (0, numChannels));

    for (auto& c : channels)
    {
        c.inPhase.setCoeffs (iCoeffs);
        c.quadrature.setCoeffs (qCoeffs);
        c.delayed = 0.0;
    }

    // Snapped, not faded: prepare is not a control movement.
    angle.prepare (sampleRate, defaultSettlingSeconds, targetDegrees * pi / 180.0);
    engaged.prepare (sampleRate, defaultSettlingSeconds, targetDegrees != 0.0 ? 1.0 : 0.0);
}

void PhaseRotator::reset()
{
    for (auto& c : channels)
    {
        c.inPhase.reset();
        c.quadrature.reset();
        c.delayed = 0.0;
    }
}

void PhaseRotator::setPhaseDegrees (double degrees) noexcept
{
    targetDegrees = std::clamp (degrees, -180.0, 180.0);

    angle.setTarget (targetDegrees * pi / 180.0);
    engaged.setTarget (targetDegrees != 0.0 ? 1.0 : 0.0);
}

void PhaseRotator::snapPhaseDegrees (double degrees) noexcept
{
    targetDegrees = std::clamp (degrees, -180.0, 180.0);

    angle.snapTo (targetDegrees * pi / 180.0);
    engaged.snapTo (targetDegrees != 0.0 ? 1.0 : 0.0);
}

bool PhaseRotator::isEngaged() const noexcept
{
    return engaged.isMoving() || engaged.getCurrent() > 1.0e-6;
}

void PhaseRotator::process (double* const* channelData, int numChannels, int numSamples)
{
    if (! isEngaged())
        return;

    const int usable = std::min (numChannels, (int) channels.size());

    // The angle spends nearly all of its life settled — it only moves for 75 ms after the
    // control does — so the cos/sin pair is worth computing once per block and holding.
    // While the smoother is actually mid-ramp the pair has a new value every sample, and
    // per-sample transcendentals for the length of a knob gesture are a fair price.
    const bool angleMoving = angle.isMoving();

    double cosTheta = 0.0, sinTheta = 0.0;

    if (! angleMoving)
    {
        const double theta = angle.getCurrent();
        cosTheta = std::cos (theta);
        sinTheta = std::sin (theta);
    }

    for (int n = 0; n < numSamples; ++n)
    {
        if (angleMoving)
        {
            const double theta = angle.next();
            cosTheta = std::cos (theta);
            sinTheta = std::sin (theta);
        }

        const double amount = engaged.next();

        for (int ch = 0; ch < usable; ++ch)
        {
            auto& c = channels[(size_t) ch];
            const double x = channelData[ch][n];

            // The in-phase branch runs a sample behind, which is what puts the two branches
            // 90 degrees apart rather than 90 degrees apart plus a sample of skew.
            const double i = c.delayed;
            c.delayed = c.inPhase.processSample (x);

            const double q = c.quadrature.processSample (x);

            const double rotated = i * cosTheta - q * sinTheta;

            // Faded against the untouched signal, so that leaving zero engages the network
            // rather than switching it in. Below one the two paths differ in phase, so this
            // is a brief comb rather than a clean crossfade — which is the honest cost of
            // being exactly transparent at zero, and it lasts 75 ms.
            channelData[ch][n] = x + amount * (rotated - x);
        }
    }
}

} // namespace lime
