/*
  ==============================================================================

    WetLowPass — the wet-only low pass, last in the processed path before the
    blend against dry.

    The process writes harmonics above the programme. Mixed in parallel, those
    are often the only part of the effect a user wants less of, so this filter
    shades them off without touching the dry path: it sits after everything the
    wet path does — trims and loudness makeup included — and immediately ahead
    of the mix arithmetic, which still sees the dry signal untouched.

    Third-order Butterworth: one real pole and a Q = 1 quadratic, both at the
    corner. 18 dB an octave is enough slope to take harmonics down decisively
    while the passband stays maximally flat, so the corner reads true — the
    knob's number is the -3 dB point, exactly, because the bilinear transform is
    prewarped there.

    The corner is smoothed here rather than in the plugin wrapper because the
    smoothed value has to drive coefficient redesign, which is a DSP concern and
    has to stay testable without JUCE — the same reason the engine owns its loop
    trims. The sections are direct form I: with coefficients moving under a
    sweep, transposed form II state is in the units of whichever coefficients
    wrote it (see the note in Biquad64.h), and DF-I's state — the physical
    input and output history — stays meaningful across every redesign.

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"
#include "ControlSmoother.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace lime
{

/** Third-order Butterworth low pass as two sections: a first-order pole stated
    as a degenerate biquad (b2 = a2 = 0), then a quadratic with Q = 1 — the
    Butterworth pole pair at 60 degrees. Bilinear with the corner prewarped so
    the -3.01 dB point lands on `cornerHz`; the corner is clamped just below
    Nyquist, as the band filters in AtypeEngine are, so any panel setting stays
    designable at any sample rate. */
inline std::array<BiquadCoeffs, 2> butterworthLowPass18 (double cornerHz, double sampleRate)
{
    constexpr double pi = 3.14159265358979323846;

    const double nyquist = 0.5 * sampleRate;
    const double clamped = std::clamp (cornerHz, 1.0e-4 * nyquist, 0.49 * sampleRate);
    const double k = std::tan (pi * clamped / sampleRate);
    const double kk = k * k;

    std::array<BiquadCoeffs, 2> out;

    {
        auto& c = out[0];
        const double a0 = 1.0 + k;

        c.b0 = k / a0;
        c.b1 = k / a0;
        c.b2 = 0.0;
        c.a1 = (k - 1.0) / a0;
        c.a2 = 0.0;
    }

    {
        auto& c = out[1];
        const double a0 = 1.0 + k + kk;

        c.b0 = kk / a0;
        c.b1 = 2.0 * kk / a0;
        c.b2 = kk / a0;
        c.a1 = 2.0 * (kk - 1.0) / a0;
        c.a2 = (1.0 - k + kk) / a0;
    }

    return out;
}

//==============================================================================
class WetLowPass
{
public:
    /** @param sampleRate       samples per second
        @param settlingSeconds  corner smoothing, seconds to 95% — see ControlSmoother
        @param numChannels      channels to allocate state for
        @param initialCornerHz  in force on the first sample, no ramp from anywhere
    */
    void prepare (double sampleRate_, double settlingSeconds, int numChannels, double initialCornerHz)
    {
        sampleRate = sampleRate_;
        corner.prepare (sampleRate, settlingSeconds, initialCornerHz);
        channels.assign ((size_t) std::max (0, numChannels), {});
        design (initialCornerHz);
    }

    /** The smoothed target, read from the parameter once per block. */
    void setCornerHz (double hz) noexcept        { corner.setTarget (hz); }

    /** Jumps, for a host loading a session: the saved corner is in force on the
        first sample instead of being swept to over the settling time. Section
        state is kept — a restore mid-stream should not click, and DF-I tolerates
        the coefficient step. */
    void snapCornerHz (double hz)
    {
        corner.snapTo (hz);
        design (hz);
    }

    /** Clears section state; coefficients and the smoothed corner stay put. */
    void reset() noexcept
    {
        for (auto& ch : channels)
        {
            ch.first.reset();
            ch.second.reset();
        }
    }

    bool isMoving() const noexcept               { return corner.isMoving(); }
    double getCornerHz() const noexcept          { return corner.getCurrent(); }

    void process (double* const* buffers, int numChannels, int numSamples) noexcept
    {
        const int limit = std::min (numChannels, (int) channels.size());

        for (int start = 0; start < numSamples;)
        {
            int piece = numSamples - start;

            // At rest the whole block runs on the coefficients already in
            // place — no tan(), which is the common case. Under a sweep the
            // block is walked in short pieces: the corner's trajectory advances
            // on the same per-sample clock as every other control in the
            // plugin, and the sections are redesigned once per piece rather
            // than once per sample, which would cost more than the filter.
            if (corner.isMoving())
            {
                piece = std::min (redesignQuantum, piece);

                for (int i = 0; i < piece; ++i)
                    corner.next();

                design (corner.getCurrent());
            }

            for (int ch = 0; ch < limit; ++ch)
            {
                auto& s = channels[(size_t) ch];
                double* samples = buffers[ch] + start;

                for (int n = 0; n < piece; ++n)
                    samples[n] = s.second.processSample (s.first.processSample (samples[n]));
            }

            start += piece;
        }
    }

private:
    void design (double cornerHz)
    {
        const auto coeffs = butterworthLowPass18 (cornerHz, sampleRate);

        for (auto& ch : channels)
        {
            ch.first.setCoeffs (coeffs[0]);
            ch.second.setCoeffs (coeffs[1]);
        }
    }

    // Short enough that one step of a 75 ms settle is far below audibility,
    // long enough that the tan() disappears into the per-sample work.
    static constexpr int redesignQuantum = 32;

    struct Sections
    {
        BiquadDirectFormI64 first, second;
    };

    ControlSmoother corner;
    std::vector<Sections> channels;
    double sampleRate = 48000.0;
};

} // namespace lime
