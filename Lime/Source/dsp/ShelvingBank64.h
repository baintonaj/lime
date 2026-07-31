/*
  ==============================================================================

    ShelvingBank64 — a Bark-spaced bank of peaking biquads that realises a spectral
    gain profile in the time domain, and inverts it exactly.

    This exists because applying gain inside the STFT cannot be undone. A modified
    spectrum with abrupt bin-to-bin gain changes is not the transform of any signal,
    so no decoder can recover the input: measured, a single-bin 24 dB notch survives
    an encode/decode round trip only to -16 dB, against -126 dB for a gradual profile.

    A biquad escapes that entirely. From

        y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]

    the input is recovered as

        x[n] = (y[n] + a1 y[n-1] + a2 y[n-2] - b1 x[n-1] - b2 x[n-2]) / b0

    which is exact sample by sample and stays exact even while the coefficients
    change, provided both directions use the same coefficients at the same sample and
    b0 is non-zero. Time variation, which defeats the STFT, costs nothing here.

    The STFT is still what decides the gains — it is the right tool for analysis. Only
    the application moves. That is also how the analog works: control derived from
    analysis, gain applied by exactly invertible networks.

    Overlapping peaking filters do not act independently, so a requested profile is not
    obtained by simply setting each band to the value wanted at its centre. The
    interaction matrix is built and inverted once in prepare(), and setTargetDb() then
    solves for the band gains that produce the requested profile.

  ==============================================================================
*/

#pragma once

#include "Bark.h"
#include "Biquad64.h"

#include <algorithm>
#include <complex>
#include <cmath>
#include <vector>

namespace lime
{

class ShelvingBank64
{
public:
    /** @param bandsPerBark  density of the bank; 1.0 gives one filter per Bark. */
    void prepare (double sampleRate, double bandsPerBark = 1.0, double q = 2.5)
    {
        currentSampleRate = sampleRate;
        filterQ = q;

        // Centres spaced evenly in Bark from 40 Hz to just below Nyquist, which is
        // where the process has anything to say.
        const double lowBark = barkFromHz (40.0);
        const double highBark = barkFromHz (0.45 * sampleRate);
        const int count = std::max (2, (int) std::round ((highBark - lowBark) * bandsPerBark));

        centreHz.assign ((size_t) count, 0.0);

        for (int i = 0; i < count; ++i)
        {
            const double t = (double) i / (double) (count - 1);
            centreHz[(size_t) i] = hzFromBark (lowBark + t * (highBark - lowBark));
        }

        sections.assign ((size_t) count, {});
        solved.assign ((size_t) count, 0.0);

        buildInteractionMatrix();

        // Start flat; solving the all-zero profile also writes every section's
        // coefficients, so the bank is ready to process immediately.
        const std::vector<double> flat ((size_t) count, 0.0);
        setTargetDb (flat.data());
    }

    void reset()
    {
        for (auto& s : sections)
        {
            s.forward.reset();
            s.inverse.reset();
        }
    }

    int getNumBands() const noexcept { return (int) centreHz.size(); }
    double getCentreHz (int band) const noexcept { return centreHz[(size_t) band]; }

    /** Requests a gain profile sampled at the band centres, in dB. Solves for the
        band gains that produce it once the bands' overlap is accounted for. */
    void setTargetDb (const double* targetDb)
    {
        const int n = getNumBands();

        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;

            for (int j = 0; j < n; ++j)
                sum += inverseMatrix[(size_t) (i * n + j)] * targetDb[j];

            // Keep the solution sane: an ill-conditioned request must not produce
            // wild band gains that ring or make b0 approach zero.
            solved[(size_t) i] = std::clamp (sum, -36.0, 36.0);
        }

        for (int i = 0; i < n; ++i)
            updateSection (i, solved[(size_t) i]);
    }

    /** Gain in dB actually assigned to a band, after the solve. */
    double getSolvedGainDb (int band) const noexcept { return solved[(size_t) band]; }

    /** Composite response of the bank in dB at a frequency. */
    double responseDbAt (double hz) const
    {
        const double omega = 6.283185307179586476925286766559 * hz / currentSampleRate;
        std::complex<double> h { 1.0, 0.0 };

        for (const auto& s : sections)
            h *= s.forward.getCoeffs().responseAt (omega);

        return 20.0 * std::log10 (std::abs (h));
    }

    void process (double* samples, int numSamples) noexcept
    {
        for (auto& s : sections)
            for (int n = 0; n < numSamples; ++n)
                samples[n] = s.forward.processSample (samples[n]);
    }

    /** Exact inverse, applied in reverse order. */
    void processInverse (double* samples, int numSamples) noexcept
    {
        for (auto it = sections.rbegin(); it != sections.rend(); ++it)
            for (int n = 0; n < numSamples; ++n)
                samples[n] = it->inverse.processInverseSample (samples[n]);
    }

private:
    struct Section
    {
        // Direct form I, because its state is the real input/output history and so the
        // pair stays exact while the coefficients move. See BiquadDirectFormI64.
        BiquadDirectFormI64 forward, inverse;
    };

    /** Cookbook peaking filter, plus its exact inverse. */
    void updateSection (int index, double gainDb)
    {
        constexpr double twoPi = 6.283185307179586476925286766559;

        const double a = std::pow (10.0, gainDb / 40.0);
        const double omega = twoPi * centreHz[(size_t) index] / currentSampleRate;
        const double alpha = std::sin (omega) / (2.0 * filterQ);
        const double cosine = std::cos (omega);

        BiquadCoeffs c;
        const double a0 = 1.0 + alpha / a;

        c.b0 = (1.0 + alpha * a) / a0;
        c.b1 = (-2.0 * cosine) / a0;
        c.b2 = (1.0 - alpha * a) / a0;
        c.a1 = (-2.0 * cosine) / a0;
        c.a2 = (1.0 - alpha / a) / a0;

        // Both hold the *same* coefficients: the inverse recurrence is applied by
        // processInverseSample rather than by swapping numerator and denominator, so
        // that the two reference identical past samples.
        sections[(size_t) index].forward.setCoeffs (c);
        sections[(size_t) index].inverse.setCoeffs (c);
    }

    /** M[i][j] is the dB contribution of band j, at unit gain, to the response at
        band i's centre. Inverted so a target profile can be solved for directly.

        This treats a band's off-centre dB response as proportional to its centre
        gain, probed at 6 dB. For a cookbook peaking section that is only an
        approximation — the skirt's shape itself moves with gain, since the
        bandwidth is held at the half-gain points — so the solve is exact at the
        probe depth and drifts as the assigned gains deepen towards the ±36 dB
        clamp: measured, about 2 dB of profile error on a 24 dB tilt. Acceptable
        for a control profile; anything needing the realised response should read
        responseDbAt() after the solve rather than trusting the request. */
    void buildInteractionMatrix()
    {
        constexpr double probeDb = 6.0;
        const int n = getNumBands();

        std::vector<double> matrix ((size_t) (n * n), 0.0);

        for (int j = 0; j < n; ++j)
        {
            updateSection (j, probeDb);

            for (int i = 0; i < n; ++i)
                matrix[(size_t) (i * n + j)] = responseAtCentreOfOneSection (j, i) / probeDb;

            updateSection (j, 0.0);
        }

        invert (matrix, n);
        inverseMatrix = std::move (matrix);
    }

    double responseAtCentreOfOneSection (int section, int atCentre) const
    {
        constexpr double twoPi = 6.283185307179586476925286766559;
        const double omega = twoPi * centreHz[(size_t) atCentre] / currentSampleRate;

        return 20.0 * std::log10 (std::abs (sections[(size_t) section].forward.getCoeffs().responseAt (omega)));
    }

    /** Gauss-Jordan with partial pivoting, in place. Falls back to the identity if the
        matrix is singular, which leaves the bank acting per-centre rather than solved
        — degraded but never wrong. */
    static void invert (std::vector<double>& m, int n)
    {
        std::vector<double> inv ((size_t) (n * n), 0.0);

        for (int i = 0; i < n; ++i)
            inv[(size_t) (i * n + i)] = 1.0;

        for (int col = 0; col < n; ++col)
        {
            int pivot = col;

            for (int row = col + 1; row < n; ++row)
                if (std::abs (m[(size_t) (row * n + col)]) > std::abs (m[(size_t) (pivot * n + col)]))
                    pivot = row;

            if (std::abs (m[(size_t) (pivot * n + col)]) < 1.0e-12)
            {
                std::fill (m.begin(), m.end(), 0.0);

                for (int i = 0; i < n; ++i)
                    m[(size_t) (i * n + i)] = 1.0;

                return;
            }

            if (pivot != col)
                for (int k = 0; k < n; ++k)
                {
                    std::swap (m[(size_t) (col * n + k)], m[(size_t) (pivot * n + k)]);
                    std::swap (inv[(size_t) (col * n + k)], inv[(size_t) (pivot * n + k)]);
                }

            const double diagonal = m[(size_t) (col * n + col)];

            for (int k = 0; k < n; ++k)
            {
                m[(size_t) (col * n + k)] /= diagonal;
                inv[(size_t) (col * n + k)] /= diagonal;
            }

            for (int row = 0; row < n; ++row)
            {
                if (row == col)
                    continue;

                const double factor = m[(size_t) (row * n + col)];

                if (factor == 0.0)
                    continue;

                for (int k = 0; k < n; ++k)
                {
                    m[(size_t) (row * n + k)] -= factor * m[(size_t) (col * n + k)];
                    inv[(size_t) (row * n + k)] -= factor * inv[(size_t) (col * n + k)];
                }
            }
        }

        m = std::move (inv);
    }

    double currentSampleRate = 48000.0;
    double filterQ = 2.5;

    std::vector<double> centreHz, solved, inverseMatrix;
    std::vector<Section> sections;
};

} // namespace lime
