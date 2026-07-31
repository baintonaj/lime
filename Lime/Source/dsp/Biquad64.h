/*
  ==============================================================================

    Biquad64 — double-precision biquad section and cascade.

    Used for the fixed main-path networks (spectral skewing and antisaturation).
    Two properties matter here and are why the main path is not an STFT
    multiplier: a second-order analog shelf maps onto a biquad with no
    approximation in structure, and a minimum-phase biquad is inverted exactly by
    swapping its numerator and denominator. That makes the decoder's de-skewing
    network exact rather than an approximation, which frequency-sampled FIR
    design cannot achieve because the inverse response is largely anti-causal.
    OverlapSave64 keeps the FIR alternative alive as a tested reference to
    compare against.

    Transposed direct form II: good numerical behaviour in double precision and
    only two state words per section.

  ==============================================================================
*/

#pragma once

#include <cmath>
#include <complex>
#include <vector>

namespace lime
{

struct BiquadCoeffs
{
    // H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2)
    double b0 = 1.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;

    /** Exact inverse. Valid when the section is minimum phase, which every
        network in this project is; b0 must be non-zero. */
    BiquadCoeffs inverted() const noexcept
    {
        BiquadCoeffs out;
        const double norm = 1.0 / b0;

        out.b0 = 1.0 * norm;
        out.b1 = a1 * norm;
        out.b2 = a2 * norm;
        out.a1 = b1 * norm;
        out.a2 = b2 * norm;

        return out;
    }

    /** Frequency response at a normalised angular frequency (radians/sample). */
    std::complex<double> responseAt (double omega) const noexcept
    {
        const std::complex<double> z1 { std::cos (-omega), std::sin (-omega) };
        const std::complex<double> z2 = z1 * z1;

        return (b0 + b1 * z1 + b2 * z2) / (1.0 + a1 * z1 + a2 * z2);
    }
};

//==============================================================================
class Biquad64
{
public:
    void setCoeffs (const BiquadCoeffs& c) noexcept { coeffs = c; }
    const BiquadCoeffs& getCoeffs() const noexcept  { return coeffs; }

    void reset() noexcept { s1 = s2 = 0.0; }

    double processSample (double x) noexcept
    {
        const double y = coeffs.b0 * x + s1;
        s1 = coeffs.b1 * x - coeffs.a1 * y + s2;
        s2 = coeffs.b2 * x - coeffs.a2 * y;
        return y;
    }

    /** What the next processSample() would return for a zero input, without
        advancing the state. The output is b0 * x + s1, so as a function of the
        current sample it is affine with slope b0 and this value as its offset —
        which is exactly what a solver inverting the section's sample map needs
        to read. */
    double outputForZeroInput() const noexcept { return s1; }

private:
    BiquadCoeffs coeffs;
    double s1 = 0.0, s2 = 0.0;
};

//==============================================================================
/** Direct form I biquad, kept alongside the transposed form for one specific reason.

    Its state is the actual history of inputs and outputs, not filter-internal sums.
    That is what makes it exactly invertible while its coefficients are *changing*:
    from y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2] the inverse
    reads x[n] = (y[n] + a1 y[n-1] + a2 y[n-2] - b1 x[n-1] - b2 x[n-2]) / b0, and both
    directions reference the same physical past samples. A transposed direct form II
    pair cannot do this: its internal states are in the units of whichever coefficients
    were in force when they were written, so once the coefficients move the forward and
    inverse states stop corresponding and the pair diverges. Measured, with gains
    changing every block, transposed form II diverged to +46 dB while this form holds.

    The cost is four state words instead of two and slightly worse numerical behaviour
    with poles close to the unit circle, neither of which matters for a gently peaking
    filter.
*/
class BiquadDirectFormI64
{
public:
    void setCoeffs (const BiquadCoeffs& c) noexcept { coeffs = c; }
    const BiquadCoeffs& getCoeffs() const noexcept  { return coeffs; }

    void reset() noexcept { x1 = x2 = y1 = y2 = 0.0; }

    double processSample (double x) noexcept
    {
        const double y = coeffs.b0 * x + coeffs.b1 * x1 + coeffs.b2 * x2
                       - coeffs.a1 * y1 - coeffs.a2 * y2;

        x2 = x1; x1 = x;
        y2 = y1; y1 = y;

        return y;
    }

    /** Recovers the input that produced `y`, using the same coefficients the forward
        filter used for that sample. */
    double processInverseSample (double y) noexcept
    {
        // Here x1/x2 hold past *inputs to the forward filter*, which are this
        // function's outputs, and y1/y2 hold past forward outputs, its inputs.
        const double x = (y + coeffs.a1 * y1 + coeffs.a2 * y2
                            - coeffs.b1 * x1 - coeffs.b2 * x2) / coeffs.b0;

        y2 = y1; y1 = y;
        x2 = x1; x1 = x;

        return x;
    }

private:
    BiquadCoeffs coeffs;
    double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
};

//==============================================================================
class BiquadCascade64
{
public:
    void setCoeffs (const std::vector<BiquadCoeffs>& list)
    {
        sections.resize (list.size());

        for (size_t i = 0; i < list.size(); ++i)
            sections[i].setCoeffs (list[i]);
    }

    /** Replaces the cascade with the exact inverse of `list`. Section order is
        irrelevant mathematically, but is reversed so a forward-then-inverse pair
        reads symmetrically. */
    void setInverseCoeffs (const std::vector<BiquadCoeffs>& list)
    {
        sections.resize (list.size());

        for (size_t i = 0; i < list.size(); ++i)
            sections[i].setCoeffs (list[list.size() - 1 - i].inverted());
    }

    void reset() noexcept
    {
        for (auto& s : sections)
            s.reset();
    }

    size_t getNumSections() const noexcept { return sections.size(); }

    double processSample (double x) noexcept
    {
        for (auto& s : sections)
            x = s.processSample (x);

        return x;
    }

    void process (double* samples, int numSamples) noexcept
    {
        for (auto& s : sections)
            for (int n = 0; n < numSamples; ++n)
                samples[n] = s.processSample (samples[n]);
    }

    /** Combined response of the cascade at a normalised angular frequency. */
    std::complex<double> responseAt (double omega) const noexcept
    {
        std::complex<double> h { 1.0, 0.0 };

        for (const auto& s : sections)
            h *= s.getCoeffs().responseAt (omega);

        return h;
    }

private:
    std::vector<Biquad64> sections;
};

} // namespace lime
