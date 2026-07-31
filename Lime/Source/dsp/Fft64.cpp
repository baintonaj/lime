/*
  ==============================================================================

    Fft64 — double-precision real FFT, Accelerate and portable backends.

  ==============================================================================
*/

#include "Fft64.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#if LIME_HAS_ACCELERATE
 #include <Accelerate/Accelerate.h>
#endif

namespace lime
{

namespace
{
    constexpr double twoPi = 6.283185307179586476925286766559;

    int reverseBits (int value, int numBits) noexcept
    {
        int result = 0;

        for (int i = 0; i < numBits; ++i)
        {
            result = (result << 1) | (value & 1);
            value >>= 1;
        }

        return result;
    }
}

//==============================================================================
/** Accelerate state, kept out of the header so vDSP types do not leak.

    vDSP's real transform packs N/2 complex values with DC in realp[0] and
    Nyquist in imagp[0], and its forward result is scaled by 2 relative to the
    mathematical DFT. Both are undone here so the class's contract holds
    regardless of backend.
*/
struct Fft64::Accelerated
{
   #if LIME_HAS_ACCELERATE
    Accelerated (int order, int halfSize)
        : real ((size_t) halfSize), imag ((size_t) halfSize)
    {
        setup = vDSP_create_fftsetupD ((vDSP_Length) order, kFFTRadix2);
        split.realp = real.data();
        split.imagp = imag.data();
    }

    ~Accelerated()
    {
        if (setup != nullptr)
            vDSP_destroy_fftsetupD (setup);
    }

    bool isValid() const noexcept { return setup != nullptr; }

    FFTSetupD setup = nullptr;
    DSPDoubleSplitComplex split {};
    std::vector<double> real, imag;
   #else
    Accelerated (int, int) {}
    bool isValid() const noexcept { return false; }
   #endif
};

//==============================================================================
Fft64::Fft64 (int orderIn, FftBackend backend)
    : order (orderIn),
      size (1 << orderIn),
      half (1 << (orderIn - 1))
{
    assert (orderIn >= 2 && "Fft64 needs at least 4 points");

    const bool wantAccelerate = accelerateAvailable()
                             && (backend == FftBackend::automatic || backend == FftBackend::accelerate);

    if (wantAccelerate)
    {
        auto candidate = std::make_unique<Accelerated> (order, half);

        if (candidate->isValid())
        {
            accel = std::move (candidate);
            activeBackend = FftBackend::accelerate;
        }
    }

    if (activeBackend != FftBackend::accelerate)
    {
        activeBackend = FftBackend::portable;

        // The tables serve only the portable transform, and the backend is fixed
        // per instance, so an instance that landed on Accelerate never touches
        // them. Cross-validation of the two backends constructs an explicitly
        // portable instance of its own rather than reaching into this one.
        packed.resize ((size_t) half);

        twiddles.resize ((size_t) std::max (1, half / 2));
        for (size_t k = 0; k < twiddles.size(); ++k)
        {
            const double theta = -twoPi * (double) k / (double) half;
            twiddles[k] = { std::cos (theta), std::sin (theta) };
        }

        untangle.resize ((size_t) half + 1);
        for (size_t k = 0; k < untangle.size(); ++k)
        {
            const double theta = -twoPi * (double) k / (double) size;
            untangle[k] = { std::cos (theta), std::sin (theta) };
        }

        const int halfBits = order - 1;
        bitReverse.resize ((size_t) half);
        for (int i = 0; i < half; ++i)
            bitReverse[(size_t) i] = reverseBits (i, halfBits);
    }
}

Fft64::~Fft64() = default;
Fft64::Fft64 (Fft64&&) noexcept = default;
Fft64& Fft64::operator= (Fft64&&) noexcept = default;

//==============================================================================
void Fft64::forward (const double* timeIn, std::complex<double>* binsOut) noexcept
{
   #if LIME_HAS_ACCELERATE
    if (activeBackend == FftBackend::accelerate)
    {
        vDSP_ctozD (reinterpret_cast<const DSPDoubleComplex*> (timeIn), 2,
                    &accel->split, 1, (vDSP_Length) half);

        vDSP_fft_zripD (accel->setup, &accel->split, 1, (vDSP_Length) order, kFFTDirection_Forward);

        // Undo vDSP's factor of 2 and unpack DC / Nyquist out of bin 0.
        constexpr double unscale = 0.5;

        binsOut[0]    = { accel->split.realp[0] * unscale, 0.0 };
        binsOut[half] = { accel->split.imagp[0] * unscale, 0.0 };

        for (int k = 1; k < half; ++k)
            binsOut[k] = { accel->split.realp[k] * unscale,
                           accel->split.imagp[k] * unscale };

        return;
    }
   #endif

    forwardPortable (timeIn, binsOut);
}

void Fft64::inverse (const std::complex<double>* binsIn, double* timeOut) noexcept
{
   #if LIME_HAS_ACCELERATE
    if (activeBackend == FftBackend::accelerate)
    {
        // Re-apply vDSP's factor of 2 and fold Nyquist back into bin 0.
        accel->split.realp[0] = binsIn[0].real() * 2.0;
        accel->split.imagp[0] = binsIn[half].real() * 2.0;

        for (int k = 1; k < half; ++k)
        {
            accel->split.realp[k] = binsIn[k].real() * 2.0;
            accel->split.imagp[k] = binsIn[k].imag() * 2.0;
        }

        vDSP_fft_zripD (accel->setup, &accel->split, 1, (vDSP_Length) order, kFFTDirection_Inverse);

        vDSP_ztocD (&accel->split, 1,
                    reinterpret_cast<DSPDoubleComplex*> (timeOut), 2, (vDSP_Length) half);

        // vDSP's real round trip scales by 2N, and we already put back a 2.
        const double scale = 1.0 / (2.0 * (double) size);
        vDSP_vsmulD (timeOut, 1, &scale, timeOut, 1, (vDSP_Length) size);

        return;
    }
   #endif

    inversePortable (binsIn, timeOut);
}

//==============================================================================
void Fft64::forwardPortable (const double* timeIn, std::complex<double>* binsOut) noexcept
{
    // Pack the real input as half as many complex values: even samples into the
    // real part, odd samples into the imaginary part.
    for (int n = 0; n < half; ++n)
        packed[(size_t) n] = { timeIn[2 * n], timeIn[2 * n + 1] };

    complexTransform (packed.data(), false);

    // Untangle. With Z the half-length transform of the packed signal,
    //   Xeven[k] = (Z[k] + conj Z[half-k]) / 2
    //   Xodd [k] = (Z[k] - conj Z[half-k]) / 2i
    //   X[k]     = Xeven[k] + exp(-2*pi*i*k/N) * Xodd[k]
    constexpr std::complex<double> overTwoI { 0.0, -0.5 };   // 1 / 2i

    for (int k = 0; k <= half; ++k)
    {
        const auto zk  = packed[(size_t) (k % half)];
        const auto zmk = std::conj (packed[(size_t) ((half - k) % half)]);

        const auto even = (zk + zmk) * 0.5;
        const auto odd  = (zk - zmk) * overTwoI;

        binsOut[k] = even + untangle[(size_t) k] * odd;
    }

    // Bins 0 and N/2 of a real signal are real by construction; strip the
    // rounding residue so downstream code can rely on it.
    binsOut[0]    = { binsOut[0].real(),    0.0 };
    binsOut[half] = { binsOut[half].real(), 0.0 };
}

void Fft64::inversePortable (const std::complex<double>* binsIn, double* timeOut) noexcept
{
    // Rebuild the packed half-length spectrum, inverting the untangle step.
    // X[k+half] = conj X[half-k] for a real signal, so
    //   Xeven[k] = (X[k] + conj X[half-k]) / 2
    //   Xodd [k] = (X[k] - conj X[half-k]) / 2 * conj(exp(-2*pi*i*k/N))
    //   Z[k]     = Xeven[k] + i * Xodd[k]
    constexpr std::complex<double> i { 0.0, 1.0 };

    for (int k = 0; k < half; ++k)
    {
        auto a = binsIn[(size_t) k];
        auto b = std::conj (binsIn[(size_t) (half - k)]);

        if (k == 0)
        {
            // Bin 0 and bin half are the real DC and Nyquist terms.
            a = { binsIn[0].real(), 0.0 };
            b = { binsIn[(size_t) half].real(), 0.0 };
        }

        const auto even = (a + b) * 0.5;
        const auto odd  = (a - b) * 0.5 * std::conj (untangle[(size_t) k]);

        packed[(size_t) k] = even + i * odd;
    }

    complexTransform (packed.data(), true);   // includes the 1/half

    for (int n = 0; n < half; ++n)
    {
        timeOut[2 * n]     = packed[(size_t) n].real();
        timeOut[2 * n + 1] = packed[(size_t) n].imag();
    }
}

//==============================================================================
void Fft64::complexTransform (std::complex<double>* data, bool doInverse) noexcept
{
    for (int idx = 0; idx < half; ++idx)
    {
        const int rev = bitReverse[(size_t) idx];

        if (idx < rev)
            std::swap (data[idx], data[rev]);
    }

    for (int len = 2; len <= half; len <<= 1)
    {
        const int halfLen = len >> 1;
        const int step = half / len;

        for (int base = 0; base < half; base += len)
        {
            for (int j = 0; j < halfLen; ++j)
            {
                auto w = twiddles[(size_t) (j * step)];

                if (doInverse)
                    w = std::conj (w);

                const auto u = data[base + j];
                const auto v = data[base + j + halfLen] * w;

                data[base + j]           = u + v;
                data[base + j + halfLen] = u - v;
            }
        }
    }

    if (doInverse)
    {
        const double scale = 1.0 / (double) half;

        for (int idx = 0; idx < half; ++idx)
            data[idx] *= scale;
    }
}

} // namespace lime
