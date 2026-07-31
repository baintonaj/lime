/*
  ==============================================================================

    Double-precision real FFT.

    juce::dsp::FFT is float-only (perform() takes Complex<float>*), so Lime
    cannot use it anywhere: every transform in this project stays in 64-bit.

    Two backends sit behind one interface. Both are always available on Apple
    platforms so that the two can be cross-checked against each other, which is
    one of the project's acceptance tests.

  ==============================================================================
*/

#pragma once

#include <complex>
#include <memory>
#include <vector>

#ifndef LIME_HAS_ACCELERATE
 #if defined (__APPLE__)
  #define LIME_HAS_ACCELERATE 1
 #else
  #define LIME_HAS_ACCELERATE 0
 #endif
#endif

namespace lime
{

enum class FftBackend
{
    automatic,   // Accelerate where available, portable otherwise
    accelerate,
    portable
};

/** Real-input FFT in double precision.

    Normalisation is fixed so that a forward followed by an inverse is the
    identity: forward() produces the unnormalised DFT

        X[k] = sum_n x[n] * exp (-2*pi*i*k*n/N),   k = 0 .. N/2

    and inverse() applies the 1/N itself. Bins 0 and N/2 of a real signal are
    purely real; inverse() ignores any imaginary part supplied in them rather
    than letting it corrupt the result.

    Not thread-safe: each instance owns scratch buffers, so give every thread
    (and every concurrent STFT path) its own.
*/
class Fft64
{
public:
    /** @param orderIn  transform size is 1 << orderIn; must be >= 2. */
    explicit Fft64 (int orderIn, FftBackend backend = FftBackend::automatic);
    ~Fft64();

    Fft64 (Fft64&&) noexcept;
    Fft64& operator= (Fft64&&) noexcept;

    Fft64 (const Fft64&) = delete;
    Fft64& operator= (const Fft64&) = delete;

    int getOrder()   const noexcept { return order; }
    int getSize()    const noexcept { return size; }
    int getNumBins() const noexcept { return size / 2 + 1; }

    FftBackend getBackend() const noexcept { return activeBackend; }

    /** N real samples in, N/2+1 complex bins out. */
    void forward (const double* timeIn, std::complex<double>* binsOut) noexcept;

    /** N/2+1 complex bins in, N real samples out, 1/N already applied. */
    void inverse (const std::complex<double>* binsIn, double* timeOut) noexcept;

    /** True if this build can offer the Accelerate backend at all. */
    static constexpr bool accelerateAvailable() noexcept { return LIME_HAS_ACCELERATE != 0; }

private:
    void forwardPortable (const double*, std::complex<double>*) noexcept;
    void inversePortable (const std::complex<double>*, double*) noexcept;

    /** In-place iterative radix-2 complex FFT over `half` points. */
    void complexTransform (std::complex<double>* data, bool doInverse) noexcept;

    int order = 0;
    int size = 0;
    int half = 0;               // size / 2, the length of the packed complex transform
    FftBackend activeBackend = FftBackend::portable;

    // Portable backend state.
    std::vector<std::complex<double>> packed;      // half-length working buffer
    std::vector<std::complex<double>> twiddles;    // exp(-2*pi*i*k/half), k < half/2
    std::vector<std::complex<double>> untangle;    // exp(-2*pi*i*k/size), k <= half
    std::vector<int> bitReverse;

    struct Accelerated;
    std::unique_ptr<Accelerated> accel;
};

} // namespace lime
