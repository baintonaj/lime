/*
  ==============================================================================

    OverlapSave64 — exact frequency-domain FIR convolution, double precision.

    A tested reference implementation, kept for the test suite and for
    comparison against the time-domain paths; nothing in the plugin's signal
    path uses it. The main path realises its fixed networks as biquad cascades
    (see Biquad64.h): biquads win when a low-order minimum-phase network must be
    inverted exactly, while overlap-save wins when a long or arbitrary FIR
    response must be applied verbatim — it is exact linear convolution, far
    cheaper than direct convolution because the block advance can be thousands
    of samples, and unlike windowed overlap-add it involves no window smearing
    a fixed per-bin multiplier into a circular-convolution approximation.

    Given an FFT of size N and a filter of M taps, each block transforms N
    samples of which the first M-1 are carried over from the previous block. The
    first M-1 output samples of the circular convolution are contaminated by
    wrap-around and discarded; the remaining L = N - M + 1 are exact. L may be
    set smaller than that maximum, which trades throughput for latency.

  ==============================================================================
*/

#pragma once

#include "Fft64.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstring>
#include <vector>

namespace lime
{

class OverlapSave64
{
public:
    /** @param fftOrder   transform size is 1 << fftOrder
        @param filterLen  number of taps, M; must satisfy M <= N
        @param advance    output samples produced per block, L; defaults to the
                          maximum N - M + 1. Latency equals L.
    */
    OverlapSave64 (int fftOrder, int filterLen, int advance = 0,
                   FftBackend backend = FftBackend::automatic)
        : fft (fftOrder, backend),
          size (1 << fftOrder),
          taps (filterLen),
          hop (advance > 0 ? advance : (1 << fftOrder) - filterLen + 1)
    {
        // The contract is isValid(): a filter longer than the transform leaves
        // no wrap-free output at all, so the default hop comes out non-positive,
        // and an oversized advance would read history that was never kept.
        // Assert it, then clamp so a release build handed bad geometry degrades
        // to legal buffer arithmetic instead of sizing buffers from a negative
        // count; valid geometry passes through untouched.
        assert (isValid());
        taps = std::min (std::max (taps, 1), size);
        hop = std::min (std::max (hop, 1), size - taps + 1);

        block.assign ((size_t) size, 0.0);
        scratch.assign ((size_t) size, 0.0);
        bins.assign ((size_t) fft.getNumBins(), { 0.0, 0.0 });
        response.assign ((size_t) fft.getNumBins(), { 1.0, 0.0 });

        inputPending.assign ((size_t) hop, 0.0);
        outputReady.assign ((size_t) hop, 0.0);
        reset();
    }

    void reset()
    {
        std::fill (block.begin(), block.end(), 0.0);
        std::fill (inputPending.begin(), inputPending.end(), 0.0);
        std::fill (outputReady.begin(), outputReady.end(), 0.0);
        pendingCount = 0;
        readyCount = hop;
    }

    int getSize()    const noexcept { return size; }
    int getNumTaps() const noexcept { return taps; }
    int getHop()     const noexcept { return hop; }
    int getNumBins() const noexcept { return fft.getNumBins(); }

    /** Latency in samples: exactly one hop, the input buffering and nothing
        else. A block's usable output is the linear convolution ending at the
        newest input sample — not, as in windowed overlap-add, a frame still
        awaiting later contributions — so a sample entering at time t has its
        convolved value computed as soon as its hop completes and emitted while
        the next hop is gathered, re-emerging at t + hop. The hop of silence
        primed in reset() covers that first interval. Confirmed by nulling
        against direct time-domain convolution. */
    int getLatencySamples() const noexcept { return hop; }

    bool isValid() const noexcept { return taps >= 1 && hop >= 1 && hop <= size - taps + 1; }

    /** Installs a filter from its impulse response. The response is zero-padded
        to the transform length and cached as a spectrum. */
    void setImpulseResponse (const double* ir, int numTaps)
    {
        taps = std::min (numTaps, size);

        std::fill (scratch.begin(), scratch.end(), 0.0);
        std::copy (ir, ir + taps, scratch.begin());

        fft.forward (scratch.data(), response.data());
    }

    /** Streams a block of any length through the filter. Input and output may
        alias. */
    void process (const double* input, double* output, int numSamples)
    {
        int done = 0;

        while (done < numSamples)
        {
            const int room = hop - pendingCount;
            const int avail = std::min (room, numSamples - done);

            std::memcpy (inputPending.data() + pendingCount, input + done, (size_t) avail * sizeof (double));

            const int emitFrom = hop - readyCount;
            std::memcpy (output + done, outputReady.data() + emitFrom, (size_t) avail * sizeof (double));

            pendingCount += avail;
            readyCount -= avail;
            done += avail;

            if (pendingCount == hop)
            {
                runBlock();
                pendingCount = 0;
                readyCount = hop;
            }
        }
    }

private:
    void runBlock()
    {
        // Slide in the new samples, keeping size - hop of history so that the
        // discarded head of the circular convolution covers the wrap-around.
        std::memmove (block.data(), block.data() + hop, (size_t) (size - hop) * sizeof (double));
        std::memcpy (block.data() + size - hop, inputPending.data(), (size_t) hop * sizeof (double));

        fft.forward (block.data(), bins.data());

        const int numBins = fft.getNumBins();

        for (int k = 0; k < numBins; ++k)
            bins[(size_t) k] *= response[(size_t) k];

        fft.inverse (bins.data(), scratch.data());

        // The last hop of the block is free of wrap-around contamination,
        // because the preceding size - hop >= taps - 1 samples are real history.
        std::memcpy (outputReady.data(), scratch.data() + size - hop, (size_t) hop * sizeof (double));
    }

    Fft64 fft;

    int size = 0, taps = 0, hop = 0;
    int pendingCount = 0, readyCount = 0;

    std::vector<double> block, scratch;
    std::vector<double> inputPending, outputReady;
    std::vector<std::complex<double>> bins, response;
};

} // namespace lime
