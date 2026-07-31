/*
  ==============================================================================

    Wola64 — windowed overlap-add STFT engine, double precision.

    Used for the dynamic side chains, where the per-bin gain changes every frame
    and windowed overlap-add is the correct tool. Fixed filters belong in
    OverlapSave64 instead, which is exact.

    Latency is exactly `size` samples. A frame covers the most recent `size`
    inputs and the oldest hop of its overlap-add contribution is complete as soon
    as that frame has been added, which alone would cost size - hop; the engine
    consumes a hop before producing one, which adds the remaining hop. Trading
    that hop back would require priming logic that breaks on partial blocks, and
    the delay lines that align the three signal paths absorb it for free.
    Callers should read getLatencySamples() rather than assuming.

  ==============================================================================
*/

#pragma once

#include "Fft64.h"
#include "Window.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstring>
#include <vector>

namespace lime
{

class Wola64
{
public:
    Wola64 (int fftOrder, int hopIn, FftBackend backend = FftBackend::automatic)
        : fft (fftOrder, backend),
          window (1 << fftOrder, hopIn),
          size (1 << fftOrder),
          hop (hopIn)
    {
        // The contract is isValid(): a hop outside [1, size] would slide the
        // history by a negative count and size the staging buffers below from a
        // negative number, and a hop that does not divide the frame cannot
        // overlap-add flat. Assert it, then clamp so a release build handed bad
        // geometry degrades to legal buffer arithmetic instead of corrupting
        // memory; valid geometry passes through untouched.
        assert (isValid());
        hop = std::min (std::max (hop, 1), size);

        history.assign ((size_t) size, 0.0);
        overlap.assign ((size_t) size, 0.0);
        frame.assign ((size_t) size, 0.0);
        bins.assign ((size_t) fft.getNumBins(), { 0.0, 0.0 });

        inputPending.assign ((size_t) hop, 0.0);
        outputReady.assign ((size_t) hop, 0.0);
        pendingCount = 0;
        readyCount = hop;   // start with a hop of silence to establish the delay
    }

    void reset()
    {
        std::fill (history.begin(), history.end(), 0.0);
        std::fill (overlap.begin(), overlap.end(), 0.0);
        std::fill (inputPending.begin(), inputPending.end(), 0.0);
        std::fill (outputReady.begin(), outputReady.end(), 0.0);
        pendingCount = 0;
        readyCount = hop;
    }

    int getSize()   const noexcept { return size; }
    int getHop()    const noexcept { return hop; }
    int getNumBins() const noexcept { return fft.getNumBins(); }

    int getLatencySamples() const noexcept { return size; }

    /** The geometry this engine requires: the buffer bookkeeping needs a hop in
        [1, size], and exact reconstruction further needs the hop to divide the
        frame — the window pair's overlap-add sum cannot be phase-invariant
        otherwise (colaDeviation() reports how close it comes). */
    bool isValid() const noexcept { return hop >= 1 && hop <= size && size % hop == 0; }

    double colaDeviation() const noexcept { return window.colaDeviation(); }

    /** Divide a bin magnitude by this to get amplitude on the same scale as the time
        domain: a full-scale sine at a bin centre then reads 1.0. */
    double magnitudeToAmplitude() const noexcept { return 2.0 / window.analysisSum(); }
    FftBackend getBackend() const noexcept { return fft.getBackend(); }

    /** Streams a block of any length through the STFT.

        `processFrame` is called once per complete frame with the spectrum
        (getNumBins() complex values) to be modified in place. Input and output
        may alias.
    */
    template <typename FrameFn>
    void process (const double* input, double* output, int numSamples, FrameFn&& processFrame)
    {
        int done = 0;

        while (done < numSamples)
        {
            const int room = hop - pendingCount;
            const int avail = std::min (room, numSamples - done);

            std::memcpy (inputPending.data() + pendingCount, input + done, (size_t) avail * sizeof (double));

            // Emit from the hop of output produced by the previous frame.
            const int emitFrom = hop - readyCount;
            std::memcpy (output + done, outputReady.data() + emitFrom, (size_t) avail * sizeof (double));

            pendingCount += avail;
            readyCount -= avail;
            done += avail;

            if (pendingCount == hop)
            {
                runFrame (std::forward<FrameFn> (processFrame));
                pendingCount = 0;
                readyCount = hop;
            }
        }
    }

private:
    template <typename FrameFn>
    void runFrame (FrameFn&& processFrame)
    {
        // Slide the analysis history along by one hop and append the new samples.
        std::memmove (history.data(), history.data() + hop, (size_t) (size - hop) * sizeof (double));
        std::memcpy (history.data() + size - hop, inputPending.data(), (size_t) hop * sizeof (double));

        const auto* wa = window.analysis();
        const auto* ws = window.synthesis();

        for (int n = 0; n < size; ++n)
            frame[(size_t) n] = history[(size_t) n] * wa[n];

        fft.forward (frame.data(), bins.data());

        processFrame (bins.data(), fft.getNumBins());

        fft.inverse (bins.data(), frame.data());

        for (int n = 0; n < size; ++n)
            overlap[(size_t) n] += frame[(size_t) n] * ws[n];

        // The oldest hop of the overlap buffer is now complete.
        std::memcpy (outputReady.data(), overlap.data(), (size_t) hop * sizeof (double));

        std::memmove (overlap.data(), overlap.data() + hop, (size_t) (size - hop) * sizeof (double));
        std::fill (overlap.end() - hop, overlap.end(), 0.0);
    }

    Fft64 fft;
    WolaWindow window;

    int size = 0, hop = 0;
    int pendingCount = 0, readyCount = 0;

    std::vector<double> history, overlap, frame;
    std::vector<double> inputPending, outputReady;
    std::vector<std::complex<double>> bins;
};

} // namespace lime
