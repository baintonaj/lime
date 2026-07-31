/*
  ==============================================================================

    SpectrumProbe — wait-free snapshot of a magnitude spectrum for the editor.

    The audio thread writes whole frames; the message thread reads whole frames.
    Two buffers and an atomic index are not enough for that: the writer publishes
    every hop (5.3 ms at 48 kHz) while the reader walks a frame at display rate,
    so a two-slot scheme laps the reader and overwrites the frame mid-copy. The
    classic fix is a triple buffer: the writer owns one slot, the reader owns one
    slot, and the third sits in an atomic middle cell they trade through. Neither
    side ever touches a slot the other holds, so both sides are wait-free and a
    torn read cannot happen.

    Each slot carries its own geometry (bin count, FFT size, sample rate), so a
    frame claimed across a sample-rate change still describes itself correctly —
    the reader never sizes anything from a getter that another thread may be
    rewriting. All slot storage is allocated once, in the constructor, to the
    largest geometry the caller will ever ask for; prepare() never allocates.

    Nothing here is in the audio path when no editor is open — but only because it
    is switched off. push() is not cheap: it converts every bin's level and gain to
    decibels, two logarithms per bin per frame in both encode side chains. So the
    reader says when it is watching, and push() returns immediately when it is not.

  ==============================================================================
*/

#pragma once

#include <atomic>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace lime
{

class SpectrumProbe
{
public:
    /** One complete frame, as claimed by the reader. The vectors are sized by
        read() to exactly the frame's own bin count. */
    struct Snapshot
    {
        std::vector<double> levelDb, gainDb;
        int bins = 0;
        int fftSize = 0;
        double sampleRate = 0.0;
    };

    /** @param maxBins  the largest bin count any later prepare() may use. All
                        storage is allocated here, once. */
    explicit SpectrumProbe (int maxBins)
        : capacityBins (std::max (1, maxBins))
    {
        for (auto& slot : slots)
            slot.data.assign ((size_t) (2 * capacityBins), -160.0);   // interleaved level, gain
    }

    /** Writer side. Safe against a concurrent reader — it only touches atomics —
        but must not run concurrently with push(), which hosts already guarantee:
        prepareToPlay and processBlock never overlap. */
    void prepare (int numBins, double sampleRate, int fftSize)
    {
        assert (numBins <= capacityBins);

        binsHint.store (std::min (numBins, capacityBins), std::memory_order_relaxed);
        rateHint.store (sampleRate, std::memory_order_relaxed);
        fftSizeHint.store (fftSize, std::memory_order_relaxed);

        hasData.store (false, std::memory_order_release);
    }

    /** Hints for display logic only. A frame's authoritative geometry travels
        inside the frame; never size a buffer from these. */
    int getNumBins() const noexcept { return binsHint.load (std::memory_order_relaxed); }
    double getSampleRate() const noexcept { return rateHint.load (std::memory_order_relaxed); }
    int getFftSize() const noexcept { return fftSizeHint.load (std::memory_order_relaxed); }

    bool isReady() const noexcept { return hasData.load (std::memory_order_acquire); }

    /** Message thread. Set while an editor exists; until then push() does nothing. */
    void setActive (bool shouldBeActive) noexcept
    {
        active.store (shouldBeActive, std::memory_order_release);
    }

    /** Audio thread. `levels` and `gains` are numBins entries; gains are linear
        transfer factors and are stored as dB. */
    void push (const double* levels, const double* gains, int numBins,
               double sampleRate, int fftSize) noexcept
    {
        if (! active.load (std::memory_order_acquire))
            return;

        auto& slot = slots[(size_t) writeSlot];
        const int count = std::min (numBins, capacityBins);

        for (int k = 0; k < count; ++k)
        {
            const double level = levels != nullptr ? levels[k] : 0.0;
            const double gain = gains != nullptr ? gains[k] : 1.0;

            slot.data[(size_t) (2 * k)] = level > 0.0 ? 20.0 * std::log10 (level) : -160.0;
            slot.data[(size_t) (2 * k + 1)] = 20.0 * std::log10 (std::max (1.0e-8, gain));
        }

        slot.bins = count;
        slot.fftSize = fftSize;
        slot.sampleRate = sampleRate;

        // Trade the finished slot for whatever is in the middle. acq_rel: the release
        // half publishes the slot's contents with its index; the acquire half takes
        // ownership of the returned slot after any reads the reader made of it.
        writeSlot = middle.exchange (writeSlot | freshBit, std::memory_order_acq_rel) & indexMask;

        hasData.store (true, std::memory_order_release);
    }

    /** Message thread, single reader. Claims the freshest published frame — or
        re-reads the frame it already holds when nothing new has been published —
        and copies it into `out`, sizing the vectors to the frame's own bin count.
        Returns false until the writer has published anything. */
    bool read (Snapshot& out)
    {
        if (! hasData.load (std::memory_order_acquire))
            return false;

        // Only trade when the middle really holds a fresh frame. The flag is set
        // solely by the writer's exchange and cleared solely by ours, so a true
        // check cannot be invalidated in between — at worst the writer swaps in
        // something even fresher first.
        if ((middle.load (std::memory_order_relaxed) & freshBit) != 0)
            readSlot = middle.exchange (readSlot, std::memory_order_acq_rel) & indexMask;

        const auto& slot = slots[(size_t) readSlot];

        if (slot.bins <= 0)
            return false;

        out.bins = slot.bins;
        out.fftSize = slot.fftSize;
        out.sampleRate = slot.sampleRate;

        out.levelDb.resize ((size_t) slot.bins);
        out.gainDb.resize ((size_t) slot.bins);

        for (int k = 0; k < slot.bins; ++k)
        {
            out.levelDb[(size_t) k] = slot.data[(size_t) (2 * k)];
            out.gainDb[(size_t) k] = slot.data[(size_t) (2 * k + 1)];
        }

        return true;
    }

private:
    struct Slot
    {
        std::vector<double> data;   // interleaved level, gain
        int bins = 0;
        int fftSize = 0;
        double sampleRate = 0.0;
    };

    static constexpr int indexMask = 0x3;
    static constexpr int freshBit = 0x4;

    Slot slots[3];
    const int capacityBins;

    int writeSlot = 0;               // audio thread only
    int readSlot = 1;                // message thread only
    std::atomic<int> middle { 2 };   // index in bits 0-1, freshness in bit 2

    std::atomic<bool> hasData { false };
    std::atomic<bool> active { false };

    std::atomic<int> binsHint { 0 };
    std::atomic<int> fftSizeHint { 0 };
    std::atomic<double> rateHint { 48000.0 };
};

} // namespace lime
