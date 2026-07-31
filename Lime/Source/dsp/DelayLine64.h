/*
  ==============================================================================

    DelayLine64 — plain integer-sample delay, double precision.

    Used only to align the three signal paths, which have different intrinsic
    latencies. No interpolation: every delay in this project is a whole number of
    samples by construction.

  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <vector>

namespace lime
{

class DelayLine64
{
public:
    void prepare (int maxDelaySamples)
    {
        capacity = std::max (1, maxDelaySamples + 1);
        buffer.assign ((size_t) capacity, 0.0);
        writeIndex = 0;
        delay = std::min (delay, capacity - 1);
    }

    /** Changing the delay keeps what the line is holding, but silences the span the read
        position steps back across when the delay grows.

        That span is where stale audio hides. At zero delay process() returns without
        writing anything, so the ring holds whatever was last put through it — and the
        engine runs at zero padding in loop mode. Switching from B to A took the delay
        from 0 to the full reported latency, and a line that simply moved its read
        position answered with a quarter of a second of audio from whenever it had last
        been in use. Zeroing exactly the samples the longer delay newly exposes keeps
        that from ever re-emerging — nothing deeper than the read position is read before
        being rewritten — without discarding the fresh audio the line was just given, so
        the dropout lasts only as long as the delay grew by, not the length of the line.
        Shortening the delay exposes nothing: the read position moves toward the write
        position and the samples in between are simply skipped. Callers set the delay
        every block with the same value, so none of this happens unless it genuinely
        moves. */
    void setDelay (int samples)
    {
        const int wanted = std::clamp (samples, 0, capacity - 1);

        if (wanted == delay)
            return;

        for (int depth = delay + 1; depth <= wanted; ++depth)
        {
            int index = writeIndex - depth;

            if (index < 0)
                index += capacity;

            buffer[(size_t) index] = 0.0;
        }

        delay = wanted;
    }

    int getDelay() const noexcept { return delay; }

    void reset()
    {
        std::fill (buffer.begin(), buffer.end(), 0.0);
        writeIndex = 0;
    }

    /** In-place delay of a block. Zero delay is a no-op, not a copy. */
    void process (double* samples, int numSamples)
    {
        if (delay == 0)
            return;

        for (int n = 0; n < numSamples; ++n)
        {
            buffer[(size_t) writeIndex] = samples[n];

            int readIndex = writeIndex - delay;

            if (readIndex < 0)
                readIndex += capacity;

            samples[n] = buffer[(size_t) readIndex];

            if (++writeIndex >= capacity)
                writeIndex = 0;
        }
    }

private:
    std::vector<double> buffer;
    int capacity = 1, writeIndex = 0, delay = 0;
};

} // namespace lime
