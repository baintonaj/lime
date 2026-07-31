/*
  ==============================================================================

    SrSideChain — staggered stage cascade for one frequency half.

  ==============================================================================
*/

#include "SrSideChain.h"

#include "DspMath.h"

#include <algorithm>
#include <cmath>

namespace lime
{

namespace
{
    /** Section 2.5: thresholds of approximately -30, -48 and -62 dB below reference
        level. High level first, then mid, then low, following the order of the
        first- and second-stage adders in section 3.1. */
    constexpr double sideChainThresholdsDb[3] = { -30.0, -48.0, -62.0 };

    /** Section 1: "somewhat over 8 dB" per stage, giving 24 dB at high frequencies
        and 16 dB at low. */
    constexpr double sideChainStageGainDb = 8.0;
}

//==============================================================================
void SrSideChain::prepare (bool highFrequency, int fftOrder, int hop, double sampleRate,
                           SlidingComposition composition, double slewFrameScale,
                           int maxBlockSamples)
{
    isHighFrequency = highFrequency;
    currentSampleRate = sampleRate;

    wola = std::make_unique<Wola64> (fftOrder, hop);

    const int numBins = wola->getNumBins();
    const int fftSize = wola->getSize();
    const double frameSeconds = (double) hop / sampleRate;

    // Power-complementary single-pole weights at the crossover, sampled onto the section
    // grid. |HP|^2 + |LP|^2 = 1 exactly, so the two halves' dB contributions add without a
    // bump. Built by the same routine that rebuilds them when the count moves, so there is
    // one formula rather than two that have to be kept in step.
    bandWeight.assign ((size_t) numBins, 0.0);
    binHz.assign ((size_t) numBins, 0.0);
    binLogPosition.assign ((size_t) numBins, 0.0);
    detectorWeight.assign ((size_t) numBins, 1.0);

    // Sized to the maximum rather than to the current count, so that changing the count
    // later never allocates and the setter stays safe to call from the audio thread.
    sectionWeight.assign ((size_t) SrSideChain::maxSections, 0.0);

    const double sectionSpan = std::log (sectionHiHz / sectionLoHz);

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * sampleRate / (double) fftSize;
        binHz[(size_t) k] = hz;

        // Bin 0 is DC, so the logarithm needs a floor; clamping to the grid's own ends is
        // what puts every out-of-band bin into the outermost section.
        const double clamped = std::clamp (hz, sectionLoHz, sectionHiHz);
        binLogPosition[(size_t) k] = std::log (clamped / sectionLoHz) / sectionSpan;
    }

    const int numStages = highFrequency ? 3 : 2;
    stages.resize ((size_t) numStages);

    for (int i = 0; i < numStages; ++i)
    {
        SrStageParams sp;
        sp.thresholdDb = sideChainThresholdsDb[i];
        sp.lowLevelGainDb = sideChainStageGainDb;

        // Each stage's threshold is quoted relative to the system input, but it sees
        // a signal already boosted by the stages ahead of it. The analog absorbs this
        // by raising the mid- and low-level stages' circuit gains; here it is
        // explicit, so a subthreshold signal gets full boost from every stage.
        sp.precedingGainDb = sideChainStageGainDb * (double) i;
        sp.highFrequency = highFrequency;
        sp.composition = composition;

        stages[(size_t) i].prepare (sp, numBins, sampleRate, fftSize, frameSeconds,
                                    bandWeight.data());
    }

    modulation.prepare (numBins, sampleRate, fftSize, frameSeconds);
    // Per-frame limits scaled to this half's frame, so both halves permit the same rate
    // of change in dB per second whatever hop they run at.
    GainSlewParams slewParams;
    slewParams.engagingDbPerFrame *= slewFrameScale;
    slewParams.releasingDbPerFrame *= slewFrameScale;

    slew.prepare (slewParams, numBins);

    lowLevelBoostDb.assign ((size_t) numBins, 0.0);

    // Fills bandWeight and lowLevelBoostDb, and pushes the weighting into the stages, which
    // were prepared above with a pointer to bandWeight. It has to come after both vectors
    // are sized: called earlier it wrote past the end of lowLevelBoostDb, which crashed on
    // the first prepare.
    rebuildBandWeights();

    // A detector corner set before prepare — or surviving a sample-rate change —
    // is honoured rather than lost, the same contract setSections keeps.
    rebuildDetectorWeights();

    // Per-bin coefficient for a smoother whose width is constant in Bark rather than
    // in bins, so the envelope is smoothed by a critical band everywhere.
    forwardCoeff.assign ((size_t) numBins, 0.0);

    for (int k = 0; k < numBins; ++k)
    {
        const double hz = (double) k * sampleRate / (double) fftSize;
        const double nextHz = (double) (k + 1) * sampleRate / (double) fftSize;
        const double barkStep = std::max (1.0e-9, barkFromHz (nextHz) - barkFromHz (hz));

        forwardCoeff[(size_t) k] = barkStep;
    }

    rebuildEnvelopeCoefficients();

    // Bin magnitudes are not on the time-domain amplitude scale: a tone of amplitude
    // A reads A * sum(window) / 2. The stage thresholds are quoted in dB below a
    // full-scale reference, so that factor has to come out or every threshold sits
    // some 52 dB too low and the compressors never release.
    magnitudeScale = wola->magnitudeToAmplitude();

    envelope.assign ((size_t) numBins, 0.0);
    // Filter-bank realisation: a Bark-spaced bank driven from the transform's
    // analysis. Deliberately undelayed — see the note in process().
    bank.prepare (sampleRate);

    bandTargetDb.assign ((size_t) bank.getNumBands(), 0.0);
    bandBin.assign ((size_t) bank.getNumBands(), 0);

    for (int b = 0; b < bank.getNumBands(); ++b)
    {
        const int k = (int) std::lround (bank.getCentreHz (b) * (double) fftSize / sampleRate);
        bandBin[(size_t) b] = std::clamp (k, 0, numBins - 1);
    }

    // Sized to the largest block the caller promised — no larger: eight of these
    // exist per stereo instance, and rounding every one up to a worst case nobody
    // asked for costs a megabyte of scratch at a typical block size. A block that
    // exceeds the promise anyway is taken in scratch-sized pieces by process().
    const auto scratchLen = (size_t) std::max (1, maxBlockSamples);

    analyserScratch.assign (scratchLen, 0.0);
    discardScratch.assign (scratchLen, 0.0);

    magnitudes.assign ((size_t) numBins, 0.0);
    running.assign ((size_t) numBins, 0.0);
    transfer.assign ((size_t) numBins, 0.0);
    total.assign ((size_t) numBins, 1.0);
    mcFeed.assign ((size_t) numBins, 0.0);
    pending.assign ((size_t) numBins, 0.0);

    reset();
}

void SrSideChain::setSections (int count)
{
    // Rounded down to a power of two rather than to the nearest, so that a caller handing
    // over a bin count gets a grid that fits inside it rather than one that overshoots.
    int rounded = minSections;

    while (rounded * 2 <= std::clamp (count, minSections, maxSections))
        rounded *= 2;

    // A rebuild is O(bins) across a handful of tables and only happens on an actual change,
    // which is nothing beside the transform it sits in front of — but a static control
    // should still cost nothing at all.
    if (rounded == sections)
        return;

    sections = rounded;
    rebuildBandWeights();
}

void SrSideChain::rebuildEnvelopeCoefficients()
{
    const auto numBins = forwardCoeff.size();

    envelopeCoeff.assign (numBins, 0.0);
    envelopeOneMinusCoeff.assign (numBins, 1.0);
    envelopePeakNormalisation.assign (numBins, 1.0);

    for (size_t k = 0; k < numBins; ++k)
    {
        const double barkWidth = envelopeBark > 0.0 ? envelopeBark / forwardCoeff[k] : 0.0;

        // Clamped to the spectrum's own width: a kernel wider than the whole bin
        // range smooths no differently from one exactly that wide, and the clamp is
        // what keeps the peak normalisation below finite in the next line — at the
        // widest permitted kernel 1 - a is about 1e-4, nowhere near underflow.
        const double width = std::min ((double) numBins,
                                       std::max (barkWidth, minimumSmoothingBins));
        const double a = width > 0.0 ? std::exp (-1.0 / width) : 0.0;

        envelopeCoeff[k] = a;
        envelopeOneMinusCoeff[k] = 1.0 - a;

        // For a two-pass one-pole with coefficient a the impulse response peaks at
        // (1 - a)/(1 + a); the reciprocal makes the kernel unity at its centre.
        envelopePeakNormalisation[k] = (1.0 + a) / (1.0 - a);
    }
}

void SrSideChain::rebuildBandWeights()
{
    const auto numBins = (int) bandWeight.size();

    // The profile evaluated once per section, at the section's geometric centre. Both
    // halves run this identical function of frequency, and the low half takes 1 - hp from
    // the same double, so |HP|^2 + |LP|^2 is 1 to the last bit at every frequency however
    // coarse the grid is. That is the property the whole cascade rests on: the two halves
    // multiply, so their decibels add, and anything less than exact here shows up as a
    // bump in the encode characteristic.
    for (int s = 0; s < sections; ++s)
    {
        const double t = ((double) s + 0.5) / (double) sections;
        const double hz = sectionLoHz * std::pow (sectionHiHz / sectionLoHz, t);

        const double r = hz / referenceCrossoverHz;
        const double highPassPower = (r * r) / (1.0 + r * r);

        sectionWeight[(size_t) s] = isHighFrequency ? highPassPower : 1.0 - highPassPower;
    }

    // No pow in the bin loop: the positions were built in prepare, so this is a multiply,
    // a truncation and a lookup per bin however many sections there are.
    for (int k = 0; k < numBins; ++k)
    {
        const int s = std::clamp ((int) (binLogPosition[(size_t) k] * (double) sections),
                                  0, sections - 1);

        bandWeight[(size_t) k] = sectionWeight[(size_t) s];
    }

    const auto numStages = (double) stages.size();

    for (int k = 0; k < numBins; ++k)
        lowLevelBoostDb[(size_t) k] = sideChainStageGainDb * bandWeight[(size_t) k] * numStages;

    for (auto& stage : stages)
        stage.setBandWeight (bandWeight.data());
}

void SrSideChain::setDetectorHighPass (double cornerHz)
{
    // A static control costs nothing: the weights rebuild only when a corner moves.
    if (cornerHz == detectorHighPassHz)
        return;

    detectorHighPassHz = cornerHz;
    rebuildDetectorWeights();
}

void SrSideChain::setDetectorLowPass (double cornerHz)
{
    if (cornerHz == detectorLowPassHz)
        return;

    detectorLowPassHz = cornerHz;
    rebuildDetectorWeights();
}

void SrSideChain::rebuildDetectorWeights()
{
    const auto numBins = (int) detectorWeight.size();

    // The two magnitudes multiply, so their decibels add: together the corners
    // band-limit the detection, and either alone is what it says it is.
    for (int k = 0; k < numBins; ++k)
    {
        const double hz = binHz[(size_t) k];
        double weight = 1.0;

        if (detectorHighPassHz > 0.0)
            weight *= butterworth3Magnitude (hz, detectorHighPassHz, true);

        if (detectorLowPassHz > 0.0)
            weight *= butterworth3Magnitude (hz, detectorLowPassHz, false);

        detectorWeight[(size_t) k] = weight;
    }
}

void SrSideChain::reset()
{
    if (wola != nullptr)
        wola->reset();

    for (auto& s : stages)
        s.reset();

    modulation.reset();
    slew.reset();
    bank.reset();

    std::fill (magnitudes.begin(), magnitudes.end(), 0.0);
    std::fill (mcFeed.begin(), mcFeed.end(), 0.0);
    std::fill (pending.begin(), pending.end(), 0.0);
    std::fill (total.begin(), total.end(), 1.0);

    slidingMc = 0.0;
    lastDominants = 0;
    haveState = false;
}

//==============================================================================
void SrSideChain::process (const double* input, double* output, int numSamples)
{
    if (realisation == GainRealisation::stft)
    {
        wola->process (input, output, numSamples,
                       [this] (std::complex<double>* bins, int numBins)
                       {
                           processFrame (bins, numBins);
                       });

        return;
    }

    // Filter-bank realisation. The scratch was sized once in prepare(), so a block
    // longer than the caller promised is taken in scratch-sized pieces rather than
    // grown for — nothing here may allocate on the audio thread.
    const auto capacity = (int) analyserScratch.size();

    for (int offset = 0; offset < numSamples; offset += capacity)
    {
        const int len = std::min (capacity, numSamples - offset);
        processFilterBankChunk (input + offset, output + offset, len);
    }
}

void SrSideChain::processFilterBankChunk (const double* input, double* output, int numSamples)
{
    // Order matters, and getting it wrong is what made the first attempt worse than the
    // STFT it replaced. Both directions must build their control state from the
    // *encoded* signal, which is the only one they hold identically — the encoder's
    // input is not available to the decoder, and the decoder's output is not available
    // to the encoder. So the gain is applied first, from state established by earlier
    // blocks, and the encoded signal is analysed afterwards. For the encoder the encoded
    // signal is its output; for the decoder it is its input.

    // Keep the encoded signal before the decoder overwrites it in place.
    if (direction == SrDirection::decode)
        std::copy (input, input + numSamples, analyserScratch.begin());

    // Solve the bank for the profile the previous blocks' analysis produced.
    for (int b = 0; b < bank.getNumBands(); ++b)
    {
        // Floored at -160 dB: a transfer this small only arises from degenerate
        // input, and the log below must not be handed zero.
        const double factor = std::max (1.0e-8, total[(size_t) bandBin[(size_t) b]]);
        bandTargetDb[(size_t) b] = 20.0 * std::log10 (factor);
    }

    bank.setTargetDb (bandTargetDb.data());

    if (output != input)
        std::copy (input, input + numSamples, output);

    // Deliberately *not* delayed to line up with the analysis. Delaying the signal makes
    // the decoder undo a block that was encoded with coefficients from one window earlier
    // while using its own current ones, and the pair stops cancelling — measured, that
    // cost 30 dB. Leaving both undelayed means both use coefficients that are equally
    // stale, which is what makes the inversion exact. The price is that the gain responds
    // one analysis window late; see the note on control lag in DESIGN.md.
    if (direction == SrDirection::encode)
    {
        bank.process (output, numSamples);
        std::copy (output, output + numSamples, analyserScratch.begin());
    }
    else
    {
        bank.processInverse (output, numSamples);
    }

    // Now analyse the encoded signal. Its synthesised output is discarded; only the
    // profile the frame function computes is wanted.
    wola->process (analyserScratch.data(), discardScratch.data(), numSamples,
                   [this] (std::complex<double>* bins, int numBins)
                   {
                       processFrame (bins, numBins);
                   });
}

void SrSideChain::processFrame (std::complex<double>* bins, int numBins)
{
    // Apply the transfer established by previous frames. Deriving it from `pending`
    // rather than from this frame's spectrum is what makes encode and decode exact
    // inverses: both build the same state from the same |u| history.
    if (haveState)
    {
        // Smooth the level envelope across critical bands before the controls see it.
        // Without this the gain profile notches individual bins and the encode cannot
        // be inverted at all; see setEnvelopeSmoothingBark.
        if (envelopeBark > 0.0 || minimumSmoothingBins > 0.0)
        {
            // Critical-band smoothing of the level envelope. Two properties matter and
            // the first attempt got the second wrong:
            //
            //   - It must be smooth, or the gain profile notches single bins and the
            //     encode cannot be inverted. See setEnvelopeSmoothingBark.
            //   - It must preserve a tone's level. An *averaging* smoother dilutes a
            //     tone's energy across the band, so the compressor under-reads it and
            //     leaves 5 dB of boost at reference level where the paper allows about
            //     one. Smoothing power with a peak-normalised kernel keeps the tone's
            //     own level intact while still spreading a skirt around it, which is
            //     also what a critical-band analysis does.
            //
            // For a two-pass one-pole with coefficient a the impulse response peaks at
            // (1 - a)/(1 + a), so scaling by its reciprocal makes the kernel unity at
            // its centre rather than unity in area.
            // The coefficients are a table, not a computation. They depend only on the
            // bin, the Bark width and the floor, none of which move between parameter
            // changes — but they used to be evaluated inside both passes, which put two
            // std::exp and a division on every bin of every frame of all four side
            // chains for a value that had not changed. See rebuildEnvelopeCoefficients.
            const double* const a = envelopeCoeff.data();
            const double* const oneMinusA = envelopeOneMinusCoeff.data();
            const double* const peakNorm = envelopePeakNormalisation.data();

            double state = pending[0] * pending[0];

            for (int k = 0; k < numBins; ++k)
            {
                const size_t i = (size_t) k;
                const double power = pending[i] * pending[i];

                state = a[i] * state + oneMinusA[i] * power;
                envelope[i] = state;
            }

            state = envelope[(size_t) (numBins - 1)];

            for (int k = numBins - 1; k >= 0; --k)
            {
                const size_t i = (size_t) k;

                state = a[i] * state + oneMinusA[i] * envelope[i];
                envelope[i] = std::sqrt (std::max (0.0, state * peakNorm[i]));
            }

            std::copy (envelope.begin(), envelope.begin() + numBins, running.begin());
        }
        else
        {
            std::copy (pending.begin(), pending.begin() + numBins, running.begin());
        }

        // The user's side-chain filters. The detection reads through them — outside
        // the corners a bin looks quieter than it is, so the compressors treat it as
        // low-level programme — while `pending` stays the raw encoded magnitude, so
        // the shared encode/decode state is untouched and a corner change never
        // desynchronises the two directions. Weighting `running` here also reaches
        // the MC feed, which in the analog is itself taken from frequency-weighted
        // main-path signals.
        if (detectorHighPassHz > 0.0 || detectorLowPassHz > 0.0)
            for (int k = 0; k < numBins; ++k)
                running[(size_t) k] *= detectorWeight[(size_t) k];

        std::fill (total.begin(), total.begin() + numBins, 1.0);

        const int numStages = (int) stages.size();
        lastDominants = 0;

        for (int s = 0; s < numStages; ++s)
        {
            auto& stage = stages[(size_t) s];

            stage.setModulationControl (slidingMc);
            stage.processFrame (running.data(), transfer.data());

            lastDominants += stage.getNumDominants();

            for (int k = 0; k < numBins; ++k)
            {
                const double factor = 1.0 + transfer[(size_t) k];

                total[(size_t) k] *= factor;
                running[(size_t) k] *= factor;
            }

            // Section 3.2 takes MC1 to MC7 from the second-stage adder output.
            if (s == std::min (1, numStages - 1))
                std::copy (running.begin(), running.begin() + numBins, mcFeed.begin());
        }

        // Only the sliding opposition signal this half consumes; the other seven
        // MC signals are a design completeness the signal path never reads. See
        // ModulationControl::processFrameEssential.
        slidingMc = modulation.processFrameEssential (mcFeed.data(), isHighFrequency);

        // Limit how fast the composite transfer may move, in place of the analog
        // overshoot suppressors. The threshold tracks the modulation control signal,
        // as section 2.4 specifies, and both directions run it on identical input so
        // complementarity is untouched.
        if (slewEnabled)
        {
            for (int k = 0; k < numBins; ++k)
                transfer[(size_t) k] = total[(size_t) k] - 1.0;

            slew.setModulationControl (slidingMc);
            slew.process (transfer.data());

            for (int k = 0; k < numBins; ++k)
                total[(size_t) k] = 1.0 + transfer[(size_t) k];
        }

        // Smooth the composite transfer across frequency. Done in dB so the shape is
        // smoothed rather than its linear amplitude, and run forward then backward so
        // the kernel is symmetric and introduces no spectral tilt.
        if (smoothingBins > 0.0)
        {
            const double a = std::exp (-1.0 / smoothingBins);

            double state = 20.0 * std::log10 (total[0]);

            for (int k = 0; k < numBins; ++k)
            {
                const double db = 20.0 * std::log10 (total[(size_t) k]);
                state = a * state + (1.0 - a) * db;
                transfer[(size_t) k] = state;
            }

            state = transfer[(size_t) (numBins - 1)];

            for (int k = numBins - 1; k >= 0; --k)
            {
                state = a * state + (1.0 - a) * transfer[(size_t) k];
                total[(size_t) k] = std::pow (10.0, state / 20.0);
            }
        }
    }

    // Both directions take their next-frame state from the *encoded* spectrum, which
    // is the one signal both sides hold bit-identically. Taking it from the source
    // spectrum instead is what the first attempt did, and it diverges: the decoder
    // would have to reconstruct the source by division, and an STFT per-bin
    // modification is only approximately invertible, so that small error feeds the
    // control state and compounds frame after frame. Measured, it dragged the null
    // from the -126 dB the transform itself allows up to -22 dB.
    // std::abs on a complex dispatches to hypot, which spends its time on the overflow
    // and underflow guards a general-purpose routine needs. A bin magnitude here is an
    // audio spectrum, tens of orders of magnitude away from either limit, so the plain
    // square root is both correct and around an order of magnitude cheaper — and unlike
    // a library call it leaves the loop open to the vectoriser.
    const auto magnitudeOf = [] (const std::complex<double>& z) noexcept
    {
        const double re = z.real(), im = z.imag();
        return std::sqrt (re * re + im * im);
    };

    if (realisation == GainRealisation::filterBank)
    {
        // The bins are left alone: the profile is applied in the time domain. State is
        // still taken from the encoded spectrum, which in this mode means the analysed
        // signal already carries whatever the bank applied on previous blocks.
        for (int k = 0; k < numBins; ++k)
            magnitudes[(size_t) k] = magnitudeOf (bins[k]) * magnitudeScale;
    }
    else if (direction == SrDirection::encode)
    {
        for (int k = 0; k < numBins; ++k)
        {
            bins[k] *= total[(size_t) k];
            magnitudes[(size_t) k] = magnitudeOf (bins[k]) * magnitudeScale;
        }
    }
    else
    {
        for (int k = 0; k < numBins; ++k)
        {
            magnitudes[(size_t) k] = magnitudeOf (bins[k]) * magnitudeScale;

            // Left as a division deliberately. A reciprocal and two multiplies would be
            // faster, but the encode multiplies by exactly this number, and x * t * (1/t)
            // is not x. The complementarity is the point of the whole engine; two divides
            // a bin is not where the time goes.
            const double t = total[(size_t) k];
            bins[k] /= (t > 0.0 ? t : 1.0);
        }
    }

    std::copy (magnitudes.begin(), magnitudes.begin() + numBins, pending.begin());

    // Publish what the display needs: the level this half saw and the transfer it
    // actually applied — the boost on encode, the reciprocal cut on decode, so a
    // decode-only instance shows its expansion below the line instead of freezing.
    // The engine wires the probe to whichever direction is running; see
    // SrEngine::process. The probe outlives every re-prepare; see setProbe.
    if (probe != nullptr)
    {
        const double* gains = total.data();

        if (direction == SrDirection::decode)
        {
            // `transfer` is scratch between frames; the applied factor is 1/total.
            for (int k = 0; k < numBins; ++k)
            {
                const double t = total[(size_t) k];
                transfer[(size_t) k] = t > 0.0 ? 1.0 / t : 1.0;
            }

            gains = transfer.data();
        }

        probe->push (magnitudes.data(), gains, numBins,
                     currentSampleRate, wola != nullptr ? wola->getSize() : 0);
    }

    haveState = true;
}

} // namespace lime
