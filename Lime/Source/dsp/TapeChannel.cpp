/*
  ==============================================================================

    TapeChannel — coefficient design, noise generation and the per-sample chain.

    Chain order follows the machine: record equalisation and the coating's
    gradual saturation, then the modulation noise the coating adds to what it has
    recorded, then the replay head's gap loss, then the replay hiss. Only the
    hiss sits after the loss, because hiss is generated at replay and is not
    subject to it.

    Every noise level is normalised through the measured power gain of its own
    shaping filter, so what is asked for in dB is what a test measures, whatever
    the shape or the sample rate.

  ==============================================================================
*/

#include "TapeChannel.h"

#include <algorithm>
#include <complex>

namespace lime
{

namespace
{

// Prefixed against the plain spellings other files in this directory use, so the
// sources stay safe to concatenate into one translation unit — some build styles
// unity-build a directory like this one, and identical file-local names are the
// first thing that breaks.
constexpr double tapePi = 3.14159265358979323846;
constexpr double tapeTwoPi = 6.283185307179586476925286766559;

//==============================================================================
/** xoshiro256** with a splitmix64 seeder, and Box-Muller on top.

    A local generator rather than <random> for two reasons: the sequence must be
    reproducible across standard library versions for the measured figures to
    mean anything, and nothing here may allocate or branch unpredictably.
*/
class NoiseGenerator
{
public:
    void seed (uint64_t s) noexcept
    {
        // splitmix64 expansion; any seed, including zero, gives a usable state.
        for (auto& word : state)
        {
            s += 0x9E3779B97F4A7C15ull;
            uint64_t z = s;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
            word = z ^ (z >> 31);
        }

        hasSpare = false;
        spare = 0.0;
    }

    /** Unit-variance Gaussian. Consumes two uniforms per two values, so the
        stream advances at a fixed rate. */
    double nextGaussian() noexcept
    {
        if (hasSpare)
        {
            hasSpare = false;
            return spare;
        }

        const double u1 = nextUniform();
        const double u2 = nextUniform();
        const double radius = std::sqrt (-2.0 * std::log (u1));
        const double angle = tapeTwoPi * u2;

        spare = radius * std::sin (angle);
        hasSpare = true;

        return radius * std::cos (angle);
    }

private:
    /** In (0, 1]: never zero, so the logarithm above is always finite. */
    double nextUniform() noexcept
    {
        return ((double) (nextBits() >> 11) + 0.5) * 0x1.0p-53;
    }

    uint64_t nextBits() noexcept
    {
        const uint64_t result = rotl (state[1] * 5ull, 7) * 9ull;
        const uint64_t t = state[1] << 17;

        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];
        state[2] ^= t;
        state[3] = rotl (state[3], 45);

        return result;
    }

    static uint64_t rotl (uint64_t x, int k) noexcept
    {
        return (x << k) | (x >> (64 - k));
    }

    uint64_t state[4] { 0, 0, 0, 0 };
    double spare = 0.0;
    bool hasSpare = false;
};

//==============================================================================
/** Bilinear prewarp tangent, with the critical frequency held inside the band so
    the tangent cannot blow up. */
double prewarpTangent (double hz, double sampleRate)
{
    const double limit = 0.49 * sampleRate;
    return std::tan (tapePi * std::clamp (hz, 1.0e-3, limit) / sampleRate);
}

/** One-pole low shelf: `gainDb` at DC, unity at Nyquist, -3 dB corner at
    `cornerHz`. 6 dB/octave, which is the slope the tape equalisations use. */
BiquadCoeffs onePoleLowShelf (double cornerHz, double gainDb, double sampleRate)
{
    const double a = std::pow (10.0, gainDb / 20.0);
    const double t = prewarpTangent (cornerHz, sampleRate);
    const double d = 1.0 + t;

    BiquadCoeffs c;
    c.b0 = (1.0 + a * t) / d;
    c.b1 = (a * t - 1.0) / d;
    c.a1 = (t - 1.0) / d;

    return c;
}

/** One-pole high shelf: unity at DC, `gainDb` at Nyquist. */
BiquadCoeffs onePoleHighShelf (double cornerHz, double gainDb, double sampleRate)
{
    const double a = std::pow (10.0, gainDb / 20.0);
    const double t = prewarpTangent (cornerHz, sampleRate);
    const double d = 1.0 + t;

    BiquadCoeffs c;
    c.b0 = (a + t) / d;
    c.b1 = (t - a) / d;
    c.a1 = (t - 1.0) / d;

    return c;
}

/** One-pole high pass. */
BiquadCoeffs onePoleHighPass (double cornerHz, double sampleRate)
{
    const double t = prewarpTangent (cornerHz, sampleRate);
    const double d = 1.0 + t;

    BiquadCoeffs c;
    c.b0 = 1.0 / d;
    c.b1 = -1.0 / d;
    c.a1 = (t - 1.0) / d;

    return c;
}

/** Two-pole Butterworth low pass, exactly -3.01 dB at `cornerHz`. */
BiquadCoeffs butterworthLowPass (double cornerHz, double sampleRate)
{
    const double limit = 0.49 * sampleRate;
    const double w0 = tapeTwoPi * std::clamp (cornerHz, 1.0e-3, limit) / sampleRate;
    const double cosW = std::cos (w0);
    const double alpha = std::sin (w0) * 0.70710678118654752440;   // Q = 1/sqrt(2)
    const double a0 = 1.0 + alpha;

    BiquadCoeffs c;
    c.b0 = 0.5 * (1.0 - cosW) / a0;
    c.b1 = (1.0 - cosW) / a0;
    c.b2 = c.b0;
    c.a1 = -2.0 * cosW / a0;
    c.a2 = (1.0 - alpha) / a0;

    return c;
}

/** Mean of |H|^2 over the whole band, which is the power gain a cascade applies
    to white noise. Midpoint rule; 2048 points puts the quadrature error far below
    the sampling error of any noise measurement made against it. */
double meanSquareResponse (const std::vector<BiquadCoeffs>& coeffs)
{
    constexpr int numPoints = 2048;
    double sum = 0.0;

    for (int m = 0; m < numPoints; ++m)
    {
        const double omega = tapePi * ((double) m + 0.5) / (double) numPoints;
        std::complex<double> h { 1.0, 0.0 };

        for (const auto& c : coeffs)
            h *= c.responseAt (omega);

        sum += std::norm (h);
    }

    return sum / (double) numPoints;
}

} // namespace

//==============================================================================
struct TapeChannel::Channel
{
    NoiseGenerator hissNoise, modulationNoise;

    BiquadCascade64 hissShape;         // replay-shaped hiss spectrum
    BiquadCascade64 preEmphasis;       // record equalisation, into the saturator
    BiquadCascade64 deEmphasis;        // its exact inverse, out of the saturator
    BiquadCascade64 headLoss;          // replay gap loss
    BiquadCascade64 modulationShape;   // band limit of the modulating process

    void reset()
    {
        hissShape.reset();
        preEmphasis.reset();
        deEmphasis.reset();
        headLoss.reset();
        modulationShape.reset();
    }
};

//==============================================================================
TapeChannel::TapeChannel() = default;
TapeChannel::~TapeChannel() = default;

void TapeChannel::prepare (double sampleRate, int /*maxBlockSize*/, int numChannels)
{
    currentSampleRate = (sampleRate > 0.0) ? sampleRate : 48000.0;

    channels.clear();
    channels.reserve ((size_t) std::max (0, numChannels));

    for (int ch = 0; ch < numChannels; ++ch)
        channels.push_back (std::make_unique<Channel>());

    // Size the coefficient lists once, so later redesigns only overwrite them.
    hissShapeCoeffs.assign (3, BiquadCoeffs {});
    emphasisCoeffs.assign (2, BiquadCoeffs {});
    headLossCoeffs.assign (1, BiquadCoeffs {});
    modulationShapeCoeffs.assign (2, BiquadCoeffs {});

    prepared = true;

    updateCoefficients();
    applyCoefficients();
    seedChannels();

    for (auto& c : channels)
        c->reset();
}

void TapeChannel::reset()
{
    for (auto& c : channels)
        c->reset();

    seedChannels();
}

void TapeChannel::setSettings (const TapeChannelSettings& newSettings)
{
    settings = newSettings;

    if (! prepared)
        return;

    updateCoefficients();
    applyCoefficients();
}

void TapeChannel::setSeed (uint64_t seed)
{
    baseSeed = seed;
    seedChannels();
}

//==============================================================================
double TapeChannel::effectiveHfLossHz() const noexcept
{
    const double speed = std::clamp (settings.tapeSpeedIps,
                                     tapeModel::minSpeedIps,
                                     tapeModel::maxSpeedIps);
    const double scaled = std::max (20.0, settings.hfLossHz)
                              * speed / tapeModel::referenceSpeedIps;

    if (! prepared)
        return scaled;

    return std::min (scaled, 0.49 * currentSampleRate);
}

//==============================================================================
void TapeChannel::updateCoefficients()
{
    const double fs = currentSampleRate;
    const double speedRatio = std::clamp (settings.tapeSpeedIps,
                                          tapeModel::minSpeedIps,
                                          tapeModel::maxSpeedIps)
                                  / tapeModel::referenceSpeedIps;

    // Record equalisation. The low corner is speed independent (3180 us at both
    // speeds); the high corner scales, which is why 30 ips has more high
    // frequency headroom than 15.
    const double recordHighHz = std::min (tapeModel::ccirHighCornerHz * speedRatio,
                                          0.45 * fs);

    emphasisCoeffs[0] = onePoleLowShelf (tapeModel::nabLowCornerHz,
                                         tapeModel::recordLowBoostDb, fs);
    emphasisCoeffs[1] = onePoleHighShelf (recordHighHz,
                                          tapeModel::recordHighBoostDb, fs);

    // Saturator scaling. The drive is divided back out afterwards, so extra drive
    // moves the onset of saturation down rather than moving the operating level
    // up: the loop stays gain-calibrated, which is what reference alignment assumes.
    const double knee = tapeModel::referenceAmplitude
                            * std::pow (10.0, tapeModel::midbandHeadroomDb / 20.0);
    const double drive = std::pow (10.0, settings.saturationDb / 20.0);

    saturationInScale = drive / knee;
    saturationOutScale = knee / drive;

    // Head and gap loss.
    headLossCoeffs[0] = butterworthLowPass (effectiveHfLossHz(), fs);

    // Hiss shape: replay lift at both ends, bandwidth limit at the top. Shape
    // moves with tape speed, total level does not, because level is what the
    // caller asked for.
    hissShapeCoeffs[0] = onePoleLowShelf (tapeModel::hissLowShelfHz,
                                          tapeModel::hissLowShelfDb, fs);
    hissShapeCoeffs[1] = onePoleHighShelf (std::min (tapeModel::ccirHighCornerHz * speedRatio,
                                                     0.45 * fs),
                                           tapeModel::hissHighShelfDb, fs);
    hissShapeCoeffs[2] = butterworthLowPass (std::min (tapeModel::hissBandwidthHz * speedRatio,
                                                       0.45 * fs), fs);

    const double hissRms = tapeModel::referenceAmplitude
                               * std::pow (10.0, settings.hissLevelDb / 20.0);
    const double hissPower = meanSquareResponse (hissShapeCoeffs);

    hissGain = (hissPower > 0.0) ? hissRms / std::sqrt (hissPower) : 0.0;

    // Modulating process: band limited, and high-passed so its effect is noise
    // around the signal rather than a slow drift of the signal's level.
    modulationShapeCoeffs[0] = onePoleHighPass (tapeModel::modulationHighPassHz, fs);
    modulationShapeCoeffs[1] = butterworthLowPass (std::min (tapeModel::modulationSpreadHz,
                                                             0.45 * fs), fs);

    const double modulationRms = std::pow (10.0, settings.modulationNoiseDb / 20.0);
    const double modulationPower = meanSquareResponse (modulationShapeCoeffs);

    modulationGain = (modulationPower > 0.0) ? modulationRms / std::sqrt (modulationPower)
                                             : 0.0;
}

void TapeChannel::applyCoefficients()
{
    for (auto& c : channels)
    {
        c->hissShape.setCoeffs (hissShapeCoeffs);
        c->preEmphasis.setCoeffs (emphasisCoeffs);
        c->deEmphasis.setInverseCoeffs (emphasisCoeffs);
        c->headLoss.setCoeffs (headLossCoeffs);
        c->modulationShape.setCoeffs (modulationShapeCoeffs);
    }
}

void TapeChannel::seedChannels()
{
    // Distinct, well-separated streams per channel and per defect, so that no two
    // are correlated and enabling one never disturbs another's sequence.
    for (size_t ch = 0; ch < channels.size(); ++ch)
    {
        const uint64_t channelSeed = baseSeed + 0x9E3779B97F4A7C15ull * (uint64_t) (ch + 1);

        channels[ch]->hissNoise.seed (channelSeed);
        channels[ch]->modulationNoise.seed (channelSeed ^ 0xD1B54A32D192ED03ull);
    }
}

//==============================================================================
void TapeChannel::process (double* const* channelData, int numChannels, int numSamples)
{
    if (! prepared || numSamples <= 0)
        return;

    const int usable = std::min (numChannels, (int) channels.size());

    // The enable flags cannot change inside a block — setSettings is the only
    // writer — so read them once rather than once per sample.
    const bool saturationOn = settings.saturationEnabled;
    const bool modulationOn = settings.modulationNoiseEnabled;
    const bool hfLossOn = settings.hfLossEnabled;
    const bool hissOn = settings.hissEnabled;

    for (int ch = 0; ch < usable; ++ch)
    {
        auto& c = *channels[(size_t) ch];
        double* io = channelData[ch];

        if (io == nullptr)
            continue;

        for (int n = 0; n < numSamples; ++n)
        {
            double x = io[n];

            // Gradual saturation, made frequency dependent by the equalisation
            // pair around it. Small signals emerge unchanged because the pair is
            // an exact inverse and the curve has unit slope at the origin.
            if (saturationOn)
            {
                const double driven = c.preEmphasis.processSample (x) * saturationInScale;
                x = c.deEmphasis.processSample (tapeSaturationCurve (driven) * saturationOutScale);
            }

            // Modulation noise: multiplicative, so it vanishes with the signal and
            // its sidebands sit around whatever the signal contains.
            if (modulationOn)
            {
                const double m = c.modulationShape.processSample (c.modulationNoise.nextGaussian());
                x += x * m * modulationGain;
            }

            if (hfLossOn)
                x = c.headLoss.processSample (x);

            if (hissOn)
                x += c.hissShape.processSample (c.hissNoise.nextGaussian()) * hissGain;

            io[n] = x;
        }
    }

    // Channels beyond those prepared are left untouched rather than zeroed.
}

} // namespace lime
