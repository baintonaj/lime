/*
  ==============================================================================

    Calibration — implementation.

  ==============================================================================
*/

#include "Calibration.h"

#include <algorithm>
#include <cmath>
#include <complex>

namespace lime
{

namespace
{

constexpr double pi = 3.141592653589793238462643383279;
constexpr double twoPi = 6.283185307179586476925286766559;

/** Rounds a duration in seconds to a whole number of samples, never zero. */
int samplesFor (double seconds, double sampleRate)
{
    const int n = (int) std::lround (seconds * sampleRate);
    return n > 0 ? n : 1;
}

/** One-pole smoothing coefficient for a time constant in seconds. */
double smoothingCoeff (double seconds, double sampleRate)
{
    if (seconds <= 0.0)
        return 1.0;

    return 1.0 - std::exp (-1.0 / (seconds * sampleRate));
}

/** First-order section (1 + s/wz) / (1 + s/wp), bilinear transformed with
    c = 2*fs, where wp and wz are already prewarped analog angular frequencies.

    Unity at DC and exactly wp/wz at Nyquist, both by construction, so a cascade
    of these has a flat passband no matter how many sections it has. Stable for
    any positive wp, since a1 = (wp - c)/(wp + c) has magnitude below one.
*/
BiquadCoeffs firstOrderShelf (double wp, double wz, double sampleRate)
{
    const double c = 2.0 * sampleRate;
    const double k = wp / wz;

    BiquadCoeffs out;
    out.b0 = k * (wz + c) / (wp + c);
    out.b1 = k * (wz - c) / (wp + c);
    out.b2 = 0.0;
    out.a1 = (wp - c) / (wp + c);
    out.a2 = 0.0;

    return out;
}

/** Prewarp: analog angular frequency that the plain bilinear with c = 2*fs maps
    onto the digital frequency `hz`. */
double prewarp (double hz, double sampleRate)
{
    return 2.0 * sampleRate * std::tan (pi * hz / sampleRate);
}

double cascadeDb (const std::vector<BiquadCoeffs>& sections, double hz, double sampleRate)
{
    const double omega = twoPi * hz / sampleRate;
    std::complex<double> h { 1.0, 0.0 };

    for (const auto& s : sections)
        h *= s.responseAt (omega);

    const double mag = std::abs (h);
    return mag > 0.0 ? 20.0 * std::log10 (mag) : -1000.0;
}

/** Gaussian elimination with partial pivoting; `b` is overwritten with the
    solution. Returns false if the system is singular or the result is not
    finite, which the caller treats as "reject this step". */
bool solveInPlace (std::vector<double>& a, std::vector<double>& b, int n)
{
    for (int col = 0; col < n; ++col)
    {
        int pivot = col;

        for (int row = col + 1; row < n; ++row)
            if (std::abs (a[(size_t) (row * n + col)]) > std::abs (a[(size_t) (pivot * n + col)]))
                pivot = row;

        if (! (std::abs (a[(size_t) (pivot * n + col)]) > 0.0))
            return false;

        if (pivot != col)
        {
            for (int k = 0; k < n; ++k)
                std::swap (a[(size_t) (col * n + k)], a[(size_t) (pivot * n + k)]);

            std::swap (b[(size_t) col], b[(size_t) pivot]);
        }

        const double diagonal = a[(size_t) (col * n + col)];

        for (int row = col + 1; row < n; ++row)
        {
            const double factor = a[(size_t) (row * n + col)] / diagonal;

            if (factor == 0.0)
                continue;

            for (int k = col; k < n; ++k)
                a[(size_t) (row * n + k)] -= factor * a[(size_t) (col * n + k)];

            b[(size_t) row] -= factor * b[(size_t) col];
        }
    }

    for (int row = n - 1; row >= 0; --row)
    {
        double sum = b[(size_t) row];

        for (int k = row + 1; k < n; ++k)
            sum -= a[(size_t) (row * n + k)] * b[(size_t) k];

        b[(size_t) row] = sum / a[(size_t) (row * n + row)];
    }

    for (int k = 0; k < n; ++k)
        if (! std::isfinite (b[(size_t) k]))
            return false;

    return true;
}

//==============================================================================
// Pinking filter design.

constexpr double pinkLowestPoleHz = 4.0;
constexpr double pinkTopPoleFraction = 0.47;
constexpr double pinkPolesPerDecade = 4.0;
constexpr int pinkExtraPoles = 4;
constexpr int pinkFitPoints = 128;
constexpr int pinkMaxIterations = 25;
constexpr double pinkCostTolerance = 1.0e-3;
constexpr double pinkShelfMinDb = 0.01;
constexpr double pinkShelfMaxDb = 30.0;

/** Cascade of first-order shelves whose magnitude follows -3 dB/octave.

    Poles are fixed on a geometric grid; the shelf depths are the unknowns, fitted
    by Levenberg-Marquardt against the exact -3 dB/octave target with a free
    overall gain. See the class comment in the header for why the asymptotic
    placement rules are not good enough.
*/
std::vector<BiquadCoeffs> designPinkCascade (double sampleRate,
                                             double fitLowHz,
                                             double fitHighHz,
                                             double& worstDeviationDb)
{
    const double topPoleHz = pinkTopPoleFraction * sampleRate;
    const int numSections = pinkExtraPoles
                          + (int) std::lround (pinkPolesPerDecade * std::log10 (topPoleHz / pinkLowestPoleHz));
    const double ratio = std::pow (topPoleHz / pinkLowestPoleHz, 1.0 / (double) (numSections - 1));

    std::vector<double> poleW ((size_t) numSections);
    std::vector<double> shelfDb ((size_t) numSections, 10.0 * std::log10 (ratio));

    for (int k = 0; k < numSections; ++k)
        poleW[(size_t) k] = prewarp (pinkLowestPoleHz * std::pow (ratio, (double) k), sampleRate);

    auto build = [&] (const std::vector<double>& depths)
    {
        std::vector<BiquadCoeffs> sections ((size_t) numSections);

        for (int k = 0; k < numSections; ++k)
        {
            const double rho = std::pow (10.0, depths[(size_t) k] / 20.0);
            sections[(size_t) k] = firstOrderShelf (poleW[(size_t) k], rho * poleW[(size_t) k], sampleRate);
        }

        return sections;
    };

    // Fit grid, log spaced, with the ideal -3 dB/octave target on it.
    std::vector<double> grid ((size_t) pinkFitPoints), target ((size_t) pinkFitPoints);

    for (int i = 0; i < pinkFitPoints; ++i)
    {
        grid[(size_t) i] = fitLowHz * std::pow (fitHighHz / fitLowHz,
                                                (double) i / (double) (pinkFitPoints - 1));
        target[(size_t) i] = -10.0 * std::log10 (grid[(size_t) i]);
    }

    auto residualsFor = [&] (const std::vector<double>& depths, double gainDb)
    {
        const auto sections = build (depths);
        std::vector<double> res ((size_t) pinkFitPoints);

        for (int i = 0; i < pinkFitPoints; ++i)
            res[(size_t) i] = cascadeDb (sections, grid[(size_t) i], sampleRate)
                              + gainDb - target[(size_t) i];

        return res;
    };

    auto costOf = [] (const std::vector<double>& res)
    {
        double sum = 0.0;

        for (auto v : res)
            sum += v * v;

        return std::isfinite (sum) ? sum : 1.0e300;
    };

    const int unknowns = numSections + 1;   // shelf depths plus the free gain

    double gainDb = 0.0;
    auto residuals = residualsFor (shelfDb, gainDb);
    double cost = costOf (residuals);
    double lambda = 1.0e-6;

    std::vector<double> jacobian ((size_t) (pinkFitPoints * unknowns));
    std::vector<double> normal ((size_t) (unknowns * unknowns));
    std::vector<double> gradient ((size_t) unknowns);

    for (int iteration = 0; iteration < pinkMaxIterations; ++iteration)
    {
        const double previousCost = cost;
        const double step = 1.0e-3;

        for (int k = 0; k < numSections; ++k)
        {
            auto perturbed = shelfDb;
            perturbed[(size_t) k] += step;
            const auto other = residualsFor (perturbed, gainDb);

            for (int i = 0; i < pinkFitPoints; ++i)
                jacobian[(size_t) (i * unknowns + k)] =
                    (other[(size_t) i] - residuals[(size_t) i]) / step;
        }

        for (int i = 0; i < pinkFitPoints; ++i)
            jacobian[(size_t) (i * unknowns + numSections)] = 1.0;

        bool improved = false;

        for (int attempt = 0; attempt < 12 && ! improved; ++attempt)
        {
            std::fill (normal.begin(), normal.end(), 0.0);
            std::fill (gradient.begin(), gradient.end(), 0.0);

            for (int i = 0; i < pinkFitPoints; ++i)
                for (int a = 0; a < unknowns; ++a)
                {
                    const double ja = jacobian[(size_t) (i * unknowns + a)];
                    gradient[(size_t) a] -= ja * residuals[(size_t) i];

                    for (int b = 0; b < unknowns; ++b)
                        normal[(size_t) (a * unknowns + b)] += ja * jacobian[(size_t) (i * unknowns + b)];
                }

            for (int a = 0; a < unknowns; ++a)
                normal[(size_t) (a * unknowns + a)] *= 1.0 + lambda;

            if (! solveInPlace (normal, gradient, unknowns))
            {
                lambda *= 10.0;
                continue;
            }

            auto trial = shelfDb;

            for (int k = 0; k < numSections; ++k)
                trial[(size_t) k] = std::clamp (trial[(size_t) k] + gradient[(size_t) k],
                                                pinkShelfMinDb, pinkShelfMaxDb);

            const double trialGain = gainDb + gradient[(size_t) numSections];
            auto trialResiduals = residualsFor (trial, trialGain);
            const double trialCost = costOf (trialResiduals);

            if (trialCost < cost)
            {
                shelfDb = trial;
                gainDb = trialGain;
                residuals = trialResiduals;
                cost = trialCost;
                lambda = std::max (lambda * 0.1, 1.0e-12);
                improved = true;
            }
            else
            {
                lambda *= 10.0;
            }
        }

        // Stop when a whole iteration no longer buys anything, rather than
        // grinding through the iteration cap for a hundredth of a dB.
        if (! improved || previousCost - cost < pinkCostTolerance * cost)
            break;
    }

    const auto sections = build (shelfDb);

    // Report the peak-to-peak deviation from -3 dB/octave on a grid much finer
    // than the one that was fitted, so ripple between fit points is included.
    double highest = -1.0e300, lowest = 1.0e300;

    for (int i = 0; i <= 4000; ++i)
    {
        const double hz = fitLowHz * std::pow (fitHighHz / fitLowHz, (double) i / 4000.0);
        const double deviation = cascadeDb (sections, hz, sampleRate) + 10.0 * std::log10 (hz);

        highest = std::max (highest, deviation);
        lowest = std::min (lowest, deviation);
    }

    worstDeviationDb = highest - lowest;
    return sections;
}

/** RMS gain of a cascade for white input: sqrt of the mean of |H|^2 over
    [0, pi].

    Integrated on a logarithmic grid. A uniform grid would be hopeless here: a
    pink response has most of its power spread over the first few decades, so the
    integrand needs resolution measured in fractions of a Hz at the bottom and
    hundreds of Hz at the top. Below the lowest evaluated frequency the cascade
    is flat at unity, so that stretch is added analytically.
*/
double meanPowerGain (const std::vector<BiquadCoeffs>& sections, double sampleRate)
{
    const double omegaMin = twoPi * 0.01 / sampleRate;
    const int steps = 8192;                       // even, for Simpson
    const double logStep = std::log (pi / omegaMin) / (double) steps;

    auto integrand = [&] (double omega)
    {
        std::complex<double> h { 1.0, 0.0 };

        for (const auto& s : sections)
            h *= s.responseAt (omega);

        // d(omega) = omega * d(log omega)
        return std::norm (h) * omega;
    };

    double sum = integrand (omegaMin) + integrand (pi);

    for (int i = 1; i < steps; ++i)
        sum += (i % 2 == 1 ? 4.0 : 2.0) * integrand (omegaMin * std::exp (logStep * (double) i));

    const double integral = omegaMin + sum * logStep / 3.0;
    return std::sqrt (integral / pi);
}

} // namespace

//==============================================================================
WhiteNoise64::WhiteNoise64 (std::uint64_t seed) noexcept
{
    setSeed (seed);
}

void WhiteNoise64::setSeed (std::uint64_t seed) noexcept
{
    seedValue = seed;
    reset();
}

void WhiteNoise64::reset() noexcept
{
    // splitmix64 expansion, so that even a seed of zero gives a usable state.
    std::uint64_t z = seedValue;

    for (auto& s : state)
    {
        z += 0x9e3779b97f4a7c15ull;
        std::uint64_t x = z;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        s = x ^ (x >> 31);
    }
}

double WhiteNoise64::nextSample() noexcept
{
    // xoshiro256++: rotl(s0 + s3, 23) + s0.
    const std::uint64_t sum = state[0] + state[3];
    const std::uint64_t result = ((sum << 23) | (sum >> 41)) + state[0];

    const std::uint64_t t = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];
    state[2] ^= t;
    state[3] = (state[3] << 45) | (state[3] >> 19);

    // Top 53 bits to a double on [0, 1), then to [-1, 1).
    const double unit = (double) (result >> 11) * 0x1.0p-53;
    return 2.0 * unit - 1.0;
}

//==============================================================================
void PinkFilter64::prepare (double sampleRate)
{
    currentSampleRate = sampleRate;
    fitLowHz = 20.0;
    fitHighHz = std::min (20000.0, 0.45 * sampleRate);

    if (fitHighHz <= fitLowHz * 2.0)
        fitHighHz = fitLowHz * 2.0;

    sections = designPinkCascade (sampleRate, fitLowHz, fitHighHz, worstDeviationDb);
    noiseGain = meanPowerGain (sections, sampleRate);

    cascade.setCoeffs (sections);
    cascade.reset();
}

void PinkFilter64::reset() noexcept
{
    cascade.reset();
}

double PinkFilter64::responseDbAt (double hz) const
{
    return cascadeDb (sections, hz, currentSampleRate);
}

//==============================================================================
void PinkNoiseGenerator::prepare (double sampleRate)
{
    filter.prepare (sampleRate);
    updateScale();
    reset();
}

void PinkNoiseGenerator::reset()
{
    white.reset();
    filter.reset();
}

void PinkNoiseGenerator::setSeed (std::uint64_t seed) noexcept
{
    white.setSeed (seed);
}

void PinkNoiseGenerator::setRmsLevel (double rms) noexcept
{
    targetRms = rms;
    updateScale();
}

void PinkNoiseGenerator::updateScale() noexcept
{
    const double gain = WhiteNoise64::rms() * filter.getNoiseGain();
    scale = gain > 0.0 ? targetRms / gain : 0.0;
}

void PinkNoiseGenerator::process (double* out, int numSamples) noexcept
{
    for (int n = 0; n < numSamples; ++n)
        out[n] = nextSample();
}

//==============================================================================
void CalibrationNoiseGenerator::prepare (double sampleRate)
{
    noise.prepare (sampleRate);

    nickLengthSamples = samplesFor (calibration::nickSeconds, sampleRate);
    nickPeriodSamples = samplesFor (calibration::nickPeriodSeconds, sampleRate);

    if (nickLengthSamples >= nickPeriodSamples)
        nickLengthSamples = nickPeriodSamples - 1;

    reset();
}

void CalibrationNoiseGenerator::reset()
{
    noise.reset();
    phase = 0;
}

void CalibrationNoiseGenerator::setSeed (std::uint64_t seed) noexcept
{
    noise.setSeed (seed);
}

void CalibrationNoiseGenerator::process (double* out, int numSamples) noexcept
{
    const int gateFrom = nickPeriodSamples - nickLengthSamples;

    for (int n = 0; n < numSamples; ++n)
    {
        // The source runs through the nick; only the output is gated.
        const double sample = noise.nextSample();

        out[n] = phase < gateFrom ? sample : 0.0;

        if (++phase >= nickPeriodSamples)
            phase = 0;
    }
}

//==============================================================================
void CalibrationToneGenerator::prepare (double sampleRate)
{
    currentSampleRate = sampleRate;

    baseIncrement = twoPi * calibration::toneHz / sampleRate;
    modulatedIncrement = twoPi * calibration::toneModulatedHz / sampleRate;

    burstSamples = samplesFor (calibration::toneBurstSeconds, sampleRate);
    periodSamples = samplesFor (calibration::tonePeriodSeconds, sampleRate);

    if (burstSamples >= periodSamples)
        burstSamples = periodSamples - 1;

    reset();
}

void CalibrationToneGenerator::reset()
{
    angle = 0.0;
    phase = 0;
}

double CalibrationToneGenerator::getCurrentFrequencyHz() const noexcept
{
    return isModulated() ? calibration::toneModulatedHz : calibration::toneHz;
}

void CalibrationToneGenerator::process (double* out, int numSamples) noexcept
{
    const int burstFrom = periodSamples - burstSamples;

    for (int n = 0; n < numSamples; ++n)
    {
        // Amplitude is a constant and phase carries over from the previous
        // sample whatever the increment does, so the envelope is exactly flat
        // through both edges of the burst.
        out[n] = amplitude * std::sin (angle);

        angle += phase < burstFrom ? baseIncrement : modulatedIncrement;

        if (angle >= twoPi)
            angle -= twoPi;

        if (++phase >= periodSamples)
            phase = 0;
    }
}

//==============================================================================
namespace
{
    constexpr double detectorFastSeconds = 0.0007;
    constexpr double detectorSlowSeconds = 0.400;
    constexpr double detectorPredictorSeconds = 0.200;

    constexpr double detectorMinRms = 0.001;                 // -60 dBFS
    constexpr double detectorGapPowerRatio = 0.01;           // -20 dB in power
    constexpr double detectorNoiseLikenessDb = -55.0;

    constexpr double detectorMinGapSeconds = 0.008;
    constexpr double detectorMaxGapSeconds = 0.045;
    constexpr double detectorMinIntervalSeconds = 1.7;
    constexpr double detectorMaxIntervalSeconds = 2.3;
    constexpr double detectorTimeoutSeconds = 2.6;
    constexpr int detectorRequiredIntervals = 2;
}

void CalibrationNoiseDetector::prepare (double sampleRate)
{
    fastCoeff = smoothingCoeff (detectorFastSeconds, sampleRate);
    slowCoeff = smoothingCoeff (detectorSlowSeconds, sampleRate);
    predictorCoeff = smoothingCoeff (detectorPredictorSeconds, sampleRate);

    minGapSamples = samplesFor (detectorMinGapSeconds, sampleRate);
    maxGapSamples = samplesFor (detectorMaxGapSeconds, sampleRate);
    minIntervalSamples = samplesFor (detectorMinIntervalSeconds, sampleRate);
    maxIntervalSamples = samplesFor (detectorMaxIntervalSeconds, sampleRate);
    timeoutSamples = samplesFor (detectorTimeoutSeconds, sampleRate);

    reset();
}

void CalibrationNoiseDetector::reset()
{
    fastPower = slowPower = 0.0;
    midPower = sumPower = crossPower = 0.0;
    history1 = history2 = 0.0;

    position = 0;
    gapStart = 0;
    lastNickStart = -1;

    nickCount = 0;
    validIntervals = 0;
    inGap = false;
    present = false;
}

double CalibrationNoiseDetector::getLevelDbFs() const noexcept
{
    return slowPower > 0.0 ? 10.0 * std::log10 (slowPower) : -1000.0;
}

double CalibrationNoiseDetector::getNoiseLikenessDb() const noexcept
{
    if (sumPower <= 0.0 || midPower <= 0.0)
        return -1000.0;

    // Residual of the least-squares fit of x[n+1] + x[n-1] = c * x[n]. Exactly
    // zero for any sinusoid, whatever its frequency and amplitude.
    const double residual = sumPower - crossPower * crossPower / midPower;

    if (residual <= 0.0)
        return -1000.0;

    return 10.0 * std::log10 (residual / sumPower);
}

void CalibrationNoiseDetector::push (const double* samples, int numSamples) noexcept
{
    const double minPower = detectorMinRms * detectorMinRms;

    for (int n = 0; n < numSamples; ++n)
    {
        const double x = samples[n];
        const double power = x * x;

        fastPower += fastCoeff * (power - fastPower);
        slowPower += slowCoeff * (power - slowPower);

        // Sinusoid predictor, centred one sample back so the triple is complete.
        const double predictorSum = x + history2;
        const double predictorMid = history1;

        midPower += predictorCoeff * (predictorMid * predictorMid - midPower);
        sumPower += predictorCoeff * (predictorSum * predictorSum - sumPower);
        crossPower += predictorCoeff * (predictorMid * predictorSum - crossPower);

        history2 = history1;
        history1 = x;

        ++position;

        const bool levelOk = slowPower > minPower;

        if (! levelOk)
        {
            // Silence, or near enough. Forget everything rather than let a gap
            // that opened before the signal died count towards an interval.
            inGap = false;
            lastNickStart = -1;
            validIntervals = 0;
            present = false;
            continue;
        }

        const bool below = fastPower < detectorGapPowerRatio * slowPower;

        if (below && ! inGap)
        {
            inGap = true;
            gapStart = position;
        }
        else if (! below && inGap)
        {
            inGap = false;

            const long long length = position - gapStart;
            const bool noiseLike = getNoiseLikenessDb() >= detectorNoiseLikenessDb;

            if (length >= minGapSamples && length <= maxGapSamples && noiseLike)
            {
                ++nickCount;

                if (lastNickStart >= 0)
                {
                    const long long interval = gapStart - lastNickStart;

                    if (interval >= minIntervalSamples && interval <= maxIntervalSamples)
                        ++validIntervals;
                    else
                        validIntervals = 0;
                }

                lastNickStart = gapStart;
            }
        }

        if (lastNickStart >= 0 && position - lastNickStart > timeoutSamples)
        {
            // A nick was expected and did not arrive.
            lastNickStart = -1;
            validIntervals = 0;
        }

        present = validIntervals >= detectorRequiredIntervals;
    }
}

//==============================================================================
void AutoCompare::prepare (double sampleRate)
{
    currentSampleRate = sampleRate;
    segmentSamples = samplesFor (calibration::autoCompareSegmentSeconds, sampleRate);
    reset();
}

void AutoCompare::reset()
{
    phase.store (0, std::memory_order_relaxed);
    running.store (false, std::memory_order_relaxed);
}

void AutoCompare::setRunning (bool shouldRun) noexcept
{
    if (shouldRun == running.load (std::memory_order_relaxed))
        return;

    running.store (shouldRun, std::memory_order_relaxed);
    phase.store (0, std::memory_order_relaxed);
}

AutoCompareSource AutoCompare::getCurrentSource() const noexcept
{
    if (! running.load (std::memory_order_relaxed))
        return AutoCompareSource::tape;

    return phase.load (std::memory_order_relaxed) < segmentSamples
             ? AutoCompareSource::tape : AutoCompareSource::reference;
}

double AutoCompare::getSecondsIntoSegment() const noexcept
{
    const int p = phase.load (std::memory_order_relaxed);
    const int intoSegment = p < segmentSamples ? p : p - segmentSamples;
    return (double) intoSegment / currentSampleRate;
}

void AutoCompare::process (const double* const* tape, const double* reference,
                           double* const* out, int numChannels, int numSamples) noexcept
{
    if (! running.load (std::memory_order_relaxed))
    {
        for (int ch = 0; ch < numChannels; ++ch)
            for (int n = 0; n < numSamples; ++n)
                out[ch][n] = tape[ch][n];

        return;
    }

    // The clock advances once per sample whatever the channel count, so every
    // channel switches source on the same sample and the segments keep their
    // configured length in stereo. A local copy, stored once at the end: the
    // member is atomic only for the front-panel indicator's benefit.
    int p = phase.load (std::memory_order_relaxed);

    for (int n = 0; n < numSamples; ++n)
    {
        const bool onTape = p < segmentSamples;

        for (int ch = 0; ch < numChannels; ++ch)
            out[ch][n] = onTape ? tape[ch][n] : reference[n];

        if (++p >= 2 * segmentSamples)
            p = 0;
    }

    phase.store (p, std::memory_order_relaxed);
}

void AutoCompare::process (const double* tape, const double* reference, double* out, int numSamples) noexcept
{
    const double* tapeChannels[1] = { tape };
    double* outChannels[1] = { out };

    process (tapeChannels, reference, outChannels, 1, numSamples);
}

} // namespace lime
