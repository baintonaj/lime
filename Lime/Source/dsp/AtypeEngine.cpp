/*
  ==============================================================================

    AtypeEngine — band design, low-level gain fit, compressors, exact decoder.

    Three things here are worth reading before the code.

    The band split. Bands 1, 3 and 4 are 12 dB/octave Butterworth sections at the
    published corners. Band 2 is not designed at all: it is formed as
    1 - band1 - band3, which is the literal reading of "complementary to bands 1
    and 3" and makes those three sum to unity in magnitude and phase, at every
    frequency and every sample rate, with no ripple to fit around.

    The low-level gain. Band 4 overlaps band 3 to give the rise above 5 kHz, and
    a 12 dB/octave high-pass is 180 degrees out below its corner, so its
    contribution subtracts there. That is physical, not a modelling artefact, and
    it means no set of band gains makes the total both exactly flat and exactly
    5 dB up at 15 kHz. The four gains are therefore fitted, weighted so the flat
    region and the 15 kHz endpoint — the two things the paper actually states —
    dominate the shape in between.

    The decoder. Given the compressors' gains, which depend only on state from
    earlier samples, the encoder's map from the current sample v to its output is

        E(v) = v + sum_i g_i * clamp (a_i v + b_i, -clip, +clip)

    with every a_i >= 0 — a rate-dependent property of the band designs that
    holds across all supported rates and is asserted where the slopes are
    formed. That makes E continuous, piecewise linear and strictly
    increasing, so it inverts in closed form: try the unclipped segment, and if
    the answer does not land there, enumerate the eight breakpoints and solve the
    segment that brackets the input. No iteration to a tolerance and no one-sample
    delay in the feedback path, so the decoder undoes the encoder to numerical
    precision.

  ==============================================================================
*/

#include "AtypeEngine.h"

#include "DspMath.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace lime
{

namespace
{

//==============================================================================
/** Bilinear design of a two-pole section from a normalised analog prototype.

    `K` is the prewarped corner, tan(pi * fc / fs), so the -3 dB point lands on
    fc. The corner is clamped just below Nyquist: at 22.05 kHz the 9 kHz band is
    still legal but a hypothetical rate below 18 kHz would not be.
*/
double prewarpedCorner (double cornerHz, double sampleRate)
{
    constexpr double pi = 3.14159265358979323846;

    const double nyquist = 0.5 * sampleRate;
    const double clamped = std::clamp (cornerHz, 1.0e-4 * nyquist, 0.49 * sampleRate);

    return std::tan (pi * clamped / sampleRate);
}

double onePoleCoefficient (double seconds, double sampleRate) noexcept
{
    if (seconds <= 0.0 || sampleRate <= 0.0)
        return 1.0;

    return 1.0 - std::exp (-1.0 / (seconds * sampleRate));
}

//==============================================================================
// Fit machinery for the four band gains.

using FitPoint = std::array<double, 4>;

constexpr size_t fitDimension = 4;
constexpr size_t fitPoints = 129;
constexpr double fitLowestHz = 20.0;
constexpr double fitTransitionWeight = 0.3;   // "rising smoothly" is not a number
constexpr double fitTopWeightFromHz = 14000.0;
constexpr double fitGainFloor = 0.05;
constexpr double fitGainCeiling = 8.0;

FitPoint clampGains (const FitPoint& p)
{
    FitPoint out {};

    for (size_t d = 0; d < fitDimension; ++d)
        out[d] = std::clamp (p[d], fitGainFloor, fitGainCeiling);

    return out;
}

/** Nelder-Mead in four dimensions: fixed schedule, no randomness, no derivative,
    so the same sample rate always yields the same gains. */
template <typename Cost>
FitPoint minimise (FitPoint start, double step, Cost&& cost)
{
    constexpr size_t vertices = fitDimension + 1;
    constexpr int maxIterations = 4000;

    std::array<FitPoint, vertices> simplex {};
    std::array<double, vertices> value {};

    for (size_t i = 0; i < vertices; ++i)
    {
        simplex[i] = start;

        if (i > 0)
            simplex[i][i - 1] += step;

        value[i] = cost (simplex[i]);
    }

    for (int iteration = 0; iteration < maxIterations; ++iteration)
    {
        size_t best = 0, worst = 0;

        for (size_t i = 0; i < vertices; ++i)
        {
            if (value[i] < value[best])  best = i;
            if (value[i] > value[worst]) worst = i;
        }

        if (value[worst] - value[best] <= 1.0e-15)
            break;

        FitPoint centroid {};

        for (size_t i = 0; i < vertices; ++i)
            if (i != worst)
                for (size_t d = 0; d < fitDimension; ++d)
                    centroid[d] += simplex[i][d] / (double) fitDimension;

        // factor 1 reflects the worst vertex through the centroid, 2 expands past
        // it, -0.5 contracts halfway back towards it.
        const auto along = [&] (double factor)
        {
            FitPoint p {};

            for (size_t d = 0; d < fitDimension; ++d)
                p[d] = centroid[d] + factor * (centroid[d] - simplex[worst][d]);

            return p;
        };

        const FitPoint reflected = along (1.0);
        const double reflectedValue = cost (reflected);

        if (reflectedValue < value[best])
        {
            const FitPoint expanded = along (2.0);
            const double expandedValue = cost (expanded);

            simplex[worst] = (expandedValue < reflectedValue) ? expanded : reflected;
            value[worst] = std::min (expandedValue, reflectedValue);
            continue;
        }

        bool beatsSomeVertex = false;

        for (size_t i = 0; i < vertices; ++i)
            if (i != worst && reflectedValue < value[i])
                beatsSomeVertex = true;

        if (beatsSomeVertex)
        {
            simplex[worst] = reflected;
            value[worst] = reflectedValue;
            continue;
        }

        const FitPoint contracted = along (-0.5);
        const double contractedValue = cost (contracted);

        if (contractedValue < value[worst])
        {
            simplex[worst] = contracted;
            value[worst] = contractedValue;
            continue;
        }

        // Nothing helped: shrink the whole simplex towards the best vertex.
        for (size_t i = 0; i < vertices; ++i)
        {
            if (i == best)
                continue;

            for (size_t d = 0; d < fitDimension; ++d)
                simplex[i][d] = simplex[best][d] + 0.5 * (simplex[i][d] - simplex[best][d]);

            value[i] = cost (simplex[i]);
        }
    }

    size_t best = 0;

    for (size_t i = 0; i < vertices; ++i)
        if (value[i] < value[best])
            best = i;

    return simplex[best];
}

} // namespace

//==============================================================================
BiquadCoeffs atypeBands::lowPass12 (double cornerHz, double sampleRate)
{
    const double k = prewarpedCorner (cornerHz, sampleRate);
    const double kk = k * k;
    const double a0 = 1.0 + k / butterworthQ + kk;

    BiquadCoeffs c;

    c.b0 = kk / a0;
    c.b1 = 2.0 * kk / a0;
    c.b2 = kk / a0;
    c.a1 = 2.0 * (kk - 1.0) / a0;
    c.a2 = (1.0 - k / butterworthQ + kk) / a0;

    return c;
}

BiquadCoeffs atypeBands::highPass12 (double cornerHz, double sampleRate)
{
    const double k = prewarpedCorner (cornerHz, sampleRate);
    const double kk = k * k;
    const double a0 = 1.0 + k / butterworthQ + kk;

    BiquadCoeffs c;

    c.b0 = 1.0 / a0;
    c.b1 = -2.0 / a0;
    c.b2 = 1.0 / a0;
    c.a1 = 2.0 * (kk - 1.0) / a0;
    c.a2 = (1.0 - k / butterworthQ + kk) / a0;

    return c;
}

double atypeBands::targetGainDb (double hz)
{
    if (hz <= flatUpToHz)
        return flatGainDb;

    if (hz >= topHz)
        return topGainDb;

    const double span = std::log (topHz / flatUpToHz);

    return flatGainDb + (topGainDb - flatGainDb) * std::log (hz / flatUpToHz) / span;
}

std::array<double, 4> atypeBands::solveBandGains (double sampleRate,
                                                  const BiquadCoeffs& band1,
                                                  const BiquadCoeffs& band3,
                                                  const BiquadCoeffs& band4)
{
    constexpr double twoPi = 6.28318530717958647692;

    // The grid stops at 0.45 of the sample rate: above that the bilinear design
    // no longer tracks the analog prototype, and at low rates the top of the
    // published curve simply does not exist.
    const double highest = std::min (topHz, 0.45 * sampleRate);
    const double lowest = std::min (fitLowestHz, 0.5 * highest);

    std::vector<std::complex<double>> r1 (fitPoints), r2 (fitPoints), r3 (fitPoints), r4 (fitPoints);
    std::vector<double> target (fitPoints), weight (fitPoints);

    for (size_t i = 0; i < fitPoints; ++i)
    {
        const double hz = lowest * std::pow (highest / lowest, (double) i / (double) (fitPoints - 1));
        const double omega = twoPi * hz / sampleRate;

        r1[i] = band1.responseAt (omega);
        r3[i] = band3.responseAt (omega);
        r4[i] = band4.responseAt (omega);
        r2[i] = 1.0 - r1[i] - r3[i];

        target[i] = targetGainDb (hz);

        // The paper pins the flat region and the 15 kHz endpoint; "rising
        // smoothly" between them is not a number, so the interpolated middle is
        // weighted down rather than allowed to pull the two stated ends about.
        weight[i] = (hz <= flatUpToHz || hz >= fitTopWeightFromHz) ? 1.0 : fitTransitionWeight;
    }

    const auto cost = [&] (const FitPoint& raw)
    {
        const FitPoint g = clampGains (raw);
        double sum = 0.0;

        for (size_t i = 0; i < fitPoints; ++i)
        {
            const std::complex<double> total = 1.0 + g[0] * r1[i] + g[1] * r2[i]
                                                   + g[2] * r3[i] + g[3] * r4[i];

            const double db = 20.0 * std::log10 (std::max (1.0e-300, std::abs (total)));
            const double error = db - target[i];

            sum += weight[i] * error * error;
        }

        return sum;
    };

    // Starting point: the gain that would give exactly flatGainDb if the four
    // bands summed in phase to unity, which bands 1 to 3 do by construction.
    const double nominal = std::pow (10.0, flatGainDb / 20.0) - 1.0;

    return clampGains (minimise ({ nominal, nominal, nominal, nominal }, 0.3, cost));
}

//==============================================================================
struct AtypeEngine::Channel
{
    Biquad64 band1, band3, band4;
    std::array<BandState, numBands> state {};
};

//==============================================================================
AtypeEngine::AtypeEngine() = default;
AtypeEngine::~AtypeEngine() = default;

void AtypeEngine::setPeakOperatingLevel (double amplitude) noexcept
{
    if (amplitude > 0.0)
        peakLevel = amplitude;
}

void AtypeEngine::setSidechainHighPass (double cornerHz) noexcept
{
    // A static control costs nothing: four magnitudes rebuild only on a move.
    if (cornerHz == detectorHighPassHz)
        return;

    detectorHighPassHz = cornerHz;
    rebuildDetectorWeights();
}

void AtypeEngine::setSidechainLowPass (double cornerHz) noexcept
{
    if (cornerHz == detectorLowPassHz)
        return;

    detectorLowPassHz = cornerHz;
    rebuildDetectorWeights();
}

void AtypeEngine::setTimeConstants (double upAttack, double upRelease,
                                    double downAttack, double downRelease) noexcept
{
    if (upAttack == upAttackMs && upRelease == upReleaseMs
        && downAttack == downAttackMs && downRelease == downReleaseMs)
        return;

    upAttackMs = upAttack;
    upReleaseMs = upRelease;
    downAttackMs = downAttack;
    downReleaseMs = downRelease;

    rebuildEffectiveTimes();
}

void AtypeEngine::rebuildEffectiveTimes() noexcept
{
    const double attackMs[2] = { upAttackMs, downAttackMs };
    const double releaseMs[2] = { upReleaseMs, downReleaseMs };

    for (int d = 0; d < 2; ++d)
    {
        const auto i = (size_t) d;

        effectiveAttackRate[i] = attackMs[d] > 0.0 && currentSampleRate > 0.0
                                   ? 1.0 / (attackMs[d] * 0.001 * currentSampleRate)
                                   : attackRate;

        effectiveRecoveryCoef[i] = releaseMs[d] > 0.0 && currentSampleRate > 0.0
                                     ? onePoleCoefficient (releaseMs[d] * 0.001, currentSampleRate)
                                     : recoveryCoef;
    }
}

void AtypeEngine::rebuildDetectorWeights() noexcept
{
    // Each band is represented by the geometric centre of its nominal range,
    // the audio band's own edges standing in where a band runs out to one: the
    // 80 Hz low-pass from 20 Hz, bands 3 and 4 up to 20 kHz. A band-pass has no
    // single frequency, but its log midpoint is where its energy is judged
    // from, and four scalars is what a four-band side chain can carry. The two
    // corners multiply, as they do per bin in the spectral side chains.
    const double centres[numBands] = {
        std::sqrt (20.0 * atypeBands::band1CornerHz),
        std::sqrt (atypeBands::band1CornerHz * atypeBands::band3CornerHz),
        std::sqrt (atypeBands::band3CornerHz * 20000.0),
        std::sqrt (atypeBands::band4CornerHz * 20000.0)
    };

    for (int b = 0; b < numBands; ++b)
    {
        double weight = 1.0;

        if (detectorHighPassHz > 0.0)
            weight *= butterworth3Magnitude (centres[b], detectorHighPassHz, true);

        if (detectorLowPassHz > 0.0)
            weight *= butterworth3Magnitude (centres[b], detectorLowPassHz, false);

        detectorWeight[(size_t) b] = weight;
    }
}

void AtypeEngine::prepare (double sampleRate, int maxBlockSize, int numChannels)
{
    (void) maxBlockSize;   // no scratch: every path is a sample recursion

    currentSampleRate = (sampleRate > 0.0) ? sampleRate : 48000.0;

    bandCoeffs[0] = atypeBands::lowPass12  (atypeBands::band1CornerHz, currentSampleRate);
    bandCoeffs[1] = atypeBands::highPass12 (atypeBands::band3CornerHz, currentSampleRate);
    bandCoeffs[2] = atypeBands::highPass12 (atypeBands::band4CornerHz, currentSampleRate);

    gains = atypeBands::solveBandGains (currentSampleRate, bandCoeffs[0], bandCoeffs[1], bandCoeffs[2]);

    threshold = peakLevel * std::pow (10.0, atypeBands::thresholdDb / 20.0);

    // Clipper level. A peak-amplitude step drives every band past its clipper at
    // once and in phase, so the worst differential component the output can see
    // is clipLevel * sum(gains); setting that to the permitted 2 dB of overshoot
    // is the paper's own criterion, taken at its worst case.
    const double sumGains = gains[0] + gains[1] + gains[2] + gains[3];
    const double permitted = std::pow (10.0, atypeBands::overshootLimitDb / 20.0) - 1.0;

    clipLevel = (sumGains > 0.0) ? peakLevel * permitted / sumGains : peakLevel;

    peakAttackCoef = onePoleCoefficient (atypeBands::peakAttackSeconds, currentSampleRate);
    peakReleaseCoef = onePoleCoefficient (atypeBands::peakReleaseSeconds, currentSampleRate);
    averageCoef = onePoleCoefficient (atypeBands::averageSeconds, currentSampleRate);

    // The slow attack is expressed as a rate rather than a coefficient because it
    // is scaled by the size of the transition on every sample; the fast limit is
    // a genuine one-pole coefficient and caps the result.
    attackRate = 1.0 / (atypeBands::attackSeconds * currentSampleRate);
    attackCoefMax = onePoleCoefficient (atypeBands::fastestAttackSeconds, currentSampleRate);
    recoveryCoef = onePoleCoefficient (atypeBands::recoverySeconds, currentSampleRate);

    // User times set before prepare — or surviving a sample-rate change — are
    // honoured rather than lost, the setSidechainHighPass contract.
    rebuildEffectiveTimes();

    channels.clear();
    channels.reserve ((size_t) std::max (0, numChannels));

    for (int ch = 0; ch < numChannels; ++ch)
    {
        auto c = std::make_unique<Channel>();

        c->band1.setCoeffs (bandCoeffs[0]);
        c->band3.setCoeffs (bandCoeffs[1]);
        c->band4.setCoeffs (bandCoeffs[2]);

        channels.push_back (std::move (c));
    }

    prepared = true;

    reset();
}

void AtypeEngine::reset()
{
    for (auto& c : channels)
    {
        c->band1.reset();
        c->band3.reset();
        c->band4.reset();

        // Quiescent state is full noise reduction: the compressors sit at unity
        // gain, which is where a decoder must start if it is to track an encoder
        // that started there too.
        for (auto& s : c->state)
            s = BandState {};
    }
}

//==============================================================================
const BiquadCoeffs& AtypeEngine::getBandCoeffs (int band) const noexcept
{
    // Band 2 has no section of its own, so anything other than band 3 or band 4
    // reads band 1's.
    const size_t index = (band == 3) ? 2 : ((band == 2) ? 1 : 0);

    return bandCoeffs[index];
}

std::complex<double> AtypeEngine::bandResponseAt (int band, double hz) const
{
    constexpr double twoPi = 6.28318530717958647692;

    const double omega = twoPi * hz / currentSampleRate;

    const std::complex<double> low = bandCoeffs[0].responseAt (omega);
    const std::complex<double> high = bandCoeffs[1].responseAt (omega);

    switch (band)
    {
        case 0:  return low;
        case 1:  return 1.0 - low - high;      // complementary to bands 1 and 3
        case 2:  return high;
        default: break;
    }

    return bandCoeffs[2].responseAt (omega);
}

std::complex<double> AtypeEngine::lowLevelResponseAt (double hz) const
{
    std::complex<double> total { 1.0, 0.0 };

    for (int b = 0; b < numBands; ++b)
        total += gains[(size_t) b] * bandResponseAt (b, hz);

    return total;
}

double AtypeEngine::lowLevelGainDb (double hz) const
{
    return 20.0 * std::log10 (std::max (1.0e-300, std::abs (lowLevelResponseAt (hz))));
}

//==============================================================================
void AtypeEngine::updateControl (BandState& s, double bandSample) const noexcept
{
    const double rectified = std::fabs (bandSample);

    // Peak and average in parallel. The paper asks for a control signal that
    // tracks a combination of the two, as a proxy for rms that stays correct when
    // the channel's phase response is badly non-linear.
    s.peak += (rectified - s.peak) * (rectified > s.peak ? peakAttackCoef : peakReleaseCoef);
    s.average += (rectified * atypeBands::averageScale - s.average) * averageCoef;

    const double level = atypeBands::peakWeight * s.peak
                         + (1.0 - atypeBands::peakWeight) * s.average;

    // Linear limiter: unity below the threshold, and above it holding the band's
    // contribution constant at the threshold, so the differential component
    // becomes negligible in proportion as the signal grows.
    const double targetGain = (level > threshold) ? threshold / level : 1.0;

    // The direction selects the user's ballistics: the up pair while encoding,
    // the down pair while decoding, the published constants wherever a time is
    // disengaged.
    const auto dir = (size_t) (mode == AtypeMode::decode ? 1 : 0);

    if (targetGain < s.gain)
    {
        // Attack. Slow for small amplitude variations, shortened in proportion to
        // the size of the transition, and capped at the fast limit; steps big
        // enough to reach the clipper are handled by the clipper meanwhile.
        const double ratio = s.gain / targetGain;

        s.gain += (targetGain - s.gain) * std::min (ratio * effectiveAttackRate[dir], attackCoefMax);
    }
    else
    {
        s.gain += (targetGain - s.gain) * effectiveRecoveryCoef[dir];
    }
}

//==============================================================================
double AtypeEngine::solveForInput (Channel& c, double y) const noexcept
{
    // Affine form of each band in the unknown sample v. A transposed direct form
    // II section emits b0 * x + s1, so the slope is the section's b0 and the
    // offset is what it would emit for a zero input, read from the state without
    // advancing it. Band 2 inherits both from the identity
    // band2 = 1 - band1 - band3.
    double slope[numBands], offset[numBands];

    slope[0] = c.band1.getCoeffs().b0;
    slope[2] = c.band3.getCoeffs().b0;
    slope[3] = c.band4.getCoeffs().b0;
    slope[1] = 1.0 - slope[0] - slope[2];

    // The closed-form inversion needs every slope non-negative so that E stays
    // increasing. The Butterworth b0's are positive by construction; band 2's
    // slope is the rate-dependent one, and stays positive at every supported
    // rate — about 0.90 at 8 kHz falling to 0.067 at 192 kHz, since the
    // high-pass b0 approaches 1 as the 3 kHz corner shrinks against Nyquist.
    assert (slope[1] >= 0.0);

    offset[0] = c.band1.outputForZeroInput();
    offset[2] = c.band3.outputForZeroInput();
    offset[3] = c.band4.outputForZeroInput();
    offset[1] = -offset[0] - offset[2];

    // Limiter output before clipping: u_i(v) = a_i v + b_i, with a_i >= 0.
    const double* g = gains.data();
    double a[numBands], b[numBands];

    for (int i = 0; i < numBands; ++i)
    {
        a[i] = c.state[(size_t) i].gain * slope[i];
        b[i] = c.state[(size_t) i].gain * offset[i];
    }

    // The unclipped segment first. That is the steady state at every level: the
    // clipper only acts on transients large enough to outrun the slow attack.
    double totalSlope = 1.0, totalOffset = 0.0;

    for (int i = 0; i < numBands; ++i)
    {
        totalSlope += g[i] * a[i];
        totalOffset += g[i] * b[i];
    }

    double v = (y - totalOffset) / totalSlope;
    bool clipping = false;

    for (int i = 0; i < numBands; ++i)
        if (std::fabs (a[i] * v + b[i]) > clipLevel)
            clipping = true;

    if (! clipping)
        return v;

    // Otherwise enumerate the segment boundaries and solve the one that brackets
    // y. E is strictly increasing, so exactly one segment does.
    double breakpoint[2 * numBands];
    int count = 0;

    for (int i = 0; i < numBands; ++i)
    {
        if (a[i] <= 0.0)
            continue;

        breakpoint[count++] = ( clipLevel - b[i]) / a[i];
        breakpoint[count++] = (-clipLevel - b[i]) / a[i];
    }

    if (count == 0)
        return v;   // no band responds to the current sample: nothing to bracket

    for (int i = 1; i < count; ++i)
    {
        const double x = breakpoint[i];
        int j = i - 1;

        while (j >= 0 && breakpoint[j] > x)
        {
            breakpoint[j + 1] = breakpoint[j];
            --j;
        }

        breakpoint[j + 1] = x;
    }

    const auto encodeAt = [&] (double x)
    {
        double e = x;

        for (int i = 0; i < numBands; ++i)
            e += g[i] * std::clamp (a[i] * x + b[i], -clipLevel, clipLevel);

        return e;
    };

    // Solves E(v) = y using the affine form that holds wherever `probe` sits.
    // Exact, provided the segment containing `probe` is the one containing v.
    const auto solveNear = [&] (double probe)
    {
        double s = 1.0, o = 0.0;

        for (int i = 0; i < numBands; ++i)
        {
            const double u = a[i] * probe + b[i];

            if (u > clipLevel)
            {
                o += g[i] * clipLevel;
            }
            else if (u < -clipLevel)
            {
                o -= g[i] * clipLevel;
            }
            else
            {
                s += g[i] * a[i];
                o += g[i] * b[i];
            }
        }

        return (y - o) / s;
    };

    if (y <= encodeAt (breakpoint[0]))
        return solveNear (breakpoint[0] - 1.0);

    for (int i = 0; i + 1 < count; ++i)
        if (y <= encodeAt (breakpoint[i + 1]))
            return solveNear (0.5 * (breakpoint[i] + breakpoint[i + 1]));

    return solveNear (breakpoint[count - 1] + 1.0);
}

//==============================================================================
void AtypeEngine::process (double* const* channelData, int numChannels, int numSamples)
{
    if (! prepared || numSamples <= 0)
        return;

    const int usable = std::min (numChannels, (int) channels.size());
    const bool decoding = (mode == AtypeMode::decode);

    for (int ch = 0; ch < usable; ++ch)
    {
        auto& c = *channels[(size_t) ch];
        double* io = channelData[ch];

        if (io == nullptr)
            continue;

        for (int n = 0; n < numSamples; ++n)
        {
            // The sample the encoder saw. On encode that is the input; on decode
            // it is recovered by inverting the encoder's map, which is what makes
            // the feedback configuration exactly complementary.
            const double v = decoding ? solveForInput (c, io[n]) : io[n];

            double band[numBands];

            band[0] = c.band1.processSample (v);
            band[2] = c.band3.processSample (v);
            band[3] = c.band4.processSample (v);
            band[1] = v - band[0] - band[2];

            // Each band: linear limiter, then the symmetrical clipper, then the
            // band's contribution to the differential component. The decoder
            // already subtracted that component in recovering v, so only the
            // encoder sums it here; both still run the bands, whose state the
            // control update below feeds on.
            if (decoding)
            {
                io[n] = v;
            }
            else
            {
                double differential = 0.0;

                for (int b = 0; b < numBands; ++b)
                    differential += gains[(size_t) b]
                                    * std::clamp (c.state[(size_t) b].gain * band[b], -clipLevel, clipLevel);

                io[n] = v + differential;
            }

            // The detectors read through the side-chain high pass; the band
            // signals themselves — and the differential built from them — do
            // not. Encoder and decoder scale identically, so the states they
            // derive from the recovered v stay in lockstep.
            for (int b = 0; b < numBands; ++b)
                updateControl (c.state[(size_t) b], band[b] * detectorWeight[(size_t) b]);
        }
    }

    // Channels beyond those prepared are left untouched rather than zeroed.
}

} // namespace lime
