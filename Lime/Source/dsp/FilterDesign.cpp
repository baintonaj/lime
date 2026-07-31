/*
  ==============================================================================

    FilterDesign — analog prototypes and their bilinear mapping to biquads.

  ==============================================================================
*/

#include "FilterDesign.h"

#include <cmath>
#include <map>
#include <mutex>

namespace lime
{

namespace
{
    constexpr double pi = 3.141592653589793238462643383280;
    constexpr double twoPi = 2.0 * pi;
    const double sqrt2 = std::sqrt (2.0);

    double dbToGain (double db) { return std::pow (10.0, db / 20.0); }
}

//==============================================================================
ShelfPrototype ShelfPrototype::lowPassShelf (double cornerHz, double shelfDb)
{
    // Attenuation g at high frequency. Placing the zeros a factor 1/sqrt(g) above
    // the poles makes H(inf) = (wp/wz)^2 = g, while K = g restores unity at DC and
    // leaves the -3 dB point exactly at the pole frequency.
    const double g = dbToGain (-std::abs (shelfDb));

    ShelfPrototype p;
    p.poleHz = cornerHz;
    p.zeroHz = cornerHz / std::sqrt (g);
    p.gain = g;

    return p;
}

ShelfPrototype ShelfPrototype::highPassShelf (double cornerHz, double shelfDb)
{
    // Mirror image: zeros a factor sqrt(g) below the poles gives H(0) = g, and
    // unity at high frequency needs no scaling.
    const double g = dbToGain (-std::abs (shelfDb));

    ShelfPrototype p;
    p.poleHz = cornerHz;
    p.zeroHz = cornerHz * std::sqrt (g);
    p.gain = 1.0;

    return p;
}

std::complex<double> ShelfPrototype::responseAt (double hz) const
{
    const std::complex<double> s { 0.0, twoPi * hz };
    const double wz = twoPi * zeroHz;
    const double wp = twoPi * poleHz;

    const auto num = s * s + sqrt2 * wz * s + wz * wz;
    const auto den = s * s + sqrt2 * wp * s + wp * wp;

    return gain * num / den;
}

BiquadCoeffs ShelfPrototype::toBiquad (double sampleRate) const
{
    const double c = 2.0 * sampleRate;
    const double cc = c * c;

    // Prewarp so the intended frequency in Hz survives the bilinear map.
    const auto prewarp = [c] (double hz) { return c * std::tan (twoPi * hz / c); };

    const double wp = prewarp (poleHz);

    const double a0 = cc + sqrt2 * wp * c + wp * wp;
    const double a1 = 2.0 * wp * wp - 2.0 * cc;
    const double a2 = cc - sqrt2 * wp * c + wp * wp;

    // The analog DC gain, H(0) = K * wz^2 / wp^2: exactly 1 for the low-pass
    // shelf and exactly the shelf attenuation for the high-pass, both by
    // construction of the prototypes above.
    const double dcGain = gain * (zeroHz * zeroHz) / (poleHz * poleHz);

    // Assembles the numerator for a given digital zero frequency, then pins the
    // gain at z = 1 to the analog DC gain. Prewarping preserves each critical
    // frequency but not their ratio, so the raw bilinear numerator lands only
    // near the analog DC value (over half a dB off for the fitted HF
    // antisaturation shelf at 48 kHz); scaling by desired * (1 + a1 + a2) /
    // (b0 + b1 + b2) makes DC exact whichever way the zero was placed.
    const auto makeCoeffs = [&] (double wz)
    {
        BiquadCoeffs out;
        const double norm = 1.0 / a0;

        out.b0 = (cc + sqrt2 * wz * c + wz * wz) * norm;
        out.b1 = (2.0 * wz * wz - 2.0 * cc) * norm;
        out.b2 = (cc - sqrt2 * wz * c + wz * wz) * norm;
        out.a1 = a1 * norm;
        out.a2 = a2 * norm;

        const double scale = dcGain * (1.0 + out.a1 + out.a2)
                                 / (out.b0 + out.b1 + out.b2);

        out.b0 *= scale;
        out.b1 *= scale;
        out.b2 *= scale;

        return out;
    };

    if (zeroHz < 0.45 * sampleRate)
        return makeCoeffs (prewarp (zeroHz));

    // The analog zeros lie above Nyquist, so no biquad at this rate can place
    // them faithfully. Two tempting shortcuts are both wrong:
    //
    //   - Mapping 90 kHz through the plain bilinear folds the zeros down to about
    //     20 kHz and flattens the response early: 4 dB of error at 15 kHz at a
    //     44.1 kHz sample rate.
    //   - Dropping the zeros gives a plain two-pole low-pass, whose bilinear form
    //     always has a double zero at z = -1. Inverting that yields a double pole
    //     on the unit circle, so the decoder becomes marginally unstable.
    //
    // The second case is exactly what section 2.6 means when it says the shelves
    // "provide phase characteristics that are essential in the decoding mode":
    // without a shelf the network is not invertible at all. So the zeros are kept,
    // and their digital frequency is solved for to best match the analog magnitude
    // across the band while staying clear of Nyquist.
    const double lowHz = 0.20 * sampleRate;
    const double highHz = 0.48 * sampleRate;

    // Log-spaced probe frequencies up to the top of the useful band. The analog
    // reference at each probe does not depend on the candidate, so it is
    // evaluated once here rather than once per candidate.
    std::vector<double> probes, analogDb;
    for (double hz = 100.0; hz < 0.40 * sampleRate; hz *= 1.15)
    {
        probes.push_back (hz);
        analogDb.push_back (20.0 * std::log10 (std::abs (responseAt (hz))));
    }

    double bestHz = highHz;
    double bestError = -1.0;

    constexpr int steps = 512;

    for (int i = 0; i <= steps; ++i)
    {
        const double t = (double) i / (double) steps;
        const double candidateHz = lowHz * std::pow (highHz / lowHz, t);
        const auto candidate = makeCoeffs (prewarp (candidateHz));

        double sumSquared = 0.0;

        for (size_t p = 0; p < probes.size(); ++p)
        {
            const double digitalDb = 20.0 * std::log10 (std::abs (candidate.responseAt (twoPi * probes[p] / sampleRate)));
            const double err = digitalDb - analogDb[p];
            sumSquared += err * err;
        }

        if (bestError < 0.0 || sumSquared < bestError)
        {
            bestError = sumSquared;
            bestHz = candidateHz;
        }
    }

    return makeCoeffs (prewarp (bestHz));
}

//==============================================================================
namespace srNetworks
{

std::vector<ShelfPrototype> skewing()
{
    return {
        ShelfPrototype::lowPassShelf  (skewHfCornerHz, skewHfShelfDb),
        ShelfPrototype::highPassShelf (skewLfCornerHz, skewLfShelfDb)
    };
}

std::vector<ShelfPrototype> antiSaturation()
{
    return {
        ShelfPrototype::lowPassShelf  (antiSatHfCornerHz, antiSatHfShelfDb),
        ShelfPrototype::highPassShelf (antiSatLfCornerHz, antiSatLfShelfDb)
    };
}

std::vector<ShelfPrototype> mainPath()
{
    auto protos = skewing();

    for (const auto& p : antiSaturation())
        protos.push_back (p);

    return protos;
}

std::complex<double> responseAt (const std::vector<ShelfPrototype>& protos, double hz)
{
    std::complex<double> h { 1.0, 0.0 };

    for (const auto& p : protos)
        h *= p.responseAt (hz);

    return h;
}

std::vector<BiquadCoeffs> toBiquads (const std::vector<ShelfPrototype>& protos,
                                    double sampleRate)
{
    std::vector<BiquadCoeffs> out;
    out.reserve (protos.size());

    for (const auto& p : protos)
        out.push_back (p.toBiquad (sampleRate));

    return out;
}

const std::vector<AntisaturationTarget>& publishedTargets()
{
    // Section 4.1: "an antisaturation effect of about 2 dB at 5 kHz, 6 dB at
    // 10 kHz, and 10 dB at 25 Hz and 15 kHz".
    static const std::vector<AntisaturationTarget> targets
    {
        {    25.0, 10.0 },
        {  5000.0,  2.0 },
        { 10000.0,  6.0 },
        { 15000.0, 10.0 }
    };

    return targets;
}

namespace
{
    // Frequencies held flat so the fixed networks do not encroach on the midband,
    // which the least treatment principle and the flat midband of Fig. 13 require.
    const std::vector<double> midbandProbes { 200.0, 400.0, 700.0, 1500.0, 2500.0 };
    constexpr double midbandWeight = 4.0;

    double cascadeDb (const std::vector<BiquadCoeffs>& coeffs, double hz, double sampleRate)
    {
        std::complex<double> h { 1.0, 0.0 };
        const double omega = twoPi * hz / sampleRate;

        for (const auto& c : coeffs)
            h *= c.responseAt (omega);

        return 20.0 * std::log10 (std::abs (h));
    }

    // Each prototype maps to its biquad independently, so a cascade can be
    // assembled from ready-made skewing sections plus fresh antisaturation
    // ones. The skewing networks are fixed constants, and the 12 kHz low-pass
    // shelf has its zeros above Nyquist at every supported rate, so converting
    // it triggers the zero-placement search — worth doing once per fit rather
    // than once per objective evaluation.
    std::vector<BiquadCoeffs> assemble (double sampleRate,
                                        const std::vector<BiquadCoeffs>& skewCoeffs,
                                        double hfHz, double hfDb,
                                        double lfHz, double lfDb)
    {
        std::vector<BiquadCoeffs> out = skewCoeffs;
        out.push_back (ShelfPrototype::lowPassShelf  (hfHz, hfDb).toBiquad (sampleRate));
        out.push_back (ShelfPrototype::highPassShelf (lfHz, lfDb).toBiquad (sampleRate));

        return out;
    }

    double objective (double sampleRate, const std::vector<BiquadCoeffs>& skewCoeffs,
                      double hfHz, double hfDb, double lfHz, double lfDb)
    {
        const auto coeffs = assemble (sampleRate, skewCoeffs, hfHz, hfDb, lfHz, lfDb);
        const double ref = cascadeDb (coeffs, midbandReferenceHz, sampleRate);

        double sum = 0.0;

        for (const auto& t : publishedTargets())
        {
            if (t.hz > 0.45 * sampleRate)
                continue;   // target above the band at this rate; nothing to fit

            const double err = (cascadeDb (coeffs, t.hz, sampleRate) - ref) + t.attenuationDb;
            sum += err * err;
        }

        for (double hz : midbandProbes)
        {
            const double err = midbandWeight * (cascadeDb (coeffs, hz, sampleRate) - ref);
            sum += err * err;
        }

        return sum;
    }
}

std::vector<BiquadCoeffs> mainPathBiquadsForRate (double sampleRate)
{
    // The fit is a pure function of the sample rate, and every caller runs on
    // the message thread or inside prepare() — never per block — so repeated
    // requests for the same rate (editor redraws, replicated prepares) are
    // served from a memo instead of re-running the descent. Keyed on 1/16th of
    // a hertz, far finer than any two rates a host distinguishes.
    static std::mutex cacheMutex;
    static std::map<long long, std::vector<BiquadCoeffs>> cache;

    const long long key = std::llround (sampleRate * 16.0);

    {
        const std::lock_guard<std::mutex> lock (cacheMutex);
        const auto it = cache.find (key);

        if (it != cache.end())
            return it->second;
    }

    // Coordinate descent with shrinking steps, started from the offline analog
    // fit. The objective is smooth and the start is close, so this converges in a
    // handful of passes; a general optimiser would be overkill.
    double hfHz = antiSatHfCornerHz, hfDb = antiSatHfShelfDb;
    double lfHz = antiSatLfCornerHz, lfDb = antiSatLfShelfDb;

    double stepHfHz = 800.0, stepHfDb = 1.5;
    double stepLfHz = 20.0, stepLfDb = 0.5;

    const auto skewCoeffs = toBiquads (skewing(), sampleRate);

    double best = objective (sampleRate, skewCoeffs, hfHz, hfDb, lfHz, lfDb);

    for (int pass = 0; pass < 60; ++pass)
    {
        bool improved = false;

        struct Axis { double* value, * step, low, high; };

        Axis axes[] = {
            { &hfHz, &stepHfHz, 2000.0, 0.40 * sampleRate },
            { &hfDb, &stepHfDb,    0.1,  20.0 },
            { &lfHz, &stepLfHz,   40.0, 400.0 },
            { &lfDb, &stepLfDb,    0.1,  20.0 }
        };

        for (auto& axis : axes)
        {
            for (double direction : { +1.0, -1.0 })
            {
                const double original = *axis.value;
                const double trial = std::min (axis.high, std::max (axis.low, original + direction * *axis.step));

                if (trial == original)
                    continue;

                *axis.value = trial;
                const double score = objective (sampleRate, skewCoeffs, hfHz, hfDb, lfHz, lfDb);

                if (score < best)
                {
                    best = score;
                    improved = true;
                    break;      // keep the improvement and move to the next axis
                }

                *axis.value = original;
            }
        }

        if (! improved)
        {
            stepHfHz *= 0.5;
            stepHfDb *= 0.5;
            stepLfHz *= 0.5;
            stepLfDb *= 0.5;

            if (stepHfHz < 0.01 && stepHfDb < 1.0e-5 && stepLfHz < 1.0e-4 && stepLfDb < 1.0e-6)
                break;
        }
    }

    auto result = assemble (sampleRate, skewCoeffs, hfHz, hfDb, lfHz, lfDb);

    const std::lock_guard<std::mutex> lock (cacheMutex);
    return cache.emplace (key, std::move (result)).first->second;
}

double worstTargetErrorDb (const std::vector<BiquadCoeffs>& coeffs, double sampleRate)
{
    const double ref = cascadeDb (coeffs, midbandReferenceHz, sampleRate);
    double worst = 0.0;

    for (const auto& t : publishedTargets())
    {
        if (t.hz > 0.45 * sampleRate)
            continue;

        worst = std::max (worst, std::abs ((cascadeDb (coeffs, t.hz, sampleRate) - ref) + t.attenuationDb));
    }

    return worst;
}

double worstMidbandErrorDb (const std::vector<BiquadCoeffs>& coeffs, double sampleRate)
{
    const double ref = cascadeDb (coeffs, midbandReferenceHz, sampleRate);
    double worst = 0.0;

    for (double hz : midbandProbes)
        worst = std::max (worst, std::abs (cascadeDb (coeffs, hz, sampleRate) - ref));

    return worst;
}

} // namespace srNetworks

} // namespace lime
