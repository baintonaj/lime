/*
  ==============================================================================

    FilterDesign — the fixed networks of the process B main path.

    Spectral skewing (section 2.6) is taken verbatim from the paper. The
    antisaturation networks (section 2.7) are described there only as "operative
    above about 5 kHz and below about 100 Hz", so their corners and depths were
    solved offline such that the combined response reproduces the antisaturation
    figures quoted in section 4.1. See the constants below for the fit result.

  ==============================================================================
*/

#pragma once

#include "Biquad64.h"

#include <complex>
#include <vector>

namespace lime
{

/** A second-order analog shelving prototype,

        H(s) = K * (s^2 + sqrt(2)*wz*s + wz^2) / (s^2 + sqrt(2)*wp*s + wp^2)

    with Butterworth damping on both quadratics. This single form expresses every
    network SR needs: a low-pass rolling off from `poleHz` and flattening to a
    shelf, or a high-pass doing the same, depending on whether the zero sits above
    or below the pole.
*/
struct ShelfPrototype
{
    double poleHz = 1000.0;
    double zeroHz = 1000.0;
    double gain = 1.0;        // K

    /** Low-pass shelf: unity at DC, `shelfDb` of attenuation at high frequency,
        with the -3 dB Butterworth corner at `cornerHz`. */
    static ShelfPrototype lowPassShelf (double cornerHz, double shelfDb);

    /** High-pass shelf: unity at high frequency, `shelfDb` of attenuation at DC,
        with the -3 dB Butterworth corner at `cornerHz`. */
    static ShelfPrototype highPassShelf (double cornerHz, double shelfDb);

    /** Analog response at a frequency in Hz. */
    std::complex<double> responseAt (double hz) const;

    /** Bilinear transform to a digital biquad.

        Critical frequencies below 0.45 * sampleRate are prewarped so they land
        on their intended value in Hz; above that the prewarp tangent blows up, so
        the frequency is mapped by the plain bilinear instead. That only ever
        affects the skewing zeros, which sit far above the audio band (the 12 kHz
        skewing shelf is not reached below Nyquist at all), and any residual
        in-band phase difference cancels exactly between encoder and decoder
        because they use inverse coefficient sets.
    */
    BiquadCoeffs toBiquad (double sampleRate) const;
};

//==============================================================================
/** The four fixed networks, as published and as fitted. */
namespace srNetworks
{
    // Section 2.6, stated verbatim.
    inline constexpr double skewHfCornerHz = 12000.0;
    inline constexpr double skewHfShelfDb = 35.0;
    inline constexpr double skewLfCornerHz = 40.0;
    inline constexpr double skewLfShelfDb = 25.0;

    // Section 2.7 gives only "above about 5 kHz" and "below about 100 Hz". These
    // were solved for so the combined response matches the section 4.1 figures
    // (2 dB at 5 kHz, 6 dB at 10 kHz, 10 dB at 25 Hz and 15 kHz) with the midband
    // held flat. Two-pole shelves fit to 0.27 dB worst case; one-pole reached only
    // 0.77 dB. The solved corners landed on the paper's own stated values without
    // being constrained to them.
    inline constexpr double antiSatHfCornerHz = 5204.628974;
    inline constexpr double antiSatHfShelfDb = 4.819522;
    inline constexpr double antiSatLfCornerHz = 95.843195;
    inline constexpr double antiSatLfShelfDb = 1.283246;

    /** Spectral skewing alone: what drives the stage circuits. */
    std::vector<ShelfPrototype> skewing();

    /** Antisaturation alone. */
    std::vector<ShelfPrototype> antiSaturation();

    /** The complete fixed main path, skewing followed by antisaturation. */
    std::vector<ShelfPrototype> mainPath();

    /** Combined analog response of a prototype list at a frequency in Hz. */
    std::complex<double> responseAt (const std::vector<ShelfPrototype>& protos, double hz);

    /** Converts a prototype list to biquad coefficients at a sample rate. */
    std::vector<BiquadCoeffs> toBiquads (const std::vector<ShelfPrototype>& protos,
                                        double sampleRate);

    /** The section 4.1 antisaturation figures: attenuation in dB relative to the
        midband. These are the paper's own numbers and are the acceptance criteria
        the main path is designed against. */
    struct AntisaturationTarget { double hz, attenuationDb; };

    const std::vector<AntisaturationTarget>& publishedTargets();

    inline constexpr double midbandReferenceHz = 1000.0;

    /** Builds the main-path biquads for a sample rate, solving for the
        antisaturation parameters so that the *digital* cascade reproduces
        publishedTargets().

        Fitting in the digital domain rather than the analog one matters. The
        skewing low-pass has its shelf zeros above Nyquist, so its digital
        realisation cannot track the analog prototype exactly near the top of the
        band — roughly 1 dB at 15 kHz at 44.1 kHz. Since antisaturation is the
        fitted half of the design anyway, letting it absorb that deviation makes
        the total, which is what section 4.1 actually specifies, come out right.

        Deterministic, and memoised per rate behind a lock, so repeated calls
        from the editor or from replicated prepares cost a lookup. Call it from
        prepare() or the message thread; it is not for the audio thread.
    */
    std::vector<BiquadCoeffs> mainPathBiquadsForRate (double sampleRate);

    /** Worst deviation of a cascade from publishedTargets(), in dB. */
    double worstTargetErrorDb (const std::vector<BiquadCoeffs>& coeffs, double sampleRate);

    /** Worst midband deviation from flat, in dB. */
    double worstMidbandErrorDb (const std::vector<BiquadCoeffs>& coeffs, double sampleRate);
}

} // namespace lime
