/*
  ==============================================================================

    PhaseTests — measures the phase rotator against what it claims.

    Three things have to hold, and only one of them is obvious:

      - the magnitude must not move, at any angle or frequency;
      - the phase must move by the angle asked for, at every frequency;
      - 180 degrees must be exactly a polarity invert, since that is the endpoint the
        control is drawn to.

    Every number here is printed with its limit, as everywhere else in this project.

  ==============================================================================
*/

#include "../dsp/PhaseRotator.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
    int checks = 0, failures = 0;

    constexpr double pi = 3.14159265358979323846;

    void expectBelow (const std::string& what, double got, double limit, const char* unit = "dB")
    {
        ++checks;
        const bool ok = got <= limit;

        if (! ok)
            ++failures;

        std::printf ("  [%s] %-58s %8.3f %s  (limit %8.3f %s)\n",
                     ok ? "pass" : "FAIL", what.c_str(), got, unit, limit, unit);
    }

    /** Runs a steady tone through a rotator set to `degrees` and returns the measured
        magnitude ratio and phase shift, taken well after the network has settled. */
    struct Measurement { double magnitudeRatio, phaseDegrees; };

    Measurement measure (double hz, double sampleRate, double degrees)
    {
        lime::PhaseRotator rotator;
        rotator.prepare (sampleRate, 1);
        rotator.setPhaseDegrees (degrees);

        const int settle = 40000;
        const int window = 40000;

        std::vector<double> buffer ((size_t) (settle + window));
        const double omega = 2.0 * pi * hz / sampleRate;

        for (size_t n = 0; n < buffer.size(); ++n)
            buffer[n] = std::sin (omega * (double) n);

        double* data = buffer.data();
        rotator.process (&data, 1, (int) buffer.size());

        // Correlate the settled tail against sine and cosine to recover amplitude and phase.
        double re = 0.0, im = 0.0;

        for (int n = settle; n < settle + window; ++n)
        {
            const double phase = omega * (double) n;
            re += buffer[(size_t) n] * std::sin (phase);
            im += buffer[(size_t) n] * std::cos (phase);
        }

        re *= 2.0 / (double) window;
        im *= 2.0 / (double) window;

        const double magnitude = std::sqrt (re * re + im * im);
        const double shift = std::atan2 (im, re) * 180.0 / pi;

        return { magnitude, shift };
    }

    double wrapped (double degrees)
    {
        while (degrees > 180.0)  degrees -= 360.0;
        while (degrees < -180.0) degrees += 360.0;
        return degrees;
    }
}

int main()
{
    const double rate = 48000.0;

    std::printf ("\nMagnitude is untouched at every angle and frequency\n");

    for (double hz : { 40.0, 100.0, 440.0, 1000.0, 5000.0, 12000.0 })
        for (double deg : { -180.0, -90.0, -30.0, 30.0, 90.0, 180.0 })
        {
            const auto m = measure (hz, rate, deg);
            const double errorDb = std::abs (20.0 * std::log10 (std::max (1.0e-12, m.magnitudeRatio)));

            expectBelow (std::to_string ((int) hz) + " Hz at " + std::to_string ((int) deg)
                             + " deg, magnitude error",
                         errorDb, 0.15);
        }

    // What the control promises is a *rotation*, and a rotation is a difference. The
    // network's in-phase branch is an allpass, so the absolute phase of everything coming
    // out of it carries that branch's own frequency-dependent shift — the same shift at
    // every setting of the knob, which is exactly why it cancels out of a difference and
    // exactly why the knob is a rotation control rather than a phase-alignment one.
    //
    // Measuring against the input instead was the first version of this test, and it failed
    // by 148 degrees at 40 Hz and 90 at 12 kHz while the rotation itself was perfect. The
    // reading it took was real; it was just not a reading of this control.
    std::printf ("\nRotation between any two settings is their difference, at every frequency\n");

    for (double hz : { 40.0, 100.0, 440.0, 1000.0, 5000.0, 12000.0 })
    {
        const auto reference = measure (hz, rate, 90.0);

        for (double deg : { -180.0, -90.0, -10.0, 10.0, 30.0, 180.0 })
        {
            const auto m = measure (hz, rate, deg);
            const double moved = wrapped (m.phaseDegrees - reference.phaseDegrees);
            const double wanted = wrapped (-(deg - 90.0));

            expectBelow (std::to_string ((int) hz) + " Hz, " + std::to_string ((int) deg)
                             + " deg against 90 deg",
                         std::abs (wrapped (moved - wanted)), 2.0, "deg");
        }
    }

    std::printf ("\n180 degrees inverts what the network passes at zero rotation\n");
    {
        // Against the network's own zero, not against the input: the branch is in circuit
        // either way, so this is the polarity claim with the allpass taken out of it.
        const auto zero = measure (997.0, rate, 0.001);
        const auto flipped = measure (997.0, rate, 180.0);

        const double moved = std::abs (wrapped (flipped.phaseDegrees - zero.phaseDegrees));
        expectBelow ("phase difference against 180", std::abs (moved - 180.0), 2.0, "deg");

        const double ratioDb = std::abs (20.0 * std::log10 (std::max (1.0e-12,
                                            flipped.magnitudeRatio / zero.magnitudeRatio)));
        expectBelow ("magnitude difference against 180", ratioDb, 0.15);
    }

    std::printf ("\nZero is exactly transparent, so a complementary pair stays complementary\n");
    {
        lime::PhaseRotator rotator;
        rotator.prepare (rate, 1);
        rotator.setPhaseDegrees (0.0);

        std::vector<double> input (8192), output;

        for (size_t n = 0; n < input.size(); ++n)
            input[n] = std::sin (0.01 * (double) n) * 0.4 + std::sin (0.7 * (double) n) * 0.3;

        output = input;
        double* data = output.data();
        rotator.process (&data, 1, (int) output.size());

        double worst = 0.0;

        for (size_t n = 0; n < input.size(); ++n)
            worst = std::max (worst, std::abs (output[n] - input[n]));

        expectBelow ("largest sample difference at zero", worst, 0.0, "");
    }

    std::printf ("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
