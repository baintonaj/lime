/*
  ==============================================================================

    WOLA window pair.

    Analysis and synthesis both use sqrt-Hann, so their product is a Hann window,
    and a periodic Hann overlap-adds to a constant at any hop that divides the
    size into two or more frames — the single cosine it is made of sums to zero
    over every such grid. Rather than relying on the closed-form normalising
    constant, the constant is measured and divided out, and the deviation from
    flatness is exposed so a test can assert perfect reconstruction directly;
    the suite measures the power-of-two hops from size/2 down to 256.

  ==============================================================================
*/

#pragma once

#include <cmath>
#include <vector>

namespace lime
{

class WolaWindow
{
public:
    WolaWindow (int sizeIn, int hopIn)
        : size (sizeIn), hop (hopIn)
    {
        constexpr double twoPi = 6.283185307179586476925286766559;

        analysisWindow.resize ((size_t) size);
        synthesisWindow.resize ((size_t) size);

        // Periodic Hann, square-rooted so that analysis * synthesis == Hann.
        for (int n = 0; n < size; ++n)
        {
            const double hann = 0.5 - 0.5 * std::cos (twoPi * (double) n / (double) size);
            const double root = std::sqrt (hann);

            analysisWindow[(size_t) n] = root;
            synthesisWindow[(size_t) n] = root;
        }

        // Measure the overlap-add sum of the window product over one hop period.
        // Every phase within a hop must give the same total for reconstruction
        // to be exact, so the largest and smallest are both tracked.
        double smallest = 0.0, largest = 0.0;

        for (int phase = 0; phase < hop; ++phase)
        {
            double sum = 0.0;

            for (int n = phase; n < size; n += hop)
                sum += analysisWindow[(size_t) n] * synthesisWindow[(size_t) n];

            if (phase == 0)
                smallest = largest = sum;
            else
            {
                smallest = std::min (smallest, sum);
                largest = std::max (largest, sum);
            }
        }

        colaSum = 0.5 * (smallest + largest);
        deviation = (colaSum > 0.0) ? (largest - smallest) / colaSum : 1.0;

        // Fold the normalisation into the synthesis window.
        if (colaSum > 0.0)
            for (auto& w : synthesisWindow)
                w /= colaSum;
    }

    int getSize() const noexcept { return size; }
    int getHop()  const noexcept { return hop; }

    const double* analysis()  const noexcept { return analysisWindow.data(); }
    const double* synthesis() const noexcept { return synthesisWindow.data(); }

    /** Relative spread of the overlap-add sum across one hop period. Zero means
        the window pair reconstructs exactly at this hop. */
    double colaDeviation() const noexcept { return deviation; }

    /** Sum of the analysis window.

        Needed to convert bin magnitudes into physical amplitude. A tone of amplitude
        A at a bin centre produces a bin magnitude of A * sum / 2, so a control
        threshold quoted in dB relative to a full-scale sine only means what it says
        once that factor is divided out. */
    double analysisSum() const noexcept
    {
        double sum = 0.0;

        for (auto w : analysisWindow)
            sum += w;

        return sum;
    }

private:
    int size = 0, hop = 0;
    double colaSum = 1.0, deviation = 0.0;
    std::vector<double> analysisWindow, synthesisWindow;
};

} // namespace lime
