# Lime — design record

A companding noise reduction system reinterpreted as a per-bin STFT process. This document
is the record of the design as built: the specification extracted from the source papers,
every decision with its reasoning, and the measurements — negative results included — that
forced each reversal.

## Context

Lime is a JUCE plugin (AU + VST3, C++20, `pluginVST3Category="Dynamics"`), begun from an
empty Projucer template. It takes the architecture of the **spectral recording** process —
the analog multiband tape noise-reduction system introduced in 1986 and described in the
1987 paper — and extends it into an **STFT
process where every bin carries its own unit**. The entire project processes in **double
precision**.

The DSP constants were taken from the published papers — **"The Spectral Recording
Process", JAES 35(3) 1987** and **"An Audio Noise Reduction System", JAES 15 1967** — which
are the source of truth. The reference hardware's service documentation was consulted
alongside them. Neither document is redistributed here; every constant taken from them is
cited in the source by section.

---

## Reference specification (from the papers)

### Dual-path differential topology — the foundation (1967 paper)

Main path passes the signal unmodified; a side chain produces a small differential
component, **added on encode, subtracted on decode**:

```
encode:  y = (1 + G(x))·x
decode:  z = y − z·G(z)   ⇒   z = y / (1 + G(z))
```

`z = x` exactly when the two operators are identical. The decoder is a **feedback**
configuration — that is what makes it complementary. Consequences: the side chain is
negligible at high level, so low distortion, low overshoot, and tolerance to channel
gain error (worst-case tracking error ≈ the gain error in dB, ~30 dB below peak).

### SR stage layout

| Property | Value |
|---|---|
| Stages | 3 HF + 2 LF, series (staggered), crossover **800 Hz** |
| Thresholds | **−30, −48, −62 dB** rel. reference |
| Gain per stage | ~**8 dB** → **24 dB total HF, 16 dB LF**; ~1 dB more above reference |
| Band-defining filters | **single-pole** (HF stages reach to ~200 Hz, LF to ~3 kHz) |
| Dynamic action range | LF −48…−5 dB, HF −62…−5 dB, ~**2:1** ratio between |
| reference level at module boundary | 388 mV (the reference module) |

Series staggering multiplies transfer functions (dB add) and **steepens spectral
discrimination**. The HF and LF band-defining filters are single-pole complements at the
same corner, so `jω/(jω+ω₀) + ω₀/(jω+ω₀) ≡ 1` — they sum to **exactly unity**.

The 800 Hz corner in that table is the reference design's and is what Lime uses. It was briefly
a user control and is a compiled-in constant again; how finely the complementary pair is sampled
across the band was then briefly a control as well, and since v1.0.0 is a compiled-in constant
too — 2048 sections. See *The crossover control became a resolution control* and *The sections
control retired* under v1.0.0.

### Action substitution (§2.2)

Fixed band and sliding band fed in parallel; output taken from the sliding band, whose
variable element is referenced to the fixed band's output:

```
Vo/Vi = F1 + F2 − F1·F2   ≡   1 − (1−F1)(1−F2)
```

Gives fixed-band behaviour on the **stop-band** side of a dominant (uniform loss, so
channel response errors are not exaggerated) and sliding-band behaviour on the
**pass-band** side (NR maintained above a dominant HF component / below a dominant LF
one). The fixed band also supplies fast frequency-independent recovery, letting the
sliding band use long time constants — low modulation distortion *and* fast recovery.

### Modulation control (§2.3)

Frequency-weighted rectified main-path signals fed **in opposition** to each stage's
control signal, so gain stops moving once the required attenuation is reached. Justified
by the phasor argument (Fig. 5): in the stop band the side chain is near quadrature to
the main path, so it barely affects total amplitude — threshold can rise, limiting can
weaken. Each MC is a **single scalar control voltage** in the analog.

| Signal | Weighting | Smoothing | Controls |
|---|---|---|---|
| MC1 | 3 kHz 1-pole HP + all-pass | — | HF sliding, steady state |
| MC2 | = MC1 | 2-stage 1 ms | HF sliding O/S |
| MC3 | 400 Hz + 800 Hz 1-pole LP | — | HF fixed, steady + transient |
| MC4 | 200 Hz 1-pole LP + all-pass | — | LF sliding |
| MC5 | = MC4 | 2-stage 2 ms | LF sliding O/S |
| MC6 | 800 Hz + 1.6 kHz 1-pole HP | — | LF fixed |
| MC7 | = MC6 | 2-stage 2 ms | LF fixed O/S, + LF sliding O/S |
| MC8 | 5 kHz HP, rect, double-diff 15 µs, peak-hold 30 ms | — | inhibits LF O/S under HF transients |

MC1–MC7 taken from the 2nd-stage adder output; MC8 from just after spectral skewing.

### Per-stage control detail

**HF stage** (Fig. 9): input via 800 Hz 1-pole HP then 400 Hz 1-pole HP; output through a
reciprocal 400 Hz network, so quiescent stage response is a single-pole 800 Hz HP.
- *Fixed band*: two control paths into a **maximum selector** — main (rectify, oppose
  MC3, 15 ms) and **pass-band** (1.6 kHz 1-pole HP, rectify, 15 ms); selector output
  further smoothed **160 ms**. Handles simple- vs complex-signal conditions.
- *Sliding band*: control from stage output via ~10 kHz 1-pole HP (differs per stage),
  rectify, oppose MC1 (the ratio forms a **sliding end-stop**), smooth ~5 ms then
  **80 ms**. Fixed-band signal combined in *opposition* with the sliding output to raise
  the sliding threshold at HF.
- Variable filter is **quasi-two-pole**: 1-pole fixed filter + variable shelf, variable
  turnover one octave into the stop band relative to the fixed cutoff.

**LF stage** (Figs. 10–11): sliding band acts *downward*; the fixed 800 Hz band-defining
LP **follows** the variable filters.
- *Fixed band*: weighting = cascaded 800 Hz + 1.6 kHz 1-pole LP; main path rectify,
  oppose MC6, 15 ms; pass-band path adds 400 Hz 1-pole LP, rectify, 15 ms; max selector
  → **300 ms**.
- *Sliding band*: control from stage output via 80 Hz 1-pole LP, rectify, oppose MC4,
  7.5 ms then **150 ms**.

### Overshoot suppression (§2.4)

Thresholds ~**10 dB above** the steady-state thresholds, so they rarely act. Unsmoothed
rectified control signals, opposed by the relevant MC, diode-coupled into the final
integrators. Primary + secondary everywhere, plus a gentle slow LF O/S below ~200 Hz.
Thresholds are *derived from the MC signals* so transient and steady-state track.

### Spectral skewing (§2.6) — main path, ahead of the stages

- HF: **12 kHz 2-pole** Butterworth-like LP with a **~35 dB shelf**.
- LF: **40 Hz 2-pole** Butterworth-like HP with a **~25 dB shelf**.

The shelves exist to provide **phase characteristics essential to decoding** — these must
be modelled as complex transfer functions, not magnitude curves.

### Antisaturation (§2.7) — fixed main-path shelves

Above ~5 kHz and below ~100 Hz. Combined with the secondary antisaturation effect of the
skewing networks plus ~1 dB wideband compensation in the stage signal combiner, the total
is ≈ **2 dB @ 5 kHz, 6 dB @ 10 kHz, 10 dB @ 25 Hz and 15 kHz**.

### Low-level target curve (Fig. 13)

Encode boost resembles the inverse **CCIR noise-weighting** curve; decode resembles the
low-level **Fletcher–Munson / Robinson–Dadson** contours. ~24 dB in the mid/upper-mid,
tapering at both extremes. SR deliberately *slightly degrades* sub-20 Hz noise.

### A-type (for the switchable mode)

Four bands — **≤80 Hz LP, 80 Hz–3 kHz BP, ≥3 kHz HP, ≥9 kHz HP** (bands 1/3/4 are
12 dB/oct; band 2 is complementary to 1 and 3). Thresholds **40 dB below peak**. Sums to
**10 dB flat, rising to 15 dB at 15 kHz**.

### Calibration (§5)

**calibration noise**: pink noise with **20 ms nicks every 2 s**, recorded **15 dB below** reference level
level. **Auto Compare**: decoder alternates tape and internal reference noise in **4 s**
segments (8 s cycle). **calibration tone** (A-type): 850 Hz, FM'd up 10% for 30 ms every 750 ms.

---

## Agreed design decisions

Settled interactively with the user:

| # | Decision |
|---|---|
| 1 | **Hybrid architecture**: per-bin staggered fixed-band stack + a *real* sliding layer |
| 2 | Sliding corners: **all local maxima above a masking spread threshold** (count unbounded) |
| 3 | Skirt shape: **faithful quasi-two-pole** (1-pole fixed + variable shelf, turnover 1 octave into stop band) |
| 4 | Skirts combine by **generalised action substitution** `F = 1 − Π(1−Fᵢ)` — **revised, see below** |
| 5 | Decoder: **true feedback from smoothed state**, one frame delayed |
| 6 | Overshoot suppression replaced by a **per-frame gain-slew clamp tracking the MC signals** |
| 7 | Main path: **STFT-domain complex bin multipliers**, dedicated **8192** transform |
| 8 | Side chains: **dual resolution split at 800 Hz** — LF 4096, HF 1024, common hop 256. The corner was made a control and then fixed again; its *resolution* followed the same arc — a control at v0.2.0, pinned at 2048 sections at v1.0.0 — see below |
| 9 | Window sizes **scale with sample rate**, snapped to powers of two, **capped** |
| 10 | FFT: **Accelerate `vDSP_*D`** with a **project-local portable double FFT** fallback |
| 11 | Modes: **Encode**, **Decode**, and **Encode → channel → Decode** loop |
| 12 | Channel (loop mode): **full tape model** — noise, saturation, HF loss, modulation noise |
| 13 | **SR and A-type switchable**; A-type is a **faithful 4-band** build on the shared core |
| 14 | Parameters **faithful and fixed** — only the reference hardware's real controls; internal constants compiled in |
| 15 | **calibration-noise generator + Auto Compare + spectral display** in scope |

### Sizes and latency

**Revised after building and measuring the transform layer — see the note below.**

At 48 kHz: main path FFT 8192 with a 4096-tap filter and a **4096-sample block advance**;
LF side chain 4096 window; HF side chain 1024 window; common hop 256.

| Path | Engine | Latency | Alignment delay |
|---|---|---|---|
| Main (skew + antisaturation) | `OverlapSave64` 8192, 4096 taps, advance 4096 | 4096 | 0 |
| LF side chain | `Wola64` 4096, hop 256 | 4096 | 0 |
| HF side chain | `Wola64` 1024, hop 256 | 1024 | 3072 |

**Total plugin latency: 4096 samples (85.3 ms at 48 kHz)** — half the 8192 the plan
originally assumed. The original figure came from imagining the main path as a 32×-overlap
WOLA at 8192, which costs 8192 samples of delay. Overlap-save (already identified below as
the correct method for a time-invariant filter) costs only its block advance, and that
advance can be set to 4096 — the largest value that still discards enough of each block to
avoid wrap-around, since `advance <= N - taps + 1 = 4097`. Setting it to exactly 4096 makes
the main path land on the same latency as the LF side chain, so neither needs a delay line
and only the HF path is padded.

Scaling is a single power-of-two exponent relative to 48 kHz,
`exponent = clamp(round(log2(sr/48000)), -2, +1)`, applied to every size at once:

| Exponent | Rates | Main | LF | HF | Hop | Latency |
|---|---|---|---|---|---|---|
| −2 | ≤ ~17 kHz | 2048 | 1024 | 256 | 64 | 1024 |
| −1 | ~17–34 kHz | 4096 | 2048 | 512 | 128 | 2048 |
| 0 | ~34–68 kHz | 8192 | 4096 | 1024 | 256 | 4096 |
| +1 | ≥ ~68 kHz | 16384 | 8192 | 2048 | 512 | 8192 |

Scaling in **both** directions is what keeps bin width constant in Hz and latency constant
in milliseconds, which is the whole point given the paper's constants are specified in Hz and
ms. The original plan clamped the exponent at 0, which `auval` exposed as wrong: it renders
at 11025 Hz and 22050 Hz, and at 11 kHz a fixed 4096-sample window means 372 ms of latency
and 2.7 Hz bins — far finer than anything the process needs. The upper clamp is kept, and
is purely a cost decision: 176.4/192 kHz reuse the 96 kHz sizes, coarsening bins to 11.7 Hz
rather than 5.9 Hz. The 40 Hz skewing corner still lands near bin 3, which is adequate, and
latency in milliseconds improves at those rates.

`Wola64`'s latency is exactly its window size, not `size - hop`. A frame's oldest hop is
complete as soon as that frame is added, which alone would cost `size - hop`, but the
engine consumes a hop before producing one and that adds the remaining hop back. Reclaiming
it would need priming logic that breaks on partial blocks, and the alignment delays absorb
it for free. Measured and asserted by test rather than assumed.

---

## Main path: reversed to time-domain biquads, forced by measurement

**The main path is implemented as cascaded double-precision biquads, not as an STFT
multiplier.** This reverses the earlier decision, and the reversal is forced rather than
preferential.

The concern raised when the decision was made was impulse-response length: applying a fixed
filter by per-bin multiplication is a circular convolution, so the response must fit inside
the window. The 40 Hz two-pole skewing high-pass fits inside 8192 (85 ms) comfortably. What
was not checked at the time is the **inverse** network, which the decoder needs.

Measured, on a 2^19 design grid at 48 kHz:

| Network | Energy in first 22 ms | Energy at negative time | Peak index |
|---|---|---|---|
| Forward (skew × antisat) | **99.30 %** | 0.70 % | 1 |
| Inverse (de-skew × de-antisat) | **20.55 %** | **78.87 %** | L−1, i.e. t = −1 sample |

The inverse response is overwhelmingly **anti-causal** when obtained by frequency sampling.
A causal 4096-tap FIR would retain about 21 % of its energy, so the encode→decode null
would be limited to roughly −7 dB. No increase in tap count helps: the residual plateaus,
because the missing energy is at negative time, not in a long tail. Enforcing causality by
minimum-phase reconstruction would discard the phase — which §2.6 states explicitly is
"essential in the decoding mode".

Biquads have none of these problems. Each network is a second-order analog shelf, so the
bilinear transform maps it to a biquad with no approximation in structure; the networks are
minimum-phase, so the inverse is obtained exactly by swapping numerator and denominator
coefficients and is stable; there is no truncation, no wrap-around, and no added latency.
The encode→decode null through the main path becomes exact to numerical precision.

Consequences: `OverlapSave64` is no longer used by the main path (it stays in the codebase,
tested, in case a genuinely fixed filter is wanted later), the dedicated 8192-point
transform disappears along with two transforms per hop, and the main path now needs a
4096-sample alignment delay to match the LF side chain. **Total latency is unchanged at
4096 samples**, since the LF side chain still sets it.

### Fitted antisaturation networks

Spectral skewing is verbatim from §2.6. Antisaturation corners and depths were solved for so
the combined fixed response reproduces §4.1's figures, with the midband held flat and the
corners pulled toward the paper's stated values. Two-pole shelves fit markedly better than
one-pole (worst error 0.27 dB versus 0.77 dB), and the solved corners landed on the paper's
own numbers without being forced there, which is good evidence the model is right:

| Network | Form | Corner | Depth |
|---|---|---|---|
| HF skewing | 2-pole Butterworth LP + shelf | 12 kHz (stated) | 35 dB (stated) |
| LF skewing | 2-pole Butterworth HP + shelf | 40 Hz (stated) | 25 dB (stated) |
| HF antisaturation | 2-pole shelf | **5204.63 Hz** ("about 5 kHz") | **4.8195 dB** |
| LF antisaturation | 2-pole shelf | **95.843 Hz** ("below about 100 Hz") | **1.2832 dB** |

Achieved against §4.1: 25 Hz −9.97 dB (target −10), 5 kHz −1.73 (−2), 10 kHz −5.96 (−6),
15 kHz −10.06 (−10). Worst error **0.27 dB**; midband flat to **0.15 dB**.

### The shelves are what make the decoder realisable

Implementing this turned up a concrete explanation for a remark in §2.6 that otherwise reads
as an aside. The paper says the skewing shelves "do not interfere with the attenuation within
the audio band, but provide phase characteristics that are **essential in the decoding
mode**."

The HF skewing network's shelf zeros sit at 90 kHz, so within the audio band the network is
indistinguishable from a plain two-pole Butterworth low-pass — it measures −5.37 dB at 15 kHz
with the shelf and −5.36 dB without. That makes dropping the zeros look free. It is not: the
bilinear form of *any* plain two-pole low-pass has a double zero at z = −1, so inverting it
puts a double pole on the unit circle and the decoder becomes marginally unstable. Measured
inverse pole radius went to exactly 1.000000 and the null degraded to −171 dB.

So the shelf is not cosmetic — without it the network cannot be inverted at all. The zeros are
kept, and because 90 kHz cannot be represented at audio rates, their digital frequency is
solved for to best match the analog magnitude in band while staying clear of Nyquist. Inverse
pole radii are then 0.476, 0.999, 0.552 and 0.996, all comfortably stable.

## THE CENTRAL CONSTRAINT: per-bin gain notches are not STFT-invertible

This is the most important finding in the project and it bounds what the architecture can
achieve. Measured on a bare pair of WOLAs with a **time-invariant** gain, so control state,
divergence and every other mechanism are excluded:

| Gain profile | Depth | Encode→decode residual |
|---|---|---|
| Smooth tilt across the band | 6 dB | **−150.4 dB** |
| Smooth tilt across the band | 12 dB | **−138.3 dB** |
| Smooth tilt across the band | 24 dB | **−125.9 dB** |
| Single-bin notch | 6 dB | −44.9 dB |
| Single-bin notch | 12 dB | −31.8 dB |
| Single-bin notch | 24 dB | **−15.9 dB** |
| Notch every 32 bins | 24 dB | **+4.0 dB** (error exceeds signal) |

A modified spectrum with abrupt bin-to-bin gain changes is **not the transform of any signal**,
so no decoder can undo it — this is the STFT consistency problem, and it is structural rather
than a tuning matter. A gradual profile is undone to −126 dB even at the full 24 dB of SR boost.

**So the per-bin *units* are fine; what cannot survive is a gain that jumps tens of dB between
neighbouring bins.** That is a real constraint on the "an SR unit per bin" idea, and it is the
thing to design around rather than against.

### What was tried, and what the measurements said

1. **Multiplicative cascade of the two halves instead of summing them.** Required: summing gives
   `v = (I + H + L)u` with H and L at different resolutions, which has no closed-form inverse,
   and iterating diverges because the side chain reaches ~15× the signal at full boost.
   Cascading also matches §2.5's own description. **Kept.**
2. **Control state taken from the source spectrum.** Diverges. The decoder would have to
   reconstruct the source by division, and that small consistency error feeds the control state
   and compounds frame after frame — measured, it dragged the null to **−22 dB**. Fixed by taking
   the state from the **encoded** spectrum, which both sides hold bit-identically: **−41 dB**.
3. **Smoothing the finished gain across bins.** Made things worse: barely helped broadband
   (−48 → −52 dB) and severely harmed tonal material (−77 → **−18 dB**), because smoothing a deep
   single-bin notch in dB spreads it into a wide deep trough — *more* total modification, not
   less. **Rejected**, kept behind a defaulted-off parameter with the numbers recorded.
4. **Critical-band smoothing of the level envelope before the controls see it.** This is the right
   place to intervene: it keeps every bin's own threshold, time constants and weighting, and makes
   the gain profile smooth by construction. Best at **1.0 Bark**. **Kept as the default.**

### The trade-off, now quantified

Two further mechanisms were found and both matter.

**A scaling bug, and a real one.** Bin magnitudes are not on the time-domain amplitude scale:
a tone of amplitude A reads `A · Σw / 2`, some 52 dB of offset at 1024 points. The stage
thresholds are quoted in dB below a full-scale reference, so until that factor was divided out
every threshold sat 52 dB too low and the compressors never released — the low-level boost
measured 15.1 dB instead of 24. Fixed via `Wola64::magnitudeToAmplitude()`.

**Compression accuracy and STFT invertibility trade against each other directly.** Accurate
level detection needs sharp spectral selectivity so a tone is read at its true level;
invertibility needs the opposite. Measured, with everything else held constant:

| Envelope smoothing | Boost at reference (want ~1 dB) | Tonal null | Broadband null |
|---|---|---|---|
| Averaging kernel | **4.97 dB** ✗ | −86.6 dB | −48.6 dB |
| Power, peak-normalised | **0.85 dB** ✓ | −68.0 dB | −41.0 dB |

An averaging kernel dilutes a tone's energy across the critical band, so the compressor
under-reads it and leaves five times the boost the paper allows at reference level. Smoothing
in power with a peak-normalised kernel keeps a tone's own level intact while still spreading a
skirt around it — which is also what a critical-band analysis actually does.

**Fidelity to the papers was chosen**, since that is the point of the project. Every published
figure now lands within about a dB (see below), and the null is recorded rather than hidden.

### Validation against the published curves — all passing

Measured on the fully assembled encoder, end to end:

| Figure | Quantity | Target | Achieved |
|---|---|---|---|
| Fig. 13 | Peak low-level boost | 24 dB | **22.2 dB** |
| Fig. 13 | Tapers at both extremes | required | **9.1 / 22.2 / 12.9 dB** ✓ |
| Fig. 12 | Subthreshold boost at 3 kHz | 24 dB | **22.9 dB** |
| Fig. 12 | Boost at reference level | ~1 dB | **0.85 dB** |
| §4.1 | Antisaturation at 25 Hz | 10 dB | **9.05 dB** |
| §4.1 | Antisaturation at 5 kHz | 2 dB | **1.82 dB** |
| §4.1 | Antisaturation at 10 kHz | 6 dB | **5.61 dB** |
| §4.1 | Antisaturation at 15 kHz | 10 dB | **9.85 dB** |

### Where the null stands, honestly

| Material | Achieved | Plan target |
|---|---|---|
| Steady tone, −20 dB | **−68.0 dB** | −100 dB |
| Broadband at −60 dB | **−60.8 dB** | −100 dB |
| Broadband at −26 dB | **−46.7 dB** | −100 dB |
| Broadband near reference | **−41.0 dB** | −100 dB |
| Loop mode, broadband | **−29 to −34 dB** | −100 dB |

**The acceptance criterion is not met.** The test limits record what is achieved so regressions
stay visible, and the shortfall is marked open in the code and tracked as a task.

The route that could satisfy both sides of the trade-off is to stop applying the gain in the
STFT at all: use the transform for **analysis and control only**, and realise the gain with a
**time-domain filter** whose coefficients come from that analysis — a cascade of shelving
biquads, roughly one per critical band, driven per frame. Encode applies the cascade, decode
applies the exact inverse by swapping numerator and denominator, as the main path already does
successfully at −230 dB. That removes the consistency problem entirely, so sharp level detection
no longer costs invertibility.

It is also structurally what the analog does: control derived from analysis, gain applied by
exactly invertible networks. The reason it was not done first is that it needs roughly 24 biquad
sections per channel per stage with per-frame coefficient interpolation, which is a substantial
piece of work; the STFT route was the cheaper thing to try and it got the process audibly right
and every published curve within a dB, which is what made the trade-off measurable in the first
place.

## Sliding layer composition: decided by measurement

Both compositions were built behind `SlidingComposition` and measured on the paper's
probe-tone setup (a dominant tone over a −70 dB probe floor, as Figs. 14–16 were made).

**Action substitution is not viable, and the reason is structural rather than tuning.**
Eq. (1) means "full boost if *either* circuit achieves it". A masking threshold excludes the
bin it is computed for, so at a dominant's own bin the spread from its neighbours is low, the
sliding circuit reports full boost, and it restores exactly what the fixed band correctly
removed. Measured boost at a 3 kHz dominant: **8.00 dB, i.e. no reduction at all** — the
encoder would overload the channel there. Feeding the sliding circuit `max(own, spread)`
repairs the dominant but then both circuits see the same thing everywhere and the composition
degenerates into a single band.

**Restriction works**: `φ_stage = φ_fixed · σ_slide`. Measured with a 3 kHz dominant —
**0.73 dB** at the dominant, **8.00 dB** recovered far below, **7.99 dB** far above, and a
local dip in between. Both modes remain in the code with the measured outcomes asserted by
test, so a later change that alters either is noticed. `restriction` is the default.

`1 − Π(1−Fᵢ)` is retained for combining several sliding skirts *with each other*, where each
genuinely is a provider. It is the fixed-band-to-sliding-layer combination that must be a
product, because that one is a restriction.

### The skirt asymmetry is the reverse of the circuit, deliberately

The analog's HF sliding band recovers boost *above* a dominant. Measured here, the profile is
broad above and narrow below (8.00 dB at 1.2 kHz, 7.64 dB at 7.5 kHz for a 3 kHz dominant) —
the opposite direction.

That is not an error. The analog recovers boost above a dominant because its fixed band is
*wideband* and has stripped boost across the whole range; there is broad loss to undo. A
per-bin fixed band never creates that loss, so there is nothing to restore. The only remaining
reasons to withhold boost away from a dominant are that masking makes it pointless and that a
smooth profile tolerates channel errors — and masking spreads upward. So the profile follows
the paper's own psychoacoustic statement, that "the only region of the spectrum that is not
boosted in gain is the region that is controlled by masking", rather than an artefact of the
circuit topology it replaced.

### Antisaturation is fitted in the digital domain, per sample rate

Because the skewing low-pass cannot track its analog prototype exactly near the top of the
band (about 1 dB at 15 kHz at 44.1 kHz), the antisaturation parameters are re-solved against
the *digital* cascade at each sample rate, letting them absorb that deviation so the total —
which is what §4.1 actually specifies — comes out right. This is a deterministic coordinate
descent from the offline analog fit, run once per sample-rate change, in
`srNetworks::mainPathBiquadsForRate()`.

**Per-bin fixed bands make MC3/MC6 partly redundant, and that is informative.** In the
analog, MC3 and MC6 exist mainly to stop a *wideband* fixed band from reducing gain in
response to out-of-band energy. A per-bin fixed band only ever sees its own bin, so it is
immune to that by construction. What remains necessary is the *other* half of modulation
control — the level-dependent **end-stop** ("after an appropriate degree of limiting has
taken place, it is unnecessary to continue the limiting when the signal level rises even
further"). So: implement MC as a **per-bin end-stop** for the fixed stack, and keep the
faithful **scalar frequency-weighted** MC1/MC4 for the sliding layer, where a skirt spans
many bins and stop-band immunity genuinely matters.

---

## Signal flow

```
                    ┌─ 8192 overlap-save ─▶ skew · antisat ─────────────┐  (main, z⁰)
                    │                                                    │
  input ─▶ split ───┼─ 4096 WOLA ─▶ LF stages ─▶ iFFT ────────────────── ┤
                    │                 ▲                                  ├─▶ ± ─▶ out
                    └─ 1024 WOLA ─▶ HF stages ─▶ iFFT ─▶ z⁻³⁰⁷² ────────┘
                                      ▲                          + on encode
                    MC1..MC8 ─────────┘                          − on decode

  each stage:  per-bin fixed stack (3 levels HF / 2 LF)  ─┐
               sliding layer (peaks → quasi-2-pole skirts)─┴─▶ 1 − Π(1−Fᵢ) ─▶ gain-slew clamp
```

Channels are processed **independently** — the 1967 paper notes gain-error tolerance
"enables the noise reduction system to operate without control signal interconnections",
and the reference specification quotes crosstalk `<−100 dB` encode/decode. No stereo linking.

---

## Files

All new work lives in `Lime/Source/dsp/`. `PluginProcessor.{h,cpp}` and
`PluginEditor.{h,cpp}` are rewritten from the template.

| File | Responsibility |
|---|---|
| `dsp/Fft64.{h,cpp}` | double real FFT; `vDSP_*D` backend + portable radix-4 fallback behind one interface |
| `dsp/Window.h` | √Hann, WOLA normalisation, COLA verification helper |
| `dsp/Wola64.h` | windowed overlap-add analysis/synthesis engine (side chains); header-only |
| `dsp/OverlapSave64.h` | fixed-filter frequency-domain convolution (main path); header-only |
| `dsp/FilterDesign.{h,cpp}` | analog `H(jω)` → bin tables: skewing, antisaturation, band-defining, MC weightings, A-type bands |
| `dsp/Bark.{h,cpp}` | Bark/ERB edges + masking spread function (−27 dB/Bark down, −6…−15 up) |
| `dsp/ModulationControl.{h,cpp}` | MC1–MC8 scalar weighted sums; per-bin end-stops |
| `dsp/FixedBand.{h,cpp}` | per-bin staggered fixed bands, dual control path + maximum selector |
| `dsp/SlidingLayer.{h,cpp}` | peak picking, corner tracking, quasi-2-pole skirts, `1 − Π(1−Fᵢ)` |
| `dsp/GainSlew.h` | per-frame gain-slew clamp (overshoot-suppression analogue); header-only |
| `dsp/SrStage.h` | one staggered stage = fixed stack + sliding layer + slew clamp; header-only |
| `dsp/SrEngine.{h,cpp}` | top level: main path, LF + HF side chains, encode/decode, alignment |
| `dsp/AtypeEngine.{h,cpp}` | faithful 4-band A-type on the shared dual-path core |
| `dsp/Calibration.{h,cpp}` | calibration noise, calibration tone, Auto Compare, calibration display level |
| `dsp/TapeChannel.{h,cpp}` | loop-mode impairments: noise, saturation, HF loss, modulation noise |
| `dsp/SpectrumProbe.h` | wait-free triple-buffered snapshot of input spectrum + composite curve for the editor; header-only |

**Projucer note — the unity translation unit was tried and abandoned.** The idea was to
`#include` the DSP `.cpp` files into a single `dsp/LimeDsp.cpp` to keep the `.jucer` diff
small. It fails on the classic unity-build hazard: an anonymous namespace is unique per
*translation unit*, not per file, so every file's private `twoPi` and `dbToGain` helpers
landed in the same namespace and collided. Working around it would mean prefixing every
file-local helper, which is a tax on every future file for no benefit. Each `.cpp` now
compiles on its own, which also matches how the standalone tests build them. Registering a
new file is a scripted edit to `Lime.jucer` followed by `Projucer --resave`.

---

## Double precision

- `supportsDoublePrecisionProcessing()` is overridden to return `true`.
- `processBlock (AudioBuffer<double>&, MidiBuffer&)` is the **real** path.
- `processBlock (AudioBuffer<float>&, …)` up-converts into a persistent double scratch
  buffer, calls the same engine, and down-converts — so behaviour is identical whatever
  precision the host chooses.
- All state, coefficient tables, windows and transforms are `double`. No `float` anywhere
  in `dsp/` except the editor's display snapshot.
- `juce::ScopedNoDenormals` is kept in both `processBlock` overloads.
- `Fft64` never falls through to `juce::dsp::FFT`, which is float-only.

---

## Build order

The system was built bottom-up, each layer gated on its tests before the next: the
transforms and reconstruction first (`Fft64` with both backends, `Window`, `Wola64`,
`OverlapSave64`, gated on the perfect-reconstruction checks below); then the
double-precision plumbing and latency reporting around a pass-through engine; the main
path, verified against the §4.1 figures before anything dynamic; the per-bin fixed-band
stack; modulation control; the sliding layer; full staggering with the dual-resolution
split and delay alignment; the gain-slew clamp; decode and the encode→decode null tests;
the modes and tape channel; the A-type process; and finally calibration and the editor
with its hardware-style panel and spectral display.

---

## Verification

Builds used the existing Xcode exporter (`Lime/Builds/MacOSX`), with the DSP checks run as
standalone test targets — the DSP headers do not depend on the plugin wrapper, so they are
exercised directly.

| Check | Target | Achieved |
|---|---|---|
| WOLA perfect reconstruction, side chain zeroed | null to **−150 dB** | **−301 dB** worst case ✅ |
| Overlap-save vs. direct time-domain FIR convolution | agree to **−120 dB** | **−284 dB** ✅ |
| `vDSP` backend vs. portable FFT, random double input | agree to **−280 dB** relative | **−297 dB** ✅ |
| Forward transform vs. Kahan-compensated naive DFT | agree to **−280 dB** | **−307 dB** ✅ |
| Window pair COLA flatness across hop phases | flat to **−280 dB** | **−302 dB** ✅ |
| `Wola64` reported latency vs. best-fit measured delay | exact match | **exact** ✅ |
| STFT engines indifferent to host block size | bit-identical | **verified** (blocks of 1, 7, 37, 511, 1024) ✅ |
| Three aligned paths reduce to a pure delay while side chains are silent | null to **−280 dB** | **exact** — below measurement at 22.05/48/96 kHz ✅ |
| Geometry, align delays and overlap-save legality at every rate tier | exact | **exact** at 11.025–192 kHz ✅ |
| AU loads, renders and passes Apple validation | `auval` succeeds | **AU VALIDATION SUCCEEDED** ✅ |
| AU + VST3 build in Debug and Release | no errors | **clean** ✅ |
| Main path vs §4.1 antisaturation figures, 44.1–192 kHz | within the paper's "about" | **0.27–0.47 dB** ✅ |
| Main path midband flatness | ≤ 0.25 dB | **0.15 dB** ✅ |
| Main path encode→decode null (biquad pair) | −100 dB | **−222 dB** worst case ✅ |
| Inverse cascade stability | all poles inside unit circle | **0.476 / 0.999 / 0.552 / 0.996** ✅ |
| Staggered fixed-band total boost | 24 dB HF, 16 dB LF | **24.00 / 16.00 dB** ✅ |
| Steepest compression ratio (§4.1 "about 2:1") | ~2:1 | **1.77:1** HF, 1.55:1 LF ✅ |
| Boost remaining at reference level | near unity | **0.12 dB** ✅ |
| Encode characteristic monotonic in level | required for invertibility | **monotonic** ✅ |
| Tape channel: all defects disabled | transparent | **bit-identical** ✅ |
| Tape hiss level accuracy | within 1 dB | **0.017 dB** ✅ |
| Tape saturation asymmetry at +6 dB | worse at extremes | **−3.75 dB @ 40 Hz, −0.86 @ 1 kHz, −4.92 @ 15 kHz** ✅ |
| Tape modulation noise | zero in silence, tracks signal | **exactly zero / −54.96 dB vs −55 asked** ✅ |
| Encode→decode null, steady tones and pink noise, channel bypassed | better than **−100 dB** |
| Encode→decode null, broadband transient material | better than **−60 dB** (feedback decode is exact only to within the control smoothing — record the actual figure) |
| Low-level subthreshold curve | matches **Fig. 13**: ~24 dB HF, ~16 dB LF, tapering at both extremes |
| Single-tone encoder characteristic | matches **Fig. 12**: action −62…−5 dB HF, −48…−5 dB LF, ~2:1 between |
| Probe-tone curves, dominants at 200 Hz / 800 Hz / 3 kHz | match **Figs. 14–16** at the stated levels |
| Noise-reduction effect, per Appendix Part 1 procedure | **>20 dB** weighted (CCIR/ARM or A + 20 kHz LP) |
| Reported latency | host null against a plain 8192-sample delay |
| Sample-rate scaling | the Fig. 13 curve is unchanged in Hz at 44.1 / 48 / 96 kHz |
| A-type low-level curve | **10 dB** flat rising to **15 dB** at 15 kHz |

The paper's own figures are the acceptance criteria — the per-bin reinterpretation should
reproduce the published static and probe-tone curves, and any deliberate divergence
(unbounded sliding corners, the slew clamp replacing overshoot suppression) should be
measured and recorded rather than assumed benign.


---

## Final state

**698 checks passing** across six standalone test binaries, plus AU and VST3 building clean
in Debug and Release and `auval` succeeding:

| Suite | Checks | Covers |
|---|---|---|
| `DspTests` | 389 | transforms, WOLA, overlap-save, geometry, main path, fixed band, modulation control, sliding layer, gain slew, engine nulls, published curves, auto gain, the probes following the active direction, the wet low pass, the side-chain detection filters, the per-half time constants |
| `CalibrationTests` | 81 | pink noise, calibration noise and tone, Auto Compare, detector |
| `PhaseTests` | 75 | phase rotator (library component) magnitude, rotation difference, inversion, transparency at zero |
| `TapeTests` | 56 | hiss, saturation asymmetry, HF loss, modulation noise |
| `BarkTests` | 52 | Bark mapping, masking spread, exact all-pairs threshold |
| `AtypeTests` | 45 | four-band A-type, its own nulls, overshoot |

### What is built

`Lime/Source/dsp/` holds the whole process, JUCE-free and double throughout, so every part is
testable standalone:

- `Fft64` — double real FFT, Accelerate plus a portable fallback that agree to −297 dB
- `Window`, `Wola64`, `OverlapSave64` — STFT engines, perfect reconstruction to −301 dB worst case
- `Biquad64`, `FilterDesign` — the fixed networks, inverting exactly (−222 dB worst case)
- `FixedBand`, `SlidingLayer`, `SrStage`, `SrSideChain` — the staggered stages
- `ModulationControl` — MC1 to MC8
- `GainSlew` — the overshoot-suppression analogue
- `Bark` — critical bands and masking spread
- `SrEngine` — encode, decode and loop, with the tape channel between
- `AtypeEngine`, `Calibration`, `TapeChannel`, `SpectrumProbe`

The plugin wraps it with a control set consolidated at v1.0.0 around the dynamics identity —
the four pass trims, drive, auto gain, the masters and mix — and a display of the composite
gain curve over the input spectrum, which is what Figs. 12 and 13 are drawings of.

### Open item — solution built and proven, integration remaining

The broadband encode/decode null (−50 to −80 dB against a −100 dB target) is the one plan
target not met. The fix is no longer a hypothesis: it is implemented, tested and in the
project as `ShelvingBank64`.

**Measured, on a 24-band Bark-spaced peaking bank realising the gain in the time domain:**

| Condition | Residual |
|---|---|
| Static 24 dB profile, encode → decode | **−235.5 dB** |
| **Gains randomised every 256 samples** | **−234.1 dB** |
| Flat request (idle) | −238.4 dB |
| Solved profile vs requested, worst error | 1.46 dB |

That is 190 dB better than the STFT route, and it holds under arbitrary time variation —
the exact case that defeats an STFT.

**Why it works, and the one subtlety that decides it.** A biquad recovers its own input:
`x[n] = (y[n] + a1·y[n-1] + a2·y[n-2] − b1·x[n-1] − b2·x[n-2]) / b0`. That stays exact while
the coefficients change, but *only in direct form I*, whose state is the real history of
inputs and outputs. A transposed direct form II pair cannot do it — its internal states are
expressed in whichever coefficients were in force when they were written, so once the
coefficients move the forward and inverse states stop corresponding. Measured, the transposed
pair diverged to **+46 dB** where direct form I holds at −234 dB. Hence `BiquadDirectFormI64`,
which exists purely for this.

Overlapping peaking filters do not act independently, so the bank builds and inverts its
24×24 interaction matrix once in `prepare()` and solves for the band gains that produce a
requested profile.

**Integration was attempted three times and is unsolved.** The component works; wiring it in
does not, and each attempt failed differently enough to be worth recording:

| Attempt | Arrangement | Result |
|---|---|---|
| 1 | Each side analyses its own input | **+1 to −18 dB** — different state, nothing cancels |
| 2 | Both analyse the encoded signal, signal delayed to match | **−31 to −41 dB** — coefficient sequences offset by the analysis latency |
| 3 | Both analyse the encoded signal, no delay | **+6 dB** — diverges; gain and measured level form a runaway loop |

Together they say the coefficient sequences must match sample for sample *and* the loop from
applied gain back to measured level must be broken. Satisfying both needs the control derived
from something neither side has to infer from the gain-modified signal — per-band envelope
followers sitting outside the gain loop, as the analog has, being the obvious candidate.

`GainRealisation::stft` remains the default, so nothing that works has regressed: the full
suite passes and the published curves are unmoved. `GainRealisation::filterBank` is present,
switchable and documented as experimental-and-worse, so the next attempt starts from measured
ground rather than from scratch.

**The decision integration would carry.** `SrSideChain` would keep its STFT
for *analysis only* — computing the per-bin profile exactly as now — and apply the result
through the bank on a delay-matched time-domain signal instead of multiplying bins. Mechanically
small. But it changes the character of the process: the **applied** gain would have 24 degrees
of freedom rather than one per bin, while the **analysis** stays per-bin. The measurements say
that is forced — a per-bin applied gain is not invertible at any depth worth having — but it is
a real narrowing of the original "an SR unit per bin" idea and is worth deciding deliberately
rather than by default. The STFT path is left as the default so the validated published curves
do not move.


---

## UI: the Aptitude Audio house style

The editor follows Apti-Q rather than inventing a look, so Lime sits beside the rest of the
range. Taken from its `PluginEditor`:

| Element | Value |
|---|---|
| Knob fill | `Colours::blue.darker()` |
| Pointer, captions, values | `whitesmoke` |
| Panel | knob fill `.brighter()` three times — the periwinkle |
| Branding | black, base type size × √2 |
| Toggles | `slategrey` off, knob fill `.brighter()` on, 45 × 45 |
| Meter | `slategrey` at 50 % behind, `green.brighter()` fill, red on clip |
| Grid pitch | 140 px, value box 42 × 21 below the knob, caption below that |

The knob is drawn the same way: a rim ellipse offset +23 and sized `d · 0.66 · knobMod`
behind the body at +20 sized `d · 0.66`, with a rounded pointer struck inward from near the
edge. The modifier constants are the range's own idiom — powers of √2 written through log10,
ratios in the same sense a decibel is: `knobMod ≈ 1.0595`, `alphaMod = 0.5`, `satMod = 2.0`.

Lime differs only where the process does. Apti-Q is sixteen knobs; Lime is four calibration
trimmers, two three-position switches and three toggles, so the freed space goes to the plot
of the encode curve over the input spectrum — which is what Figs. 12 and 13 are drawings of.
The selectors are rows of square toggles rather than menus, because SR/Off/A is a
three-position switch on the real reference hardware and that is the vocabulary the range already uses.

### The flicker, and why it happened

The first version repainted the whole editor at 30 Hz purely to animate the meter, which
redrew every knob, caption and button underneath it thirty times a second. Fixed by making
animation the business of the components that actually change:

- `MeterBar` owns the meter, has its own timer, and repaints **only when the bar would
  move** — a quiet passage costs nothing.
- `SpectrumDisplay` runs at 15 Hz and repaints only when the engine has published a new
  frame, painting its empty state once rather than continuously.
- The editor's `paint()` is now static, just the panel colour.
- The editor's timer runs at 8 Hz and only sets label text and toggle states.

Layout is built entirely by removing slices from the bounds, so no two controls can overlap
and vertical baselines are shared by construction rather than by arithmetic.

The Standalone format was enabled alongside AU and VST3, which is what makes it possible to
run and look at the panel without a host.

### A real bug, not a styling matter

The first working build showed a default-grey panel with the desktop visible *through* the
plot, dark red knobs and invisible captions. That was the **static initialisation order
fiasco**: the palette was written as namespace-scope `inline const juce::Colour` objects
derived from `juce::Colours::blue`, which is itself a global in another translation unit, and
the order between them is undefined. When ours won the race the colours came out garbage or
transparent — and a transparent `background` means `fillAll` paints nothing, hence the
desktop. They are now functions returning by value, in `LimeStyle.h`, which has no
initialisation order to lose.

Two smaller traps worth recording, both found by the compiler rather than by inspection:
`LimeGraphics.h` and `LimeLookAndFeel.h` each needed the other, so the palette moved into
`LimeStyle.h` to break the cycle; and `R"(...)"` raw strings are terminated early by the `)"`
inside an SVG `fill="url(#id)"`, so the artwork uses the `R"SVG(...)SVG"` delimiter.

### Depth: SVG artwork with one light source

The flat first version was solid fills with no light in them. Depth now comes from authored
SVG in `LimeGraphics.h`, with the markup **built at runtime from the palette**, so colour has
one source of truth rather than hex scattered through path fills.

- **Knob** — a cast shadow, a skirt ring lit across its upper left, a domed cap on a radial
  gradient offset toward the light, a contact shadow where cap meets skirt, and a soft
  specular. The pointer is separate artwork, authored pointing up and rotated into place, so
  its own shading stays consistent at every angle.
- **Buttons** — raised when off (lit top edge, shaded bottom, casting a shadow) and genuinely
  *inset* when on, with the shading inverted and a shadow thrown inward from the top edge, so
  an engaged switch reads as pressed into the panel rather than merely recoloured.

No SVG filter elements are used. JUCE's parser covers gradients well but its filter support is
partial, so a `feGaussianBlur` shadow would be silently dropped; soft edges are stacked
gradient stops instead. Artwork is rasterised once per size and cached, since knobs redraw on
every mouse move while the artwork never changes.

Large surfaces are gradients rather than SVG, in `LimeSurfaces.h`: a stretched SVG rectangle
gains nothing over a gradient and costs a rasterisation on every resize. The panel gets a
shallow vertical gradient, a vignette and lit/shaded edges; the plot and the meter sit in
**recessed wells** — darkest along the top and left inner edges where the lip casts into them,
with a highlight along the bottom and right. Getting that inverted is exactly what makes a
recess look like a raised slab, so both wells come from one helper and cannot disagree.

### The panel photographs itself

Nothing in this session could screen-capture the running plugin: `screencapture` fails with
"could not create image from display" and osascript has no accessibility permission, so every
UI change was being made blind, from first principles, with the user's screenshots as the only
feedback. That is a slow loop and it hid two real bugs for two rounds.

So the panel now takes its own photograph. Setting `LIME_SNAPSHOT` to a path makes the editor
write a PNG of itself about a second and a half after it opens:

```
LIME_SNAPSHOT=/tmp/ui.png Lime.app/Contents/MacOS/Lime
```

It goes through `createComponentSnapshot`, which runs the same `paint()` calls the screen does,
into a `juce::Image` — which on macOS is CoreGraphics-backed exactly as the window is. What it
shows is therefore what is on screen, not an approximation of it. Fifteen lines, environment
gated, and it turned an unverifiable part of the project into a checkable one.

### Two JUCE SVG traps, both of which silently render something else

**`Drawable::drawWithin` must not be used to scale parsed SVG.** It fits `getDrawableBounds()`
to the target, and for a parsed SVG that is the union of the *drawn shapes*, not the viewBox.
The knob pointer is deliberately a 3.6 x 22 strip near the top of a 100 x 100 box; `drawWithin`
stretched it to fill 101 x 101, so every knob was a white rounded slab with the body invisible
underneath it. That is why the knobs looked square — twice. `CachedSvg` now scales by the
viewBox parsed out of the markup, so where a shape sits in its box survives.

**Radial gradient radii have to be percentages.** For the default
`gradientUnits="objectBoundingBox"`, JUCE scales `cx`/`cy` by the bounding box but computes the
radius with `getCoordLength (r, width)`, which only scales a value written as a percentage. So
`r="0.78"` means 0.78 *pixels* — a sub-pixel dot that leaves the whole shape flooded with the
gradient's last stop. `r="78%"` is handled correctly. Every gradient coordinate in the artwork
is now a percentage so the two spellings cannot drift apart.

### Lime, in one hue family

The panel is named Lime, so it is lime: one hue family from 76 to 100 degrees, with saturation
and brightness carrying every difference. A panel painted in one flat bright green would be
both garish and unreadable, so the surfaces are dark, the text is light, and the brightest most
saturated lime is spent only on the few things that must be found at a glance — an engaged
switch, the level column, the encode curve.

| Role | Hue | S | B |
|---|---|---|---|
| Knob body | 88 | 0.60 | 0.46 |
| Panel | 96 | 0.34 | 0.21 |
| Plot well | 100 | 0.34 | 0.17 |
| Text | 80 | 0.05 | 0.97 |
| Secondary text | 84 | 0.12 | 0.72 |
| Branding | 84 | 0.62 | 0.74 |
| Switch, off | 94 | 0.22 | 0.34 |
| Switch, on | 82 | 0.80 | 0.62 |
| Accent | 78 | 0.85 | 0.86 |

Branding moved from black to bright lime: black reads as engraved on a pale panel and as absent
on a dark one. An engaged switch takes dark text, because near-white on the brightest lime on
the panel is barely legible.

### Legibility, which was mostly not a colour problem

- **The switches carry words.** "set up", "check tape", "bypass", "encode", "decode", "loop".
  They were "set", "chk", "byp", "enc", "dec" — abbreviations that saved space the panel did not
  need and cost the reader the meaning. Buttons are now sized to their text through
  `widthFor()`, measured in the same font the look-and-feel draws them with, so a label can
  never be wider than the slab under it. `buttonSvg` is authored per pixel size rather than
  stretched, so the corner radius stays circular at any width.
- **The trims read "0.0", not "-0.00".** The default float-parameter formatting gave two
  decimals nobody can set by hand and a minus sign on a value that is zero, because the
  normalised round trip lands a hair below it. `trimAttributes()` rounds first, then clears the
  sign bit, and adds a `+` above zero.
- **The meter is scaled and captioned** — "out" above it, marks every 10 dB, numbers every 20 —
  because a bare coloured column says how much of something without saying what or how much.
  The readout carries its unit: "-87 dB", not "-87".
- **Axis units moved out of the corner.** "dB" and "Hz" sat within a few pixels of each other
  and of the "0" at the bottom right, which read as one crowded smudge.

### Names

No product names anywhere in the interface. The three-position process switch is **Off / A / B**
— A being the four-band process, B the staggered spectral one — and the calibration signals are
"calibration noise at -15 dB" and "calibration tone". The source comments still cite the
published papers by name and section, because that is the design record and what makes the
DSP checkable against its references; it is the *interface* that carries no attribution.

---

## Control smoothing: one pole, 75 ms

Every control reaches the signal path through `ControlSmoother`, a one-pole with a settling time
of 75 ms. A control read once per block and applied as a constant puts a step into the signal
every time it moves, and at a 64-sample block that is one step every 1.3 ms during a knob drag.

**The settling time is quoted as the time to cover 95% of the distance**, which for a one-pole
is three time constants, so tau = 25 ms. Saying tau = 75 ms instead would take 225 ms to arrive,
which no longer tracks a hand on a knob. Measured at 44.1, 48, 96 and 192 kHz: 0.950 of the
distance at 75 ms in every case, monotone throughout, and no overshoot — which matters, because
a gain that rings past its target is a gain that clips.

**A three-position switch cannot be interpolated.** Off, A and B are different processes, not
three values of one thing, and the same goes for encode against decode. So the switch is thrown
while the wet path is silent: `wetAmount` fades out over 75 ms, the configuration changes, and
it fades back. What it fades *to* is a latency-matched dry copy of the input rather than
silence, so the programme is present the whole way across.

That same crossfade is what bypass now is, and it is better than what it replaced. The old
bypass switched the engine's own process to off, which put a discontinuity inside the engine.
Now the engine keeps running its real configuration — coming out of bypass finds the side chain
already tracking the programme instead of starting cold — and at `wetAmount == 0` the output is
the input, delayed, exactly rather than merely equivalently.

Set Up fades its calibration signal in over the programme rather than cutting to it, on the same
one-pole, so engaging it down a live line does not put a step into it.

Applied per sample, with the smoother advanced once per sample and every channel given the same
value, so a stereo pair cannot drift apart mid-ramp. When a smoother has arrived the per-sample
loop is skipped entirely — one multiply by a constant, or nothing at all at unity.

## State: built, validated, installed

**455 checks passing** across the five suites, `auval` succeeding, and AU + VST3 + Standalone
building clean in Release:

| Suite | Checks |
|---|---|
| `DspTests` | 227 |
| `CalibrationTests` | 75 |
| `TapeTests` | 56 |
| `BarkTests` | 52 |
| `AtypeTests` | 45 |

```
AU VALIDATION SUCCEEDED.      aufx Lim1 APTI
```

Both plug-ins are installed and ready to run in a host:

```
~/Library/Audio/Plug-Ins/Components/Lime.component
~/Library/Audio/Plug-Ins/VST3/Lime.vst3
```

Latency is 232 ms at 88.2 kHz, reported to the host, so a DAW compensates it.

---

## The crossover became a control

Decision 14 was "parameters faithful and fixed — only the real hardware's controls, every
internal constant compiled in". The 800 Hz split was one of those constants. It is now a
control, and the reasoning for the change is worth recording because it is a departure.

The split is used in exactly one place: `SrSideChain::prepare` built a table of
power-complementary per-bin weights from it, `|HP|² + |LP|² = 1`. Everything else about the
staggering reads that table. So making it adjustable was not an architectural change, it was a
matter of rebuilding one table and the three per-stage transfer scales derived from it —
`FixedBand::setBandWeight`, `SrStage::setBandWeight`, `SrSideChain::rebuildBandWeights`, all
allocation free, safe to call from the audio thread, and skipped entirely when the value has
not moved.

**Range 320 Hz to 2 kHz, logarithmic, symmetric in ratio about 800 Hz** — sqrt(320 × 2000) =
800 exactly, which is what puts the default straight up on the knob. The first attempt used
200–2000 Hz, whose geometric centre is 632 Hz, so the pointer sat at 60% of the arc at the
default and looked wrong.

The bounds are the DSP's, not taste. The two halves are analysed at different resolutions
deliberately: short window above for time resolution, long window below for frequency
resolution. At 48 kHz the short window's bins are 47 Hz apart, so 320 Hz is about seven bins
and anything much lower stops being resolvable; above 2 kHz the long window is too slow to
track the material it would be handed.

Measured, on the assembled encoder with a tone at −85 dB:

| Crossover | 20 Hz | 500 Hz | 5 kHz |
|---|---|---|---|
| 320 Hz | 3.13 | **21.94** | 20.71 |
| 500 Hz | 3.02 | **20.39** | 20.72 |
| 800 Hz | 2.98 | **18.67** | 20.72 |
| 1300 Hz | 2.96 | **17.45** | 20.70 |
| 2000 Hz | 2.96 | **16.89** | 20.58 |

500 Hz moves 5 dB as the split sweeps past it, monotonically, while 20 Hz and 5 kHz move less
than 0.2 dB — which is the check that the rebuilt weights are still power complementary. All
asserted by test, because a rebuild-in-place into tables the stages already hold pointers to is
exactly the kind of change that compiles and does nothing.

**An encode and its decode must agree on it.** Both directions in one instance are set
together, so Both mode is consistent by construction; across two instances it is on the user,
like level. That is stated in the manual rather than defended in code.

## The crossover control became a resolution control

**The resolution control was itself retired at v1.0.0 — the count is a compiled-in constant
of 2048 sections; see *The sections control retired* under v1.0.0. This section and the ones
under it are kept as the record of the resolution decision and its measurements.**

The control above lasted one version. The section is kept because the reversal only means
anything with the original decision still visible beside it.

The corner is a compiled-in constant again, at the reference system's 800 Hz. What the user
sets now is **how finely that corner is resolved**, and the reasoning came from re-reading why
800 Hz is 800 Hz in the first place. Section 5.6 and Fig. 5.10 of the reference hardware
manual put the split there because the machine has exactly
**two** band-defining networks, which is what could be built: the corner is a consequence of the
count, and the count is a consequence of the hardware. A transform is not limited that way. The
allocation between the two halves can be resolved as finely as the bins allow, so the count is
the parameter with something to say and the corner is a constant of the reference design.

`SrSideChain` therefore samples the same power-complementary single pole onto **N log-spaced
sections spanning 20 Hz to 20 kHz**, and N is the control — a power of two, **4 to 16384**,
default 256. At the bottom of the range the whole allocation between the halves arrives in a few
coarse steps; at the top it is past one section per bin in the long-window half at every rate,
and the allocation is continuous for any purpose that can be measured. Where those two bounds
came from is argued below, and neither of them is where the control started.

Where "one section per bin" lands moves with the rate, because the long window does — up to
the 2x cap on the geometry, which holds the 176.4/192 kHz tiers at the 88.2/96 kHz sizes:

| rate | long window | bins | one section per bin at |
|---|---|---|---|
| 44.1 / 48 kHz | 4096 | 2049 | 2048 |
| 88.2 / 96 kHz | 8192 | 4097 | 4096 |
| 176.4 / 192 kHz | 8192 | 4097 | 4096 |

The list does **not** track the rate, and that is the deliberate choice. The setting has to
mean the same thing in every session — an encode and its decode must agree on it, and a
control whose maximum moved with the rate could not be written down on a tape box. The price
is that the last three positions ask a 48 kHz session for a grid finer than it has bins for.
That is inert rather than wrong, for exactly the reason the second constraint below gives, and
it is measured rather than assumed: at 48 kHz, 16384 sections nulls at −67.5 dB, the same as
2048 to the tenth.

### The knob shows the exponent, not the count

N is always a power of two, so the parameter is the **exponent**: **Exp**, an integer from
**2 to 14**, with N = 2^Exp. The caption on the panel is `exp`, lower case, matching `send in`
and `mix`.

| exp | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sections | 4 | 8 | 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |

Two reasons, and neither is cosmetic. **2 to 14 is a scale someone can hold in their head and
4 to 16384 is not** — a reader who sees "9" knows immediately where on the arc it sits, where
"512" has to be counted out against 16384 before it means anything. And the interesting part of
the range is **evenly spread across the exponent while it is bunched at one end of the count**:
everything the measurements below turn on happens between 4 and 256 sections — the first seven
of the thirteen positions, but a small fraction of the count. A control whose useful travel is
all in the first few degrees of its arc is a badly scaled control.

The manual still quotes the count, because the count is what the process does. The exponent is
how the setting is named, not what it means.

### Where the two bounds came from

The range began as 1 to 8192 and neither end survived.

**The floor is 4 because the two positions below it were both traps, in opposite directions.**
At 1 section there is a single weight across the whole band and no frequency-dependent
allocation at all, and it nulled at **−67.5 dB, about 13 dB better than anything else measured
at the time**. It scored best because it does nothing: a flat weight is the smoothest profile
there is, with no bin-to-bin step anywhere for the decode to fail on, which is the central
constraint of this project seen from the other end. A position that wins the measurement by not
working is worse than useless on a control someone is going to set by looking at numbers. The
smoothed opposition has since raised the plateau itself to −67.5 dB, so the flat weight would
no longer even win — which settles the removal rather than reopening it. At 2 sections the
opposite: it is the **worst** setting measured, −24.0 dB, because that is where the entire 8 dB
of allocation arrives as one step.

Dropping 2 has a real cost and it is not going to be glossed. **2 sections was the reference
machine's own arrangement** — one crossover at 800 Hz, three stages above and two below — and it
is now unreachable. In a project whose whole premise is a faithful reinterpretation of that
machine, removing its exact configuration is a loss, and the honest statement of it is that the
configuration is no longer offered because it is a *broken* round trip rather than a *period* one:
−24.0 dB, against −29.9 for the coarsest setting that remains as both measured at the time, and
−67.5 for the plateau today. That is
a judgement about what a control should let someone choose, not a measurement, and it should be
reopened if the round-trip null is ever closed.

**The ceiling is 16384 for a panel reason, and it is labelled as one.** 16384 is past one section
per bin at every rate the plugin runs at, so the top position does nothing the position below it
does not already do. It exists to make the number of positions **odd** — thirteen rather than
twelve — so that the default sits straight up on the knob. That is the same motivation as the
skewed taper that put 800 Hz at twelve o'clock on the old crossover control, recorded under
*The crossover reaches 200 Hz, and the taper is skewed to pay for it*. It costs nothing because
the position is inert rather than wrong, and it is called a panel decision here so that nobody
later reads it as a claim about the process.

### Two constraints that are load-bearing

**The section grid is defined in hertz, not in bins.** The two halves run at different
resolutions on purpose, so a grid in bins would quantise the profile differently in each of
them and `|HP|² + |LP|² = 1` would stop being exact — the one identity the whole cascade rests
on. Defined in hertz, both halves evaluate the identical function of frequency and the sum
stays exact whatever the count. Measured: the halves stay power-complementary to **−319 dB** at
every count.

**The count is not clamped against a half's own bin count.** The temptation is obvious —
sections finer than a half can resolve look like waste. They are not waste, they are free: the
extra boundaries fall inside a bin and change nothing. Clamping per half would hand the two
halves different grids, which is the hertz-versus-bins failure arriving by another route and
breaking the same complementarity. A saved loop iteration is not worth that.

### Measured

Encode/decode null against the section count, separate encoder and decoder instances, broadband
noise at −26 dB:

| sections | 4 | 8 | 64 | 128 | 256 | 2048 | 16384 |
|---|---|---|---|---|---|---|---|
| exp | 2 | 3 | 6 | 7 | **8** | 11 | 14 |
| null | −30.0 dB | −39.1 dB | −64.0 dB | −66.9 dB | **−67.2 dB** | −67.5 dB | −67.5 dB |

**The null improves as the grid gets finer and plateaus at the 256-section default.** Past that
the 8/N dB steps the staircase puts into the encode curve are already below whatever else is
limiting the null, so the upper half of the range is free — while 64 sections, which once sat on
the plateau, is now 3.5 dB short of the finest grid. The finding that shapes the bottom of the
range is at the other end: coarse settings are **not** a cheaper approximation of fine ones. The
coarsest reachable setting costs about **37.5 dB of null** against the plateau, which is audible.
The two settings coarser than it are no longer offered, and what they measured is the argument
for that — see *Where the two bounds came from* above.

**The default is 256 sections, exp 8.** When it was chosen the measurement on its own did not
get there: the null as then measured argued for 128, the first power of two clear of the knee,
and what moved it one detent further was the odd-position requirement — thirteen positions to
put the default straight up on the knob — a panel argument, recorded as one. On the smoothed
system the null itself arrives at the same detent: 128 sits 0.6 dB short of the plateau and 256
begins it, so the panel argument and the measurement now agree.

Low-level encode boost, tone at −85 dB, through the assembled encoder:

| sections | 20 Hz | 500 Hz | 5 kHz |
|---|---|---|---|
| 4 | 2.53 | 16.77 | 19.92 |
| 8 | 2.51 | 17.63 | 19.94 |
| 64 | 2.53 | **18.80** | 19.94 |
| 128 | 2.53 | 18.38 | 19.95 |
| 256 | 2.53 | **18.22** | 19.95 |
| 2048 | 2.53 | 18.23 | 19.95 |
| 16384 | 2.53 | 18.22 | 19.95 |

**The default traces the continuous profile exactly, and that is a coincidence.** At 256
sections the boost at 500 Hz reads 18.22 dB against the finest grid's 18.22 — **a match to the
hundredth**, where 128 reads 18.38, 0.16 dB out. So the detent the panel wanted happens also to
be the detent that reproduces the reference characteristic. It is flagged as luck rather than left to look like a derivation, because a
reader who took it for one would conclude the two arguments agree and would not notice if a
later change to either broke that.

500 Hz is **not** monotone in the count, unlike the sweep in the section above where the corner
itself was moving. Nothing is wrong with it: what 500 Hz reads depends on where it falls inside
its own section's span, and that shifts every time the grid is redrawn, so the figure walks
around within a couple of decibels instead of travelling in one direction. It is recorded here
because a reader meeting the 64-section row cold would take it for a fault. 20 Hz and 5 kHz move
by a few **hundredths of a decibel** across the entire range, which is the check that the quantising did not
disturb the complementarity.

### A wrong assertion about the flat-profile case, caught by measurement

While the control still reached one section, the first version of the test for that row asserted
that the boost at 20 Hz and at 5 kHz would land on the **same value**, on the reasoning that with
one section the side-chain weight is flat, so every frequency must receive the same treatment.
**That was wrong, and the measurement said so.**

What goes flat at one section is the *side chain's allocation between the two halves*, not the
encode curve. The fixed main path is still there and still frequency-dependent — spectral
skewing and antisaturation, the same networks that make the plotted curve taper at both
extremes — and it accounts for the rest of the spread. Measured, the 20 Hz-to-5 kHz spread
narrows from 17.7 dB at the finest grid to 11.7 dB at one section: **the side chain's own 8 dB
collapsing, and nothing else moving.**

The setting has since gone from the control, but the note stays. The mistake is the instructive
kind: it is easy to write "the side chain is flat" and finish the sentence as "so the response is
flat", and the assembled encoder is not the side chain. The same slip is available anywhere in
this project where a component's behaviour gets quoted as the system's.

### Moving the control while audio runs

It is a stepped control, so every move is a jump, and a per-bin gain changing between frames is
exactly what `GainSlew` exists to limit. `testSectionTransitionIsSmooth` measured a bare live
move rather than assuming anything about it: the largest step in output level between
consecutive hops, on a −60 dB tone, while the parameter moved.

| move | worst hop-to-hop step |
|---|---|
| coarsest to finest | 1.4 dB |
| finest to coarsest | 1.7 dB |
| one detent, 128 to 256 | 0.6 dB |

Every step stayed inside the slew's own limits — a plain per-frame clamp, 3 dB engaging and
6 dB releasing, relaxed where modulation control calls for it — so a bare move did not click.
What the measurement could not make true is that a bare move is *correct*: the section count
reshapes the encode characteristic itself, and an encode and its decode must agree on that
characteristic, so a decoder running across the move would be undoing a curve the encoder had
already abandoned. The design as it stood at release preparation therefore split the routing
by mode. In loop mode both directions live in one instance and move together, so the agreement
held through the move and the control applied live, slew-limited, on the strength of the
measurement above. In encode and decode modes the other half of the pair is a separate
instance that cannot follow a mid-stream change, so there the control routed through the same
1500 ms fade-to-dry transition as Process and Mode: the wet path fades out, the table is
rebuilt while it is inaudible, and the process fades back. The measurement stands as evidence
that the rebuild itself is benign — a control that rebuilds a weight table the stages already
hold pointers to is precisely the kind of thing that could click without anyone noticing until
a session. The routing, and the control it routed, were retired before release — see *The
sections control retired* under v1.0.0.

### Panel

The knob keeps the slot the crossover had, in the middle of the top row between the sends and
the returns, and is stepped through thirteen detents reading "2" to "14" — the exponent, not the
count — instead of running continuously. Thirteen is odd so that the default is the middle
detent and the pointer at the default points straight up; that is what the ceiling of 16384 is
for. The single lime line the plot drew at the crossover is gone; in its place are faint section
boundary ticks, **drawn only when N ≤ 64**, which is exp 6 and below. Denser than that the
boundaries fall closer together than the pixels between them, and what the eye gets is an even
accent wash over the whole plot — which tells a reader less than drawing nothing would. With the
default at 256 that means **the default draws no ticks at all**, which is the honest answer: at
that density there is no structure on the plot to read off.

## An apparent discontinuity in the process that was one in the plot

The user saw a step of about 12 dB at 800 Hz in the encode curve and asked whether the response
could be splined smooth. Measured through the assembled encoder first, rather than smoothed:

| Hz | 700 | 750 | **800** | 850 | 900 | 1000 |
|---|---|---|---|---|---|---|
| boost | 19.88 | 20.14 | **20.38** | 20.61 | 20.82 | 21.18 dB |

Smooth and monotone. The step was the display's: `gather()` spliced the long window's probe
below the crossover and the short window's above, and **at the crossover each half on its own is
only doing part of the work**, so showing one of them there shows about half the boost.

Fixed by combining rather than choosing. The halves are cascaded in the engine, so their boosts
multiply and their decibels add; the input levels are power summed, since each half sees the
programme through its own band split. Built on a 480-point log grid so two different analysis
resolutions land on one continuous trace.

Interpolating across the join would have drawn a smooth curve through the wrong values, and
nobody looking at it could have told. **Measure before smoothing** is the lesson, and it is the
second time in this project that a plausible cosmetic fix would have hidden a real answer.

The same pass found that the plot was showing only the side chain's boost, so it read a flat
24 dB out to 20 kHz where the assembled encoder delivers about 13 dB at 15 kHz and 9 dB at
30 Hz. The fixed main path — skewing and antisaturation — is now evaluated per sample rate from
`srNetworks::mainPathBiquadsForRate` and added in, so the curve on the panel is the encode
response rather than a component of it, and it tapers at both extremes as Fig. 13 does.

## Check Source, which did nothing

`checkTapeParam` was read into a member and never used. The control existed on the panel and
had no effect — worse than not having it.

On the hardware it is a monitor switch: hear the return *undecoded*, so a problem on the
channel can be told from a problem in the decoding. Implemented as exactly that:
`SrEngine::setCheckTape` suppresses the decode pass wherever there would be one, and the
latency padding is now computed from how many passes were actually taken rather than from the
mode, so nothing shifts in time:

```
passesTaken = (encode ? 1 : 0) + (decode ? 1 : 0)
padding     = (2 - passesTaken) * passLatencySamples
```

Which reproduces every previous case exactly and adds the two new ones. Renamed to **check
source**, since nothing else on the panel says tape. The parameter identifier stays
`checkTape`: it is invisible, and changing it would silently drop the setting out of every
session already saved. Same for the trims, whose identifiers remain `recIn`/`playOut` while
their names became **send in** and **return out**.

## Three attempts at telling the user what a control does

1. **A sentence across the panel** saying what the process was currently doing. Removed: it
   took a whole row, it changed under the reader, and it answered a question nobody was asking
   at the moment they were looking at a knob. It also carried the one non-ASCII string literal
   in the panel, which rendered as `â€` — `juce::String (const char*)` treats its input as
   ASCII, not UTF-8, so an em dash written as `—` arrives as three Latin-1 characters.
   Worth knowing: `operator+` does not have this problem, which is why the same dash rendered
   correctly two lines away.
2. **Tooltips.** Removed at the user's request. They only appear for someone already hovering
   the thing they wanted explained.
3. **A manual.** `MANUAL.md`, rendered to `docs/Lime-Manual.pdf` with pandoc and XeLaTeX,
   twelve pages, with the panel photographed by the plugin itself. What a control does can be
   read once instead of rediscovered on every hover. The panel now carries no explanatory text
   at all — only captions.

## State

**460 checks passing**, `auval` succeeding, AU + VST3 + Standalone clean in Release:

| Suite | Checks |
|---|---|
| `DspTests` | 232 |
| `CalibrationTests` | 75 |
| `TapeTests` | 56 |
| `BarkTests` | 52 |
| `AtypeTests` | 45 |

## Set Up left the panel

Asked whether Set Up and Check Source earn their place in a DAW-only plugin, the honest answer
was: Check Source does, Set Up does not.

**Check Source stays** because in Both mode the channel is inside the instance and there is
otherwise no way to hear what it did before the decoder cleans it up. Nothing else on the panel
can do that.

**Set Up is off the panel.** It exists to align an unknown channel — send calibration noise,
listen to what comes back, trim until they match. Across a bit-transparent path between two
plugin instances there is nothing to align, so it is a switch for a case a DAW user may never
hit.

The parameter and the entire calibration path remain: the noise and tone generators, the nick
detector, Auto Compare, all 75 checks of `CalibrationTests`. It is reachable from a host's
parameter list, and the Auto Compare indicator still lights if something engages it. So the
capability is intact and only the panel space is reclaimed — which is the right trade when the
alternative is a control most users would have to be told to ignore.

---

## v0.1.1 — cost, controls, and the states the panel was not admitting to

### The CPU, and where it actually was

One stereo instance in `both` at 44.1 kHz cost about 40% of a core. An audit of the whole
per-hop path found the FFTs were not the problem — they are already Accelerate `vDSP_*D`, with
their tables built once — and neither was the build, which was already `-O3` with LTO and
universal arm64. The cost was in two places nobody had looked at:

- **around 35 million scalar `libm` calls a second** in the per-bin control loops, and
- **about 270 MB/s of `memset`**, from clearing the Li Chao envelope tree once per stage per
  frame.

Both were fixed without changing a single measured number. The whole of `DspTests` came out
byte-identical after each step.

| Change | What it was |
|---|---|
| `CompressionLaw::shapeRatio` | `std::pow (ratio, kneeExponent)` with the exponent a member rather than a literal, so a full runtime `pow` ran per bin per stage — and, worse, made the frame loops opaque to the vectoriser |
| `SrSideChain::rebuildEnvelopeCoefficients` | two `std::exp` and a division per bin, in each of two passes, for a value that only changes when a parameter does |
| `SpreadThreshold::generation` | 131 KB of `Line` cleared per LF stage per frame, when an insert only ever touches O(log n) nodes. Version-stamped instead |
| bin magnitude | `std::abs` on a complex dispatches to `hypot`, whose overflow guards are pointless on an audio spectrum |
| `GainSlew` | a `log10`/`pow` round trip per bin, on a limiter deliberately inert for most material. State kept linear; only a bin that actually trips converts |
| `SpectrumProbe::setActive` | two `log10` per bin published whether or not an editor existed. The header claimed this was free. It was not |

That took 40% to 35%. The rest was structural.

### The LF half was analysed sixteen times over

`SrGeometry` gave both halves one hop of 256. Against a 4096-point window that is sixteen-fold
overlap where the short window ran four — and the LF half carries 2049 of the 2562 bins, so it
was about three quarters of the engine's cost. The redundancy bought nothing: a 4096-point
window already smears across 85 ms, so the extra frames cannot resolve anything the window has
already blurred.

`WolaWindow` is sqrt-Hann and COLA-exact for any hop dividing `size/2`, so the LF half now runs
**hop 1024**. It is still a whole multiple of the HF hop, so the two halves' frames land
together and the frame-alignment invariant holds; latency is unchanged. `GainSlewParams` are
quoted per frame and are scaled for the longer one, so both halves still permit the same rate
of change in dB per second.

**Result: 40.3% of a core to 15.5%, a factor of 2.6.** Every null figure held:

| | before | after |
|---|---|---|
| engine-off pure delay | −999 dB | −999 dB |
| loop-mode null, 44.1 / 48 / 96 kHz | −28.8 / −29.9 / −34.0 dB | unchanged |
| main path forward-then-inverse | −230.6 dB | unchanged |
| up/down round trip, steady tone | −68.0 dB | −67.7 dB |

### 32-bit floats were considered and not adopted

Float would double the SIMD width and open faster vForce paths, and the side chain arguably
does not need double — it computes a control signal, and complementarity needs the encoder and
decoder to agree, not to be accurate. But the engine's contract is exact double behaviour
measured in the −280 dB range, and after the two stages above there was no longer a cost
problem to justify the risk. It stays 64-bit throughout.

### Switching processes replayed a quarter of a second of old audio

`DelayLine64::process` returns without writing anything at zero delay, and the engine runs at
zero padding in loop mode. So switching from B to A, which takes the padding from 0 to the full
reported latency, had the line answer with whatever it had last been holding — B's audio, from
whenever the delay was last in use, arriving underneath A.

Changing the delay now discards the line. Separately, a change of process or mode resets the
engine outright: the overlap-add buffers, the side-chain histories and the tape state all
belong to the configuration that is being left, and the switch happens while the wet path is
silent, so clearing them costs nothing audible.

### 800 ms transitions, and why not a one-pole

The mechanism was right — fade out to the latency-matched dry, throw the switch while it is
silent, fade back — and the curve was wrong. A one-pole is steepest at the moment it starts,
which is what makes it feel responsive under a knob and exactly what makes a switch sound like
a switch. It also never arrives, so the code had to wait for `> 1.0e-4` and accept whatever
time that took.

`lime::Crossfade` is a finite raised-cosine ramp: zero slope at both ends, and 400 ms means
400 ms. Bypass keeps its own 75 ms one-pole, because a bypass that takes 800 ms to answer is
broken for the job people press it for.

### Check Source went after all

The earlier entry above argued Check Source earned its place. Using it disagreed. In **down**,
which is where a DAW user reaches for it, suppressing the only pass in circuit leaves delay and
nothing else — it is bypass with extra steps. It was only ever distinct in **both**, and a panel
carrying two switches that do the same thing in the common case is worse than one without it.

### Return trims that did nothing

In **both**, the trim pair was chosen by `mode == decode`, so the send pair was live and the
return pair was inert — and `send out` was applied after the decoder, when on the hardware it
sits before the channel. `SrEngine::setLoopTrims` puts the three that fall inside the loop at
their real positions: after the encoder, after the modelled channel, and after the decoder.

### New controls

- **input** and **output**, ±24 dB, outside everything including bypass.
- **mix**, 0 to 200%. Past 100 the same line extrapolates rather than interpolates, so the
  difference between processed and dry is applied again — the quickest way to hear what a
  subtle setting is doing.
- **phase**, 0 to 180°, on the processed path before the mix. A true rotation via a 90-degree
  phase-difference network: flat magnitude, the same shift at every frequency. Bypassed
  outright at zero, because an allpass in circuit on an encode that the decode does not undo
  would break the round trip even at "no rotation".
- **delay**, 0 to 40 ms, same path, third-order Lagrange so hundredths of a millisecond mean
  something. A single utility line — no feedback.
- **polarity**, after the mix, folded into the output gain's sign so it ramps through zero.
- **crossover** is continuous again. The quarter-hertz dead zone in `setCrossoverHz` existed
  because callers passed whole hertz; with a continuous control it would hold the engine a
  fraction away from what the panel prints, which is the bug the dead zone was introduced to
  fix in the first place.

### What the plot was not admitting

The trace kept animating under bypass, and kept showing a stale B curve with the process on
**off** or **a** — in both cases the probes are not running at all, because `SrEngine::process`
pays the latency as pure delay and never reaches a side chain. It now says so: a flat line at
0 dB for **off**, and a frozen, desaturated last frame with a caption for bypass and for **a**.
Building a real curve for **a** from `AtypeEngine`'s four band gains is outstanding work.

### Panel

Ten knobs in two rows of five — the hardware's controls on top, the plugin's own below in
signal order — with every switch in one justified grid to the right, sized to fill the height
of both rows. The footer's "64-bit / latency / sample rate" readout is gone: the precision
never changed, the rate is the session's, and the latency is reported to the host and
compensated automatically.

### Panel, after using it

The switch block went from three ragged rows sized to their own words to a justified grid: each
row spans the column exactly, so every switch shares one of three widths and the rows divide the
full height of the two knob rows beside them. Legends are set in capitals at a root of two above
the old size, less 7% — capitals have no descenders and fill more of a line than the same size
in lower case, so the larger figure crowded the slab.

The crossover moved to the middle of the top row rather than the end of it. That puts the sends
to its left and the returns to its right, and lands it above the delay, so the two controls that
are neither a level nor a switch share a column.

It reads to a tenth of a hertz rather than a hundredth. A hundredth is far below anything the
side chain resolves and cost the readout a character it had to be widened for.

### The crossover reaches 200 Hz, and the taper is skewed to pay for it

**Superseded at v0.2.0 and kept for the record only.** The range, the
200 Hz floor and the two-rate taper are all gone with the control they belonged to. The corner
is fixed at 800 Hz again and the knob in that slot now sets how finely the split is resolved.
See *The crossover control became a resolution control*.

The floor was 320 Hz, chosen because at 48 kHz the short window's bins are 47 Hz apart and 320 Hz
is about seven of them. It is 200 Hz now, which is four. That is a real loss of resolution in the
side chain and it is the user's trade to make rather than one the control makes for them; the
process still works throughout the range.

The cost lands on the panel rather than in the DSP. The old range was symmetric in ratio about
its own geometric mean — sqrt(320 * 2000) = 800 — which put the reference system's own split
straight up on the knob for free. sqrt(200 * 2000) is 632, so a plain log taper would have left
the default noticeably left of centre. The taper is now logarithmic at two different rates: the
lower half of the arc covers 200 to 800 and the upper half 800 to 2000. The knob is coarser
below the default than above it, which is the right way round — that is where the range was
extended into, and where the settings are further apart in effect.

### The window is a fixed size

It was resizable, and the handle was a lie. `resized()` lays out in pixels — it removes slices
of a fixed height and steps the knobs along a fixed grid pitch — so dragging the corner never
scaled anything. Made larger the panel grew empty margin; made smaller the knobs silently
dropped off the right-hand end one at a time, because each is placed only
`if (row.getWidth() >= gridPitch)`. A handle that quietly deletes controls is worse than no
handle.

Scaling properly means every metric in LimeStyle becoming a ratio and the artwork
re-rasterising per size, which is real work for a panel whose one job is to be read. Until
that is done the window is the size it was designed at.

### Three things first use turned up

**The plot said "no signal" against a live engine.** The editor took its probe pointers once,
in its constructor, from the engine's per-channel side chains — and `SrEngine::prepare`
rebuilds those from scratch. So the pointers were null if the host built the editor before it
prepared the processor, and dangling if it prepared again afterwards. They are re-fetched on
every timer tick now, which is idempotent and costs nothing.

The symptom was diagnostic in itself: switching the process to **off** made the trace appear,
and switching back to **b** made it vanish again. The off state draws a flat line without
consulting a probe at all, so it was the one state that could not fail.

**`off` in loop mode skipped the modelled channel.** `SrEngine::process` early-returned to a
pure delay for anything that was not `spectralRecording`, and that return came before
`tape.process`. So in **both** — the mode whose entire purpose is to put the channel inside the
instance — `off` was a clean delayed signal and `b` was a clean processed signal, and the
comparison the panel is built around was inaudible. Measured on a 220 Hz burst at 48 kHz, the
gap now sits at −85.0 dB with the process off and −98.5 dB with **b**.

`DspTests` gained a check for both halves of this: encode with the process off is still an
exact delay to −999 dB, and loop with the process off is audibly not one. The existing
pure-delay test now names the mode it means rather than relying on the default.

Worth recording as still outstanding: **a** does not reduce noise in **both**. `AtypeEngine`
has one instance per plugin and the wrapper selects encode whenever the mode is not decode, so
in loop mode A encodes and never decodes. Completing it means a second A-type instance with the
channel between them.

**Mono was checked rather than assumed.** 8.49% of a core against stereo's 17.01% at 48 kHz —
exactly 2.00x, so a channel that is not there costs nothing. `SrEngine::prepare` allocates side
chains per channel and every loop is bounded by `usable`, so there was nowhere for waste to
hide, but it is now measured rather than argued.

**Transitions are 1500 ms**, up from 800.

### The plot's captions were mojibake

"process a — no spectral analysis" drew as "process a â€ no spectral analysis". Not a typo in
the source: the bytes were correct UTF-8. `juce::String`'s `const char*` constructor reads its
input as `CharPointer_ASCII`, one character per byte, so the em dash's three bytes became three
characters. Anything non-ASCII drawn on this panel has to go through `String::fromUTF8`, which
is what the degree sign on the phase readout already does — and why that one was fine.

Both captions are plain ASCII now. Punctuation that needs no escaping is the safer answer for
a caption than typography that depends on remembering a constructor's encoding rule.

---

## v0.2.0 — how finely the split is resolved

### The crossover control became a resolution control

The one substantive change of this version, recorded in full beside the decision it reverses —
see *The crossover control became a resolution control* above, and the section before it for
what was reversed. In short: the corner is a compiled-in constant again at 800 Hz, and the
control in that slot now sets how many log-spaced sections the power-complementary pair is
sampled onto, **4 to 16384 in powers of two, default 256**. The encode/decode null improves as
the grid gets finer and then plateaus, so the top of the range is free; the coarsest reachable
setting costs tens of decibels of null, so the bottom of it is not a cheaper approximation of
the default but a different and worse round trip.

The grid is in hertz rather than bins and is not clamped against either half's bin count.
Both of those are complementarity constraints rather than conveniences, and both are argued
where they were decided.

### The parameter is the exponent

The count is always a power of two, so the host-visible parameter is **Exp**, an integer
**2 to 14** giving 4 to 16384 sections, captioned `exp` on the panel in the lower case the other
captions use. 2 to 14 is a scale a reader can hold in their head where 4 to 16384 is not, and the
part of the range the measurements turn on is spread evenly across the exponent where on the
count it is bunched at one end. The manual still quotes the count, because the count is what the
process does.

### Both ends of the range moved, and only one end moved for a measured reason

It began as 1 to 8192 with the maximum as the default. All three of those numbers changed.

**The floor went to 4**, because the two positions below it were traps in opposite directions.
One section is a flat weight, no frequency-dependent allocation at all, and it nulled at
**−67.5 dB — about 13 dB better than anything else then measured** — for the reason that it does
nothing: a flat profile has no step anywhere for the decode to fail on, which is the central
constraint of this project seen from the other end. Two sections is the reverse, the **worst**
setting measured at −24.0 dB, because the whole 8 dB of allocation arrives as one step. A
control offering a position that scores best by not working, one detent from the position that
scores worst, is a control that misleads whoever reads its numbers.

**Losing 2 costs something real and it is not glossed anywhere.** Two sections was the reference
machine's own arrangement, and it is now unreachable. Removing the exact configuration of the
machine this project reinterprets is a genuine loss of fidelity; what settled it is that the
setting is a broken round trip rather than a period-correct one, at −24.0 dB against −29.9 for
the coarsest setting that remains. That is a judgement, not a measurement, and it should be
reopened if the round-trip null is ever closed.

**The ceiling went to 16384 for a panel reason**, stated as one. It is past one section per bin
at every rate the plugin runs at, so it does nothing the position below it does not. It exists to
make the position count **odd** — thirteen — so the default sits straight up on the knob, which
is the same motivation as the skewed taper the old crossover control used to put 800 Hz at twelve
o'clock.

### The default moved off the maximum, and then off the measurement

**256 sections, exp 8.** The first version of this control defaulted to the top of the range on
the grounds that the top was free. It is free, but so was most of the range: the null as then
measured plateaued at about 64 sections, and 128 was the first power of two clear of that knee,
so the measurement argued for **128**. What put the default at 256 was the odd-position
requirement above — a panel argument, not a process one.

And 256 then turned out to trace the reference characteristic exactly: 18.22 dB of boost
at 500 Hz against 18.22 dB at the finest grid, **a match to the hundredth**, where 128 reads
18.38 and is 0.16 dB out. **That is a coincidence and is flagged as one.** Two independent arguments landing
on the same detent is the kind of thing that reads afterwards as a derivation, and if either
argument moves, nobody who believed that would notice the agreement had gone.

### Transitions are measured now

`testSectionTransitionIsSmooth` checks that moving the control while audio runs does not click.
Worst step in output level between consecutive hops, on a −60 dB tone: **1.4 dB** coarsest to
finest, **1.7 dB** finest to coarsest, **0.6 dB** for a single detent from 128 to 256. The
existing gain slew already covers a parameter move as well as it covers programme material, so
nothing new was added — but it was verified rather than assumed, which for a control that
rebuilds a live weight table is the difference between knowing and hoping.

### A wrong assertion caught by measurement

Recorded in full above under *A wrong assertion about the flat-profile case*. In brief: the first
test written for the one-section case asserted that the boost at 20 Hz and at 5 kHz would land on
the same value, since the side-chain weight there is flat. Wrong — what goes flat is the side
chain's allocation, not the encode curve, and the fixed main path is still frequency-dependent.
Measured, the 20 Hz-to-5 kHz spread narrows from 17.7 dB at the finest grid to 11.7 dB at one
section: the side chain's own 8 dB collapsing and nothing else moving. The setting is gone from
the control now; the note is kept because the slip generalises.

### Panel

The knob is stepped through thirteen detents reading "2" to "14" and keeps the slot the
crossover had. Thirteen is odd so the default is the middle detent. The lime line the plot drew
at 800 Hz went with the control that fed it; faint section boundary ticks are there instead,
drawn only up to 64 sections — exp 6 and below — because past that they read as a wash rather
than as structure. With the default at 256 that means the default draws no ticks.

### JUCE 8.0.14 to JUCE 9.0.0

Done inside this version rather than deferred, because the toolchain moved under it and a
plugin that does not compile against the installed JUCE is not in a releasable state.

`BREAKING_CHANGES.md` lists five entries for 9.0.0, and because the project was on **8.0.14** —
the last of the 8.x line — those five are the whole gap rather than an accumulation across
versions. Two are Windows multi-touch and Linux EGL, neither of which a macOS AU/VST3 meets.
Three concern `Drawable` and SVG, and **exactly one of them reaches this code**:
`Drawable::createFromSVG (const XmlElement&)` is gone, because SVG support is now the lunasvg
library and lunasvg does its own XML parsing. `LimeGraphics::parse` hands it the string
directly through `createFromSVGString` instead, which is one step shorter than what it
replaced — `juce::parseXML` had no other use in the project and went with it.

Checked and found not to apply, so that nothing was changed on suspicion: `Drawable::draw` with
an opacity and a transform still exists, so `CachedSvg::render` is untouched; `Drawable` no
longer inheriting `Component` is irrelevant here because `CachedSvg` owns its drawable and
rasterises it rather than parenting it; and the `DrawableShape` stroke-type signature changes
are to a class this project does not use. The Projucer resave was mandatory regardless —
`JuceHeader.h` carries a version guard that `#error`s when the generated glue is older than the
modules — and it brought in eight new translation units, one of them `lunasvg`.

**The parser swap changed how the knobs look, and the change is an improvement.** The whole
panel is SVG authored as a string at runtime, so a new parser is a real risk to the artwork and
not only to the build. Measured by photographing the panel before and after at the same size and
settings and differencing it:

| region | drawn with | max channel delta |
|---|---|---|
| plot well, meter, background | JUCE `ColourGradient` | **0** — bit-identical |
| buttons | SVG | 2 |
| knobs | SVG | **79** |

All the knob difference is in the skirt annulus between the dome edge and the outer rim, and
all of it on the upper-left side:

| skirt colour, by angle (0 = up) | JUCE 8 | JUCE 9 |
|---|---|---|
| 315 degrees, upper left | 99, 114, 81 | 129, 148, 107 |
| 0 degrees, top | 49, 68, 27 | **109, 137, 77** |
| 45, 135, 180, 225 degrees | — | unchanged |

The old parser was compressing the lit end of the `skirt` linear gradient into a small patch;
lunasvg spreads it across the whole upper-left half. That is what the markup asks for and what
this file's own one-light-source rule says it should do, so **the artwork is left alone and
`docs/panel.png` is re-photographed** rather than the gradient being retuned to reproduce a
parser bug that no longer exists. The alternative was considered and rejected: matching the old
render would have meant carrying compensation numbers that nobody could later re-derive.

The two parser traps recorded at the top of `LimeGraphics.h` are kept and re-labelled as
history. The radial-radius-must-be-a-percentage one was a JUCE parser bug and lunasvg does not
have it; the never-use-`drawWithin` one was not re-tested, because scaling by the declared
viewBox is correct whatever a parser believes the drawn bounds to be.

Nothing under `Source/dsp` was touched, and its 583 checks are unchanged — which is the point of
keeping it JUCE-free: a major version bump in the framework cannot reach the process.

### State

**583 checks passing** across the six suites, 0 failures:

| Suite | Checks |
|---|---|
| `DspTests` | 280 |
| `CalibrationTests` | 75 |
| `PhaseTests` | 75 |
| `TapeTests` | 56 |
| `BarkTests` | 52 |
| `AtypeTests` | 45 |

`DspTests` is four checks up on where this version started. The halves stay power-complementary
to **−319 dB** at every count.

Standalone building clean in Release **against JUCE 9.0.0**, and the panel photographed from it
— both at the default and, to check the section ticks and the staircase actually draw, at a
temporary default of 8.

**The plugin copy step is off, and `Lime - All` builds again.** For most of this version it
did not:

```
Running rm -rf "/Library/Audio/Plug-Ins/VST3/Lime.vst3"
rm: /Library/Audio/Plug-Ins/VST3/Lime.vst3: Permission denied
Command PhaseScriptExecution failed with a nonzero exit code
```

The Install Target phase added when the build started installing system-wide writes to
`/Library`, which needs privileges the build does not have. The compile was always clean —
`libLime.a` linked — so this was a packaging step taking down the targets either side of it.
`enablePluginBinaryCopyStep="0"` removes both phases; the plug-ins are left in
`Builds/MacOSX/build/Release` and installed by hand. A copy that fails should not be able to
fail a build, and the smaller mechanism is the one that cannot.

**The `auval` trap is worth keeping on the record even though the build is fixed.** With
nothing installing automatically, `auval` reports on whatever component was last copied into
the search path, which need not be the one just built — and during this work it once returned
a confident pass for the *previous* version's binary. That claim was withdrawn rather than
quietly corrected, because the failure mode generalises: a validation that does not name the
binary it validated is not evidence. The check is the binary's timestamp and size before
believing the result.

`PhaseTests` reads 75 here where the README had recorded 51. That file is untouched by this
version's work; the old figure was stale and is corrected rather than caused.

---

## v1.0.0 — the release

### The sliding layer's opposition reformulated, then smoothed

The modulation-control opposition inside the sliding layer was rewritten to the paper's own
arrangement — §3.3's **weighted-RMS normalisation with ratio opposition** — in place of the
approximation it had carried since the layer was built. The raw weighted-RMS opposition stepped
at the frame rate, and that stepping was itself a modulation path into the gain; smoothing the
opposition with the stage's own two poles — the same poles the control signal already passes
through — removed the pumping and dropped the null floor by about 13 dB. This is the change
with a system-level measurement attached: the broadband encode/decode null through separate
instances at −26 dB improved from the subtraction era's **−46.6 dB to −67.2 dB** at the default,
once the ratio form landed and its opposition was smoothed by the same poles as the control, and
every section-count figure quoted in the resolution-control sections above was re-measured on
the new plateau. The full set now reads: steady tone at −20 dB **−80.2 dB**, broadband at
−26 dB **−67.2 dB**, broadband at −60 dB **−60.6 dB**, broadband near reference **−50.9 dB**;
loop-mode broadband nulls **−49.7 / −52.6 / −55.7 dB** at 44.1 / 48 / 96 kHz. The layer is not
inert at nonzero modulation control and MC = 0 is bit-identical to opposition disabled, both
asserted by test.

### The antisaturation shelves normalised at DC

The shelf sections of the fixed networks were normalised at DC rather than at their
asymptotes, which had left the low-frequency antisaturation **0.53 dB off at DC**. With the
normalisation corrected, the analog prototypes reproduce the paper's §4.1 figures — 10, 2, 6,
10 dB at 25 Hz, 5 kHz, 10 kHz, 15 kHz — to **0.035, 0.274, 0.041 and 0.055 dB**: within
0.28 dB everywhere, and within 0.06 dB at three of the four published points.

### The sections control: live in loop, transitioned across a pair

**Superseded during release preparation and kept for the record.** The mode-split routing
below shipped briefly and was then retired with the control itself — see *The sections
control retired* further down.

Recorded in full under *Moving the control while audio runs* above: a bare live move measured
clean (worst hop-to-hop step 1.7 dB, inside the slew clamp), and the count reshapes the encode
characteristic an encode/decode pair must agree on. The control therefore briefly split by
mode: in loop mode the encoder and decoder share the instance and move together, so it applied
live, slew-limited; in encode and decode modes a paired instance cannot follow a mid-stream
change, so there it routed through the same 1500 ms fade-to-dry transition as Process and Mode
rather than stepping mid-signal.

### The sections control retired — the resolution is a compiled-in constant

The control did not survive release preparation. Re-measured on the smoothed system, the null
against the section count improves to a plateau that begins at 256 sections and tops out at
**−67.5 dB** — and from the 256-section default upward every position measured the same and
could not be told apart by ear. The positions that could be told apart were the coarse ones,
and they were told apart by being worse: the coarsest surviving setting cost about 37.5 dB of
null against the plateau, which is audible. A control whose entire audible range is ways to
make the round trip worse offers nothing worth a knob, so the knob was removed — from the
panel and from the parameter set — and the count compiled in at **2048 sections**: on the
plateau, at the best figure the range produced, and never past one section per bin at any
rate the plugin runs at.

The retirement also took the mode-split routing recorded above, which shipped only briefly:
with no control to move there is no move to route, and an encode and its decode now agree on
the count by construction, across instances, sessions and rates — the one agreement the
control demanded of the user is no longer possible to get wrong.

The DSP layer is untouched: `SrSideChain` keeps the full `setSections` interface, and every
section-resolution and section-transition test still runs against it. Only the plugin wrapper
pins the value. The plot's faint section-boundary ticks — drawn only at 64 sections and
below — went with the control; at 2048 there is nothing to draw.

### Phase and delay left the panel

The phase and delay knobs were removed, parameters included. Both were utility controls that
moved the processed path against the dry one ahead of the blend, and both proved peripheral
to the dynamics identity the release presents: choosing a blend's phase relationship and its
alignment is another product's job, and on this panel they were two more things to explain
away from the process itself. `PhaseRotator` and `FractionalDelay64` remain in the library
with their test suites — the 75 checks of `PhaseTests` and the fractional-delay checks in
`DspTests` — no longer wired into the plugin.

### A drive knob — built, measured twice, and cut

A one-knob tape-style colour was tried in the freed panel space and removed before
release; the two laws it went through are worth the record. The first form was
slope-normalised, y = (K/g)·tanh(g·x/K) with a +6 dB knee: level-matched by construction,
and inaudible for exactly that reason — programme material at −18 dBFS met the curve
hundredths of a dB from the identity, so the knob did nothing a listener could find. The
second was the peak-normalised curve y = tanh(k·x)/tanh(k), k = 0..4 across the knob —
full scale pinned, up to 12 dB of small-signal lift, with the negative half its exact
inverse atanh(x·tanh(k))/k, and the pair asserted to cancel below −280 dB. That one was
audible, and it was cut anyway: colour is not this unit's job, and a saturator sitting on
the processed path of a system whose whole claim is exact complementarity invited exactly
the kind of misuse the bit-exact-at-zero rule existed to contain. The panel settled at
seven knobs, four over three with the lower row centred — an inverted trapezoid — and
character stays the business of the modelled channel in both mode, where it belongs.

### Auto gain

An **AUTO** toggle (`autoGain`) joined POLARITY and BYPASS: loudness-matching makeup on the
wet path before the blend. Both sides' powers run through a 3 s integrator; the correction is
sqrt(dry/wet), clamped to ±12 dB, and its application is smoothed over 400 ms. The correction
holds through silence below −70 dBFS rather than drifting after a level that is not there,
and returns to unity when the switch is disengaged. It measures against the delayed dry copy,
so bypass and mix comparisons are level-matched — an A/B between two loudnesses is a loudness
test, and this is the switch that removes that variable.

### The trims renamed for function

Send In / Send Out / Return In / Return Out became **Up In / Up Out / Down In / Down Out** —
host-visible names and panel captions ("up in" and so on, in the captions' lower case), with
the parameter identifiers unchanged so saved sessions load, exactly as when the same trims
left the hardware's record/play vocabulary. The new names are the modes' own vocabulary, and
they carry the functional reading the transmission framing hid: the in-trims set how hard the
programme drives each pass — how much of it falls into the action region, which makes them
this architecture's threshold control — and the out-trims are that pass's level compensation.

### The panel is a four-by-four grid, and the plot follows the direction

With exp, phase and delay gone and drive arrived, the knobs settled into two rows of four:
**up in, up out, down in, down out** across the top, **input, drive, mix, output** across the
bottom, with the process and mode switches beside them and three toggles — POLARITY, BYPASS,
AUTO — on the third line. The footer reads "Aptitude Audio | Lime | v1.0.0".

The plot's gain axis became symmetric — ±28 dB, labelled every 8 dB from −24 to +24, around
an emphasised 0 dB centre line — and its legend reads **gain** rather than **boost**, because
in down mode the plot now shows the applied cut below the line. That fixed a real bug: only
the encode side chains published frames, so in down mode the plot froze on its last frame.
The probes follow the active direction now, asserted by test.

### Corrections, each with its test

- **Fractional delay floors its read position at one sample.** An engaged delay below one
  sample would ask the interpolator to read history it does not have; the floor holds
  amplitude flat to 0.01 dB and leaves a 0.5-sample request's energy entirely inside the
  interpolator's eight taps, below measurement. Exactly zero remains a bit-exact bypass.
- **The delay line repositions instead of flushing.** A `setDelay` keeps every sample the new
  geometry can still address — shrinking loses nothing, growing gaps only the samples that
  never existed — where a flush put silence into a running path.
- **Host bypass drives the plugin's own parameter.** A host-initiated bypass now lands on the
  same delay-matched, crossfaded bypass the panel button uses, instead of the host's abrupt
  substitution.
- **Oversized host blocks are sub-chunked to the promised size**, so `process()` never
  allocates after `prepare()` — asserted with an allocation trap across oversize blocks and a
  live section move, in both gain realisations.
- **Auto Compare runs on a per-sample clock**, so its four-second alternation no longer
  quantises to the host block length.
- **State restore snaps.** A restored session's settings are in force from the first sample —
  the smoothers are seeded rather than ramped — and the state carries a version stamp so a
  later format can tell what it is reading.
- **The display probes live for the engine's lifetime** and are handed to the editor through a
  wait-free triple buffer; 44898 reads raced against 10000 pushed frames found zero torn
  snapshots. `SpectrumProbe` is header-only and owns no lock.
- **The main-path fit is memoised per rate**, dropping prepare-time fitting from about 416 ms
  to under 3 ms, which is the difference between a rate change stalling the host and not.

### How the release presents the product

The release presents Lime as a **spectral companding dynamics processor** first — per-bin
compression one way, its exact complement the other, usable single-ended as an insert: upward
compression in encode, the matching expansion in decode — with the paired use around a lossy
channel stated as the published noise-reduction duty of the architecture rather than the
product's only purpose. The README and the manual lead with the dynamics identity; this record
keeps its historical framing, since it is a record of what was built and why.

### Build

The project builds with **CMake** alongside the Projucer exporter, and CI runs the six suites
on every push — on macOS through the Accelerate backend, on Linux and Windows through the
portable paths the suites cross-validate against.

### State

**652 checks passing** across the six suites, 0 failures:

| Suite | Checks |
|---|---|
| `DspTests` | 343 |
| `CalibrationTests` | 81 |
| `PhaseTests` | 75 |
| `TapeTests` | 56 |
| `BarkTests` | 52 |
| `AtypeTests` | 45 |

---

## v1.0.1 — the side-chain filter pair

The feature was built twice. The first form was a **low-pass on the wet path** — third-order
Butterworth in direct-form-I sections because the corner sweeps, prewarped to put the −3 dB
point on the knob's number — shading the process's harmonics out of a parallel blend. It
worked as specified and was replaced at the user's direction by the better placement: filter
what the detection hears, not what the listener does. `WetLowPass` stays in the library with
its tests, wired to nothing — the precedent the drive knob set.

The shipped form is a pair of knobs — `scHpf` and `scLpf`, "Sidechain HPF" and "Sidechain
LPF" to the host, **s/c hpf** and **s/c lpf** on the panel — a **high pass and a low pass
on the compressors' level detection**: third-order Butterworth *magnitudes* —
`butterworth3Magnitude` in `DspMath` — applied as weights rather than run as filters, so
there is no bilinear warp and each corner is −3.01 dB exactly. The two weights multiply —
their decibels add — so together the pair band-limits what the compressors hear. The
spectral side chains weight per bin, `SrSideChain::setDetectorHighPass` and
`setDetectorLowPass` fanned out by `SrEngine::setSidechainHighPass` and
`setSidechainLowPass` exactly as `setSections` is — the weight lands on the `running`
levels and leaves `pending` untouched, so shared state stays bit-identical — and the
a-type engine weights each of its four band detectors at the band's geometric-centre
frequency. Outside the corners the programme reads quieter than it is and the process
treats it as low-level material; the plot follows, because the gains it draws are the
gains that moved.

Both filters rest **off** — `scHpfOn` and `scLpfOn`, false by default. The decision is
resting neutrality: a fresh instance detects with the full band, and the knobs' corners —
80 Hz and 6 kHz, each at noon — are positions held in readiness rather than filters
already in circuit, so a filter enters the detection path by deliberate engagement and
never by default. Off, the engines receive no corner at all — unity weighting, identical
to the pre-filter behaviour — and a throw is paced by the detection's own time constants
and the transfer slew limit, exactly as a corner move is.

The property the wet-path form could not have: **any corners leave the round trip exact**,
because encode and decode weight their detection identically and so compute the same gains
to apply and to invert. Measured, the loop-mode broadband null is −51.6 dB with the
high-pass corner at 6 kHz — the same figure as without it — and −57.2 dB with both corners
engaged at 200 Hz and 6 kHz; the a-type round trips null at −272.6 and −279.8 dB. The
filters join the section count as something a pair must *agree* on — switch state and
corner alike — rather than something that breaks the pair. `testSidechainHighPass` and `testSidechainLowPass` in `DspTests`
cover the prototype magnitudes, the loop nulls with the corners engaged, and the
detection-weighting behaviour of both engines.

The pair fills the bottom row the drive knob's departure left short: **s/c hpf, s/c lpf,
input, output**, four wide again, the half-pitch centring inset gone; the switch column
grew a fourth line — captioned **s/c filters**, the HPF and LPF switches side by side —
keeping every switch in the one column, now on four shared baselines, the last level with
the knob row it switches. The second position
was mix's, and **mix left the panel** to make the room — the parameter remains,
host-visible and automatable, default 100%, the blend arithmetic and the bit-exact-dry
zero unchanged. The precedent is Set Up's: a control whose real-two-instance setting is
its default, kept off the panel deliberately with the capability intact, reachable from
the parameter list by the sessions that want it.

---

## v1.0.2 — the per-half time constants

The feature was designed twice before it was built once. The first design was a single
linked pair — one attack and one release governing both halves at once — and it became
per-half independence at the user's direction: four knobs, `upAttack`, `upRelease`,
`downAttack` and `downRelease` — "Up Attack" and so on to the host — so the encode and
decode halves each carry their own global ballistics. Attack runs 0.1 to 1000 ms and
release 10 ms to 5 s, each skewed to put its default — 5 ms and 100 ms — at noon, the
side-chain corners' precedent; readouts change units at a second ("5.0 ms", "1.20 s"),
and typed entry takes a bare number as milliseconds and "1.2 s" or "2s" as seconds.

What the times override is deliberately only the **final** smoother of each spectral
stage's detection — the fixed band's 160 ms (HF) / 300 ms (LF) pole after the maximum
selector, and the sliding band's 80/150 ms finals. The staggered fast trackers ahead of
them — the 15 ms fixed-band mains and pass-band, the 5/7.5 ms sliding firsts — keep their
published values, because the stagger *is* the architecture: the character of the process
comes from fast trackers at deliberately different speeds feeding a slower pole, and a
knob that flattened that stagger would be a different process wearing the same panel.
What the user shapes is the smoother that dominates how the detection settles, and
nothing else. Engaged, a final runs asymmetric — the attack time while its input is
rising, the release time while it falls — so one smoother carries both knobs of a pair.

The sliding band's modulation-control opposition is smoothed by the same finals, and it
deliberately **follows the user pair**. The code's own invariant is that the opposition
must be as smooth as what it opposes; letting the per-bin control run user times while
the MC reference kept the published pole would let the two sides of the end-stop ratio
settle at different speeds, transiently mis-forming the ratio on every level change. So
the reference moves with the control, and the invariant holds at any setting.

Engagement is a per-value **zero sentinel**: zero for either value keeps that direction's
published tau, and with a switch off both values are zero, where the sentinel reduces to
the identical arithmetic — measured bit-identical in both engines, asserted at tolerance
zero, so `Up Times On` and `Down Times On` resting off is the published unit exactly. The
knobs hold times in readiness, the arrangement the side-chain filter switches
established. Each pair is also provably confined to its own direction — encoding under
down-pair changes and decoding under up-pair changes both measure a difference of
exactly 0.

The a-type engine maps the pairs onto its own ballistics: the up pair replaces the
steady-state attack rate — published 0.1 s, proportional to the transition — and the
33 ms recovery while encoding, the down pair while decoding. The peak/average detector
constants and the 1 ms fast-attack cap stay published, so a user attack below 1 ms is
bounded by the cap — the knob's bottom decade is a limit honoured, not a control ignored.

The headline property is a trade, stated plainly: **the round trip is exact only while
the up and down halves agree** — same switch states, same times. Measured, matched pairs
engaged at 5 ms/100 ms null at −41.3 dB in loop mode (limit −26) and −277.3 dB through
the a-type round trip (limit −100); halves deliberately mismatched — up at 5/100, down at
50/1000 — null at −19.1 dB, asserted as a property so a later change that "improves" it
gets noticed. Splitting the pairs is a legitimate creative choice that forfeits the null;
matching them, or leaving both off, preserves it. The agreement discipline the section
count and the s/c corners established gains an intra-instance clause — and the times
carry a caveat the s/c weights did not: because they change the smoothers' stored state
trajectories rather than a memoryless weighting, two instances that change a time at
different moments drift until the poles settle, so a *moving* time knob across an
encode/decode pair does not null until it lands.

The panel grew a third knob row between the pass trims and the masters, each half's pair
sitting in the columns of its own trims — levels above, ballistics below, per half — and
the masters row is now **input** and **output**, centred on the four-column grid; the
switch column's fourth line is recaptioned **time constants**, an UP and a DOWN switch
level with the knob row they govern, and the panel is one knob-row taller.

`testGlobalTimeConstants` in `DspTests` measures the lot: the matched and mismatched
nulls above, the exact-0 disengaged and cross-direction differences, and the audible
asymmetries — a 20 ms up release restores boost 14.7 dB sooner than a 2 s one within the
measurement window (7.4 dB in the a-type), and a 1 s up attack sheds boost 8.4 dB later
than a 1 ms one.

### State

**698 checks passing** across the six suites, 0 failures:

| Suite | Checks |
|---|---|
| `DspTests` | 389 |
| `CalibrationTests` | 81 |
| `PhaseTests` | 75 |
| `TapeTests` | 56 |
| `BarkTests` | 52 |
| `AtypeTests` | 45 |

---

## v1.0.3 — always-on times, and a smaller window

The `Up Times On` and `Down Times On` switches shipped in v1.0.2 and were removed within
a day, and the reason deserves recording plainly because it is the only kind of reason
that decides such things: the panel's own author looked at the UP and DOWN buttons on the
**time constants** line and could not tell what they did. Not a hypothetical user, not a
usability review — the person who designed the line, the day after shipping it. A control
its own author cannot read has failed the only measurement of clarity that matters, and
the failure was in the idea rather than the caption: two direction-named switches
engaging four knobs on a different row is a relationship a caption can assert but a
glance cannot recover. So the fix was not a better caption. The times were made
**always-on**, which removes the ambiguity at the root — a knob's position is now always
the constant in force, never a time held in readiness behind a switch elsewhere — and the
two boolean parameters left the plugin, which now publishes 21 host parameters.

What always-on changes about the stock sound, and what it does not. A fresh instance now
runs 5 ms attack / 100 ms release on both halves — the knobs' defaults *are* the unit's
ballistics — where v1.0.2's resting state was the published constants exactly. Those
published finals (160/300 ms fixed band, 80/150 ms sliding) are no longer reachable from
the plugin: they remain a capability of the DSP library, whose zero sentinel still
reduces to the identical arithmetic and is still asserted bit-identical by
`testGlobalTimeConstants` — the published unit survives at the library level; the plugin
simply never sends the zero any more. The round-trip discipline is unharmed, because both
halves default identically: a fresh pair agrees out of the box, and matched times at the
defaults measure −41.3 dB in loop mode and −277.3 dB through the a-type round trip — the
v1.0.2 figures, standing, along with the mismatched −19.1 dB and the audible
asymmetries.

The one cost lands on existing sessions and is stated rather than hidden: a v1.0.2
session loads, but a session that had a switch off now runs its stored knob times — the
switch it relied on no longer exists — so its detection ballistics change from the
published constants to whatever its knobs held in readiness. A session that never touched
the time knobs moves from the published finals to 5/100, which is the same change a fresh
instance embodies; a session that parked exotic times behind an off switch gets those
times in force.

The panel follows the parameters: the switch column is back to three rows — process,
mode, polarity/auto/bypass — while the knob grid keeps its three, the pass trims over
their time constants over the centred masters.

The window now draws at **75% of its design size**. The third knob row had made the
full-size panel taller than smaller screens comfortably give it, and the
resizable-corner note's argument still holds — scaling properly means every metric in
LimeStyle becoming a ratio and the artwork re-rasterising per size, real work for a
panel whose one job is to be read. So the panel does not do that. One
`AffineTransform::scale (windowScale)` — `windowScale = 0.75f` in `LimeStyle.h` — on a
content component that holds the entire panel presents the whole window smaller while
the layout keeps thinking in design pixels: every metric scales together through the
one transform, so nothing can fall out of step with anything else, which is precisely
the failure mode the old resize handle had when it moved the frame and left the metrics
behind. The transform lives on that child rather than on the editor because the
editor's own transform belongs to the host — JUCE asserts on any other owner, since a
host rescaling for DPI would silently obliterate it. One consequence of three switch
rows dividing three knob rows' height: the switches would have stretched into slabs, so
each now caps its height at its own width and renders square at any column height.

A pre-release performance sweep closed the version out, and its one real find predated
v1.0.0. `shapeRatio` special-cases the squared knee precisely so no `std::pow` runs — but
inside `FixedBand::processFrame`'s auto-vectorised loop the compiler speculated all three
of its branches, executing the general `pow` unconditionally per bin and blending the
result away: about two million discarded libm calls a second, invariant across the
section count, and verified in the emitted assembly. The frame loop is now parameterised
on the knee shape and chooses once per frame, so the squared knee's loop contains no
`pow` at all; the arithmetic is identical bit for bit, which the suite confirms — every
measured figure, the nulls included, is unchanged to the digit. Measured cost on Apple
silicon fell about 13%: the full loop round trip runs at 14.2% of one core stereo at
48 kHz, a single encode pass at 7.0%, and the a-type at 0.1%. The same sweep found the
scaled content component defeating JUCE's opaque-child culling — the panel gradients
were repainting under every meter tick and plot frame, because the culling walk skips
transformed children — so the content component now paints the panel itself, opaquely,
one level down where its children are untransformed and the dirty rects clip it to
nothing. The engine-level corner and time setters gained the unchanged-value guards
their chain-level halves already had (with prepare re-applying through unguarded apply
functions, so freshly built channels are not starved by their own equality checks), and
the allocation test arms all three setters alongside the section move it already
exercised.

The time ranges narrowed after listening rather than measuring. The attack knobs
originally ran from 0.1 ms to a second, and their bottom decade did nothing audible: the
published fast trackers ahead of the governed final smoother — the sliding band's 5 and
7.5 ms firsts, the fixed band's 15 ms mains — bound the detection's response whatever
the final is set to, so 0.1 ms and 5 ms differed as numbers and not as sound. The
principle that settled it: a range must not offer what cannot be heard. The final spans,
chosen by ear, are attack 1 to 100 ms and release 10 ms to 2.5 s — tight enough that
both ends of each knob are audible on dynamic material, the attack floor landing exactly
on the a-type's 1 ms fast-attack cap so the cap caveat disappeared with the dead zone,
and both stock defaults back at noon where this panel likes them.

### State

Check counts are unchanged — **698 checks passing** across the six suites, `DspTests` at
389 — and `auval` passes at 1.0.3.

---

## v1.0.4 — the words

A documentation release; no code changed. Its one finding is worth its own record: asked
whether "upward compressor" was the wrong name for the up pass, the answer came from the
engine rather than from taxonomy debate. Measured on a steady 1 kHz tone, the up pass's
gain is always positive — +20.5 dB at −80 dB input, +1.1 dB at −10 — and a 70 dB input
range leaves as 50.6 dB: gain that rises, range that shrinks, which is upward
compression by definition. The down pass mirrors it exactly — gain always negative,
70 dB in becomes 89.6 dB out — expansion whose gain moves downward. The nouns name what
happens to the range; the modifiers name where the gain goes. So the up pass kept its
name, and the down pass gained the direction word it had been missing — **downward
expander** — everywhere it is described. The remaining changes reconciled the manual and
README with the shipped state: the one-command `make all` flow, the attack figure
attributed to the library's test span rather than to unreachable knob positions, the
tools preflight claimed only for the targets that run it, and the `WetLowPass` row in
the architecture table.

### State

**698 checks passing** across the six suites, `DspTests` at 389; `auval` passes at
1.0.4. This is the final release.
