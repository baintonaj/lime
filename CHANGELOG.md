# Changelog

## 1.0.2

- Four new knobs — **up attack**, **up release**, **down attack** and **down release** —
  give the two halves of the pair independent global attack and release time constants:
  one pair for the up (encode) pass, one for the down (decode) pass. Attack runs 0.1 to
  1000 ms and release 10 ms to 5 s, each skewed so its default — 5 ms and 100 ms — points
  straight up, the side-chain corners' precedent; readouts switch to seconds above one
  ("5.0 ms", "1.20 s"), and typed entry takes a bare number as milliseconds and "1.2 s" or
  "2s" as seconds. What an engaged pair overrides is deliberately only the **final**
  smoother of each spectral stage's detection — the fixed band's 160 ms (HF) / 300 ms (LF)
  pole after the maximum selector, the sliding band's 80 and 150 ms finals — run
  asymmetric: the attack time while the level rises, the release time while it falls. The
  fast published trackers ahead of them (the 15 ms fixed-band mains and pass-band, the 5
  and 7.5 ms sliding firsts) keep their published values, so the staggered character
  survives; in the a-type process the up pair replaces the steady-state attack rate and
  the 33 ms recovery while encoding and the down pair while decoding, with the detector
  constants and the 1 ms fast-attack cap untouched.
- Each half has its own switch — `upTimeOn` and `downTimeOn`, both resting **off** — on
  a switch line now captioned **time constants**. Off, the ballistics are the published
  ones exactly — measured bit-identical in both engines, a difference of exactly 0 — and
  the knobs' positions are times held in readiness rather than constants in force, the
  arrangement the side-chain filter switches introduced. Each pair provably never reaches
  the other direction's detection, also at exactly 0.
- **The round trip is exact only while the up and down halves agree** — same switch
  states, same times. Measured, matched pairs engaged at 5 ms/100 ms null at −41.3 dB in
  loop mode and −277.3 dB through the a-type round trip; halves deliberately mismatched
  (up at 5/100, down at 50/1000) null at −19.1 dB — the price of independence, asserted as
  a property rather than treated as a bug. Splitting the pairs is a legitimate creative
  choice that forfeits the null; matching them, or leaving both off, preserves it. The
  agreement discipline the section count and the s/c corners set gains an intra-instance
  clause — and because the times change the smoothers' stored state, a time knob *moving*
  across a two-instance pair does not null until it lands.
- The **s/c hpf** and **s/c lpf** knobs and their switches left the panel to make the
  room. All four parameters remain — Sidechain HPF at 80 Hz, Sidechain LPF at 6.00 kHz,
  both switches resting off — host-visible and automatable, the arrangement mix and Set Up
  have; the filters themselves are unchanged.
- The knob grid is now three rows: the pass trims, a new row of their time constants —
  each half's pair in the columns of its own trims — and a masters row of **input** and
  **output**, centred on the four-column grid. The switch column's fourth line carries the
  UP and DOWN switches, level with the knob row they govern, and the panel is one knob-row
  taller.
- `DspTests` grew `testGlobalTimeConstants` — the nulls above, the bit-identity of the
  disengaged state, and the audible asymmetries: a 20 ms up release restores boost 14.7 dB
  sooner than a 2 s one within the measurement window (7.4 dB in the a-type), and a 1 s up
  attack sheds boost 8.4 dB later than a 1 ms one. `DspTests` is now 389 checks and the
  suites total 698.

## 1.0.1

- Two new knobs — **s/c hpf** and **s/c lpf** — the classic side-chain filters, realised
  here as a high pass and a low pass on the compressors' level detection rather than
  filters in any signal path. Each is a third-order Butterworth magnitude (18 dB per
  octave, −3.01 dB exactly at its corner) applied as a weight: per bin to the detection
  levels of all four spectral side chains, and per band to the four a-type detectors, each
  at its band's geometric-centre frequency. Outside the corners the programme reads
  quieter than it is, so the process treats it as low-level material; the two weights
  multiply — their decibels add — so together the pair band-limits what the compressors
  hear; the plot reacts because the computed gains themselves move. Both corners run
  20 Hz to 20 kHz on skewed knobs whose defaults point straight up — 80 Hz for the high
  pass, the classic corner, and 6 kHz for the low pass.
- Each filter has its own switch — `scHpfOn` and `scLpfOn`, both resting **off** — so a
  fresh instance detects with the full band and the knobs' corners are positions held in
  readiness rather than filters already in circuit. With a switch off the engines receive
  no corner at all — unity weighting, identical to the pre-filter behaviour — and a throw
  is paced by the detection's own time constants and the transfer slew limit, like a
  corner move.
- Because encode and decode weight their detection identically, **any corner settings
  leave the encode/decode round trip exact** — measured, the loop-mode broadband null is
  −51.6 dB with the high-pass corner at 6 kHz — the same figure as without it — and
  −57.2 dB with both corners engaged at 200 Hz and 6 kHz; the a-type round trips null at
  −272.6 and −279.8 dB. The discipline is the section count's: an encode and its decode
  must agree on the filters — switch states and corners both — so across two instances
  they are set the same at both ends.
- The bottom row is four knobs again — **s/c hpf**, **s/c lpf**, **input**, **output** —
  restoring the four-wide master row; the half-pitch inset that centred the three-knob row
  is gone. The switch column grew a fourth line for the pair — captioned **s/c filters**,
  the HPF and LPF switches side by side, level with the knob row it switches. Mix gave up
  its knob to the low pass and left the panel; the parameter remains, host-visible and
  automatable, default 100%, with the blend arithmetic and the bit-exact-dry zero
  unchanged — the arrangement Set Up has always had, deliberately off the panel with the
  capability intact.
- A wet-path low-pass was built and measured first and replaced by the side-chain
  placement; `WetLowPass.h` stays in the DSP library with its tests, wired to nothing. The
  weighting adds a `butterworth3Magnitude` helper to `DspMath`, and `DspTests` grew
  `testSidechainHighPass` and `testSidechainLowPass` — prototype magnitudes, the loop
  nulls with the corners engaged, and the detection-weighting behaviour of both engines.

## 1.0.0

- A CMake build and continuous integration were added, so the plug-ins and every test
  suite build and run from one command on a clean checkout. The Projucer project remains
  as a secondary route.
- The host's own bypass now drives the plugin's bypass parameter and runs the same
  delay-matched crossfade as the panel switch, so bypassed audio honours the reported
  latency instead of jumping 213 ms ahead of the session.
- Real-time-safety hardening: the display probes are engine-lifetime and hand frames to
  the editor through a wait-free triple buffer; oversized host blocks are processed in
  promised-size pieces instead of growing buffers on the audio thread; the remaining
  editor-visible state is atomic. An allocation-counting engine test now stands guard.
- The sliding layer's modulation-control opposition was reformulated as the ratio the
  source paper describes: MC signals are weighted RMS levels on the per-bin magnitude
  scale, smoothed by the same two poles as the control they oppose, and the level raises
  the compression-law reference. The old per-bin subtraction of a whole-spectrum sum
  silenced the layer entirely on broadband material and varied with the transform size.
  The measured round-trip nulls improved with it: steady tone −68 → −80 dB, broadband at
  −26 dB −55 → −67 dB.
- The fractional delay no longer reads a wrapped tap at sub-sample settings — engaged
  delays carry a documented one-sample floor — and the shelf prototypes now pin their
  digital DC gain to the analog value exactly in both conversion branches, which moved
  the measured curves *closer* to the published figures (antisaturation now within
  0.06 dB at three of the four published 10/2/6/10 points and 0.28 dB at the fourth).
- Auto Compare advances its segment clock once per sample rather than once per channel,
  so a stereo pair switches source together and segments keep their length. The section
  resolution control briefly split its routing by mode during release preparation — live
  in loop mode, where encoder and decoder move together; fade-to-dry across a pair — and
  was then retired outright with the panel consolidation below. Restored sessions snap
  into force instead of fading in over the previous settings, and the saved state carries
  a version stamp.
- Mode switches reposition the padding delay instead of flushing it, ending a
  213 ms dropout on every switch.
- Performance: the main-path filter fit dropped from ~416 ms to under 3 ms per rate (and
  is memoised); the per-bin dB conversions are vectorised; settled controls no longer pay
  per-sample transcendentals; the decode path lost three redundant biquad copies per
  sample.
- The exp control is gone from the panel and the parameter set, and the band split's
  resolution is a compiled-in constant of 2048 sections — on the measured null plateau, at
  its best figure (−67.5 dB), and never past one section per bin at any rate. From the old
  default upward every setting measured and sounded the same; the settings that differed
  were the coarse ones, and they differed by being worse, and a control offering only ways
  to make the round trip worse was not worth its knob. The DSP layer keeps the full
  `setSections` interface and its tests; only the wrapper pins it.
- The phase and delay knobs are gone — two utility controls that moved the processed path
  against the dry one, peripheral to a dynamics unit. `PhaseRotator` and
  `FractionalDelay64` remain in the library, still tested, no longer wired into the
  plug-in.
- An auto-gain switch: loudness-matching makeup on the wet path before the blend — both
  paths' powers through a 3 s integrator, correction √(dry/wet) clamped to ±12 dB and
  applied through a 400 ms smoother, holding through silence below −70 dBFS. It measures
  against the delayed dry copy, so bypass and mix comparisons are level-matched.
- The four trims are renamed for what they do: Send In / Send Out / Return In / Return Out
  became Up In / Up Out / Down In / Down Out — the in-trims set how hard the programme
  drives each pass, which is this architecture's threshold control, and the out-trims are
  that pass's level compensation. Parameter identifiers are unchanged, so saved sessions
  load.
- The panel settled into seven knobs, four over three with the lower row centred — the
  pass trims above; input, mix and output beneath — with polarity, auto and bypass beside
  them. A tape-drive knob was built, measured through two laws, and cut before release:
  colour is not this unit's job. The plot's gain axis is now
  symmetric at ±28 dB about an emphasised 0 dB centre, labelled every 8 dB, its legend
  reads "gain" rather than "boost", and in down mode it shows the applied cut below the
  line — fixing a freeze where only the encode side chains published frames.
- Documentation corrected and expanded — the latency table now matches the code — and the
  suites grew to 652 checks, including regressions for every fix above.

## 0.2.0

- The crossover control became a resolution control. The split between the two staggered
  halves is a compiled-in constant again at the reference design's 800 Hz; what the user
  sets is how finely that split is resolved, onto 4 to 16384 log-spaced sections in powers
  of two. The knob shows the exponent, 2 to 14, defaulting to 8 — 256 sections. Measured,
  the encode/decode null improves as the grid gets finer and plateaus at about 64 sections;
  the coarsest setting costs about 17 dB of null, so it is a different round trip rather
  than a cheaper one.
- The section grid is defined in hertz rather than bins and is never clamped against a
  half's own bin count, keeping the two halves power-complementary — measured to −319 dB at
  every count. Moving the control while audio runs is covered by the existing gain slew and
  is verified by test.
- Moved from JUCE 8.0.14 to JUCE 9.0.0. SVG parsing is now lunasvg; the one affected call
  site was updated, the panel was photographed before and after and differenced, and the
  knob-skirt gradient now renders as the artwork asks rather than as the old parser's bug.
  Nothing under `Source/dsp` changed.
- The build no longer installs plug-ins as a build phase, so the umbrella target builds
  again; artefacts are left in the build folder and installed by hand. The `.jucer` was
  normalised back to LF and the shared Xcode scheme is tracked.
- The modelled channel now stays in circuit in loop mode whatever the process is set to,
  so the off-versus-b comparison the panel is built around is audible: −85.0 dB of hiss
  untreated against −98.5 dB treated, on a 220 Hz burst at 48 kHz.
- The plot no longer reads "no signal" against a live engine: probe pointers are
  re-fetched every timer tick instead of captured once at construction.
- The window is a fixed size; the resize handle that silently dropped controls is gone.
- Plot captions are plain ASCII, avoiding a string-constructor encoding trap that turned
  an em dash into mojibake.
- Process and mode transitions lengthened to 1500 ms.
- 583 checks passing across six suites.

## 0.1.1

- CPU cut by 2.6× — 40.8% of a core to 15.5%, stereo loop mode at 44.1 kHz — with every
  measured null figure held. The scalar math-library calls in the per-bin loops were
  removed or tabled, the envelope tree is version-stamped instead of cleared every frame,
  and the low half's analysis hop was raised to 1024, ending sixteen-fold overlap where
  four-fold was carrying the information.
- New master controls: input and output at ±24 dB, mix from 0 to 200% (past 100 the
  wet/dry difference is applied again), a true ±180° phase rotation on the processed path,
  a 0–40 ms utility delay with fractional interpolation, and polarity.
- Three real bugs fixed: switching processes no longer replays stale audio from an unused
  delay line; all four calibration trims now sit at their hardware positions in loop mode;
  and the plot no longer animates under bypass or shows a stale curve when nothing is
  publishing frames.
- Check Source was removed — in decode-only mode it was bypass with extra steps. The
  crossover became continuous and reached down to 200 Hz, with the taper skewed to keep
  800 Hz straight up.
- Process and mode switches blend over a raised-cosine ramp rather than a one-pole;
  bypass keeps its fast 75 ms crossfade.
- Built plug-ins were installed to the system-wide plug-in folders (reverted in 0.2.0 when
  the copy step proved able to fail the whole build).

## 0.1.0

- First version. A two-pass companding noise reduction system for a lossy channel,
  reinterpreting the analog multiband design introduced in 1986 and described in its
  1987 paper as an STFT process in which
  every bin carries its own stack of level-dependent stages. AU and VST3, 64-bit double
  precision throughout.
- Two processes: A, four fixed bands; B, five staggered stages across two spectral halves.
- The fixed main-path networks are exactly invertible; the round-trip shortfall of the
  spectral gain is measured and documented rather than hidden.
- An adjustable crossover, 320 Hz to 2 kHz, symmetric in ratio about 800 Hz.
- Every control smoothed by a one-pole, 95% of any change inside 75 ms.
- The DSP layer is framework-free and exercised by five standalone suites, 460 checks.
- The full design record lives in the repository alongside the user manual.
