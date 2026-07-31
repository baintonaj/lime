# Lime

**Spectral companding dynamics processor — user manual** · Aptitude Audio · v1.0.1

# What Lime is

Lime is a **spectral companding dynamics processor** — per-bin compression one way, its
exact complement the other — built on a classic noise reduction architecture.

The forward pass compresses the signal's dynamic range from both ends at once: quiet detail
is lifted, and the extremes near saturation are eased. The reverse pass is its exact
complement, and expands both back out. Run one pass alone and Lime is an insert: **up** is
a many-band upward compressor, raising low-level detail against everything around it, and
**down** is the matching expander, pushing it back away.

The lift is not a fixed EQ curve and not a single compressor. The spectrum is divided into
bands; each band carries its own stack of level-dependent stages; each stage acts only where
the signal is quiet enough to need it. A loud trumpet gets almost no treatment while the room
around it gets 24 dB of lift, in the same band, in the same instant.

That principle — **treat only what needs treating** — is what the plot on the panel shows you
directly. The curve sits high wherever the signal is quiet, and dips at and around whatever
is loud.

## The architecture's origin

The design Lime follows is a published **companding noise reduction system**, and running
the two passes as a pair around a noisy channel — tape, a long line, a lossy codec, a
transmission path — is the duty it was published for. On the way **out** to the channel,
the compression lifts quiet detail up out of the channel's noise floor. On the way
**back**, the expansion puts that detail down where it belongs — and the channel's noise
goes down with it. The channel never sees the quiet parts at their true level, so it never
gets the chance to bury them. **both** mode is that round trip demonstrated inside one
instance, modelled channel included.

## One pass or two

A single pass needs no channel and no partner. Point **up** at anything whose quiet detail
should sit higher — a decaying tail, room sound, the small events a dense mix is sitting
on — and use **mix** to set how much; **down** is the same tool leaning the other way.

The paired use is for the cases where something in the path genuinely adds noise you cannot
remove afterwards:

- tape, real or emulated
- an analogue insert loop, a console, an outboard chain
- a lossy codec round trip
- radio, satellite, or any transmission link
- anything summed down to a low bit depth on purpose

If your channel is a bit-transparent 32-bit float bus, the pair has nothing to defeat and
will only give you back its own residual. The round trip is a tool for a lossy path; the
single pass does not need one.

## Two processes

| | Bands | Lift at low level | Character |
|---|---|---|---|
| **A** | four fixed bands | 10 dB, rising to 15 dB at 15 kHz | Gentle, simple, forgiving |
| **B** | five staggered stages across two spectral halves | 24 dB high, 16 dB low | Far more reduction, far more selective |

**A** is the simpler, earlier arrangement: four bands, one threshold each, a fixed amount of
lift. It is unfussy and hard to hear working.

**B** stacks five stages at three thresholds — approximately −30, −48 and −62 dB relative to
reference — and splits them across two frequency halves, three above the 800 Hz split and two
below. Because the stages are staggered rather than parallel, their effects multiply, and the
system discriminates far more sharply between what is loud and what is quiet. **B** is the
default.

# Installing

Two formats are provided, and they install here:

```
~/Library/Audio/Plug-Ins/Components/Lime.component     Audio Unit
~/Library/Audio/Plug-Ins/VST3/Lime.vst3                VST3
```

Both are the same processor with the same interface. Restart your host after installing so it
rescans. The AU passes Apple's `auval` validation.

Lime is **stereo or mono in, matching out**, and mono costs exactly half of stereo — nothing
is allocated or computed for a channel that is not there. Channels are processed independently
and never linked — the process tolerates channel-to-channel gain error well enough not to need it, and
linking would mean one channel's transients modulating the other's noise floor.

# The panel

![The Lime panel](docs/panel.png)

Reading order, top to bottom: what the process is doing to the signal, said as a curve; then
two rows of knobs with every switch in one column beside them; then identity and level.

The top row of knobs belongs to the two passes — **up in** and **up out**, then **down in**
and **down out**, each pair around its own pass. Beneath them is the plugin's own row of
four: **s/c hpf** and **s/c lpf**, which band-limit what the level detection hears, then
the channel trims **input** and **output**. Mix is no longer on the panel — the parameter
remains, reached from the host, and has its own section below.

Every switch is in one column to the right of both rows: **process** and **mode** on the top
two lines, **polarity**, **auto** and **bypass** on the third, and the **s/c filters** line
last — the HPF and LPF switches, level with the knob row they switch. The line under it all
reads *Aptitude Audio | Lime | v1.0.1*.

## process — off / a / b

Which system is in circuit.

- **off** — no processing. The signal passes through, delayed to match the reported latency,
  so switching to **off** does not shift anything in time relative to the rest of the
  session. In **both** the modelled channel is still in circuit — see below.
- **a** — the four-band process.
- **b** — the five-stage staggered process. Default.

## mode — up / down / both

- **up** — the forward pass: the upward compressor, and the encode half of the pair. On
  its own it is the insert; paired, use it **going to** the channel.
- **down** — the reverse pass: the matching expander, and the decode half. It reverses the
  lift exactly; paired, use it **coming back**.
- **both** — up, through a modelled tape channel, and back down again, all in one instance.

**both** is how you *hear the architecture's original duty*. It runs the complete round trip
internally, with a modelled channel in the middle: hiss, saturation, high-frequency loss and
modulation noise.

**In both, the modelled channel stays in circuit whatever process is set to.** That is the
point of the mode: the channel is the thing being demonstrated and the process is the thing
being toggled. So **off** is that channel untreated, and **b** is the same channel with the
noise reduction across it. Play something, stop, and listen to the gap — measured during
development on a 220 Hz burst at 48 kHz, the hiss sat at −85 dB with the process off and
−98.5 dB with **b**.

Note that **a** does not reduce noise in **both**. The four-band process encodes here but has
no decode pass of its own in this mode, so what you hear is the encoded signal through the
channel rather than a completed round trip. Use **b** for the demonstration, or two instances
for **a**.

To work across a real channel you want two instances — one **up**, one **down**.

## The trims — up in, up out, down in, down out

Four trims, ±10 dB, one at each side of each pass. They come in pairs named for the pass they
serve, and only one pair is in circuit at a time:

| Trim | In circuit when | Position |
|---|---|---|
| **up in** | mode is **up** or **both** | before the up pass |
| **up out** | mode is **up** or **both** | after the up pass |
| **down in** | mode is **down** or **both** | before the down pass |
| **down out** | mode is **down** or **both** | after the down pass |

In **both** all four are live at once, because in **both** the whole chain is inside the one
instance: the up pass, the modelled channel, and the down pass. **up out** sets what is
printed to that channel, **down in** sets what comes back off it, and the pair between them
is a level error you can introduce deliberately and then hear the decoder cope with — or
fail to. In **up** and **down** only the pair belonging to that direction does anything, since
the other half of the round trip is in the other instance.

**The in trims are this architecture's threshold control.** The whole system is
level-dependent: each pass decides how much to act by measuring how quiet the signal is
against fixed thresholds. So an **in** trim sets how hard the programme drives its pass — how
much of the signal falls into the action region — which is the nearest thing on the panel to
a compressor's threshold knob, and the matching **out** trim is that pass's level
compensation, restoring what the in trim moved. Drive the pass harder with **up in**, take
the level back with **up out**; the same reading holds for the down pair.

In the paired duty their job is calibration: to make the process see, on the way back, the
same level it produced on the way out.

**The calibration matters more than it looks.** If the return arrives 3 dB hot, the
decoder concludes that every quiet passage was louder than it really was, and takes back the
wrong amount in the wrong bands. The result is not a level error you can trim away
afterwards — it is a *spectral* error that moves with the programme.

If nothing in your path changes level, leave all four at **0.0**.

## The resolution of the band split — a design note

**Not a control.** Earlier versions carried a knob here — **exp**, for the exponent it read —
that set how finely the split between the two spectral halves is resolved. It is now a
compiled-in constant of **2048 sections**, and this section stays as the record of why: both
why the number is what it is, and why it is no longer yours to set.

The split itself is fixed at **800 Hz**, three stages above and two below. The system Lime
follows meets its two halves there because two band-defining networks is what the hardware
could be built out of, not because the process wants exactly two. A transform is under no
such limit, so the profile that decides how much of each half's lift a given frequency
receives is sampled onto log-spaced **sections** spanning 20 Hz to 20 kHz — and the count of
those sections is the design parameter with something to say.

Dividing the 8 dB between the two halves' lifts into N sections puts 8/N dB steps into the
encode curve, and an abrupt step from one bin to the next is the one thing this system cannot
decode exactly. Measured through separate encoder and decoder instances on broadband noise at
−26 dB at 48 kHz, the null against the section count runs:

| sections | 4 | 8 | 64 | 128 | 256 | **2048** | 16384 |
|---|---|---|---|---|---|---|---|
| null | −30.0 dB | −39.1 dB | −64.0 dB | −66.9 dB | −67.2 dB | **−67.5 dB** | −67.5 dB |

The null improves as the grid gets finer, plateaus from 256 sections, and never gets better
than **−67.5 dB**. That table is the whole argument for retiring the knob: from the old
default of 256 upward, every setting measured the same and sounded the same, and the settings
that stood apart were the coarse ones, which stood apart by being worse — the coarsest cost
about 37.5 dB of null against the plateau, which is audible. A control whose entire audible
range is ways to make the round trip worse is not worth panel space. So the count is fixed at
**2048**: on the plateau, at the best figure the range produced, and never asking for more
sections than the low half's analysis has bins at any rate Lime runs at.

### Two settings that were never offered, kept as part of the record

The retired knob's own range was the product of measurement, and two coarser settings had
already been removed from it before the knob went. They are still worth recording.

**One section** meant a single weight across the whole band — no frequency-dependent
allocation at all. It nulled at −67.5 dB, and it did so *precisely because it does nothing*:
a flat weight is the smoothest profile there is, with no step anywhere for the decode to fail
on. A position that wins the measurement by not working is a trap.

**Two sections** was the reference machine's own arrangement — one split at 800 Hz, three
stages above and two below — and it is the **worst** setting ever measured here, at −24.0 dB,
because that is where the entire 8 dB of allocation arrives as a single step. So the
reference machine's exact configuration is not reachable in Lime, and that is stated plainly
rather than glossed: what you cannot do is hear precisely what two band-defining networks
sound like. What decided it is that the setting was not so much a faithful option as a broken
one, on a round trip that reaches −67.5 dB.

### What the fixed count means in use

Nothing, and that is the point. An encode and its decode must agree on the resolution exactly
as they must agree on levels, and a constant agrees by construction — across instances,
across sessions, at every sample rate. The one setting the old control demanded you keep
written down with an encoded file is now impossible to get wrong.

The plot used to mark section boundaries with faint vertical ticks, drawn only at 64 sections
and below — denser than that, the boundaries land closer together than the pixels between
them and read as a wash rather than as structure. At the fixed 2048 the plot draws none.

## set up — not on the panel

**Set Up is not a panel control.** It is a parameter, reachable from your host's parameter
list or an automation lane, and it is left out of the interface deliberately: it exists to
align an *unknown* channel, and across a bit-transparent path between two instances there is
nothing to align. If you are working entirely inside a DAW you will never need it.

It is documented here because it is still there and still works, and because if you do put
hardware in the middle it is the fastest way to match the two ends.

Set Up sends a **calibration signal** in place of the programme, so an unknown channel can be
aligned.

- With **b**: shaped noise, 15 dB below reference level. Deliberately not at reference — at
  reference level the noise would risk saturating the extremes of the channel, and what came
  back would not tell you the truth about its response.
- With **a**: a steady tone at reference level.

On the way back — **down** mode with **set up** engaged — Lime recognises its own calibration
noise and starts **auto compare**: every four seconds it alternates between the channel return
and an internally generated reference. Two sounds, alternating; your job is to make them sound
the same using the trims. The panel tells you which one you are hearing at any moment.

The internal reference is *uninterrupted* and the calibration noise carries brief nicks, which
is how the two are told apart — by ear and by the detector.

## input and output

Master gains, ±24 dB, outside everything the plugin does. **input** is the first thing in the
signal path and **output** the last, so both are heard in bypass as well — they are channel
trims rather than a drive control on the process.

They are not the calibration trims and should not be used as them. The trims exist to make
the process see the same level on the way back that it produced on the way out; these exist
to match the plugin to whatever is either side of it in the session.

## s/c hpf and s/c lpf

A high-pass and a low-pass filter on the **level detection** — the classic side-chain
pair, each with its own switch on the **s/c filters** line. They filter what every
compressor stage *hears*, never the audio: no signal, wet or dry, passes through either.
Each weighting is a third-order Butterworth response, 18 dB per octave, −3 dB exactly at
its corner, and both corners run from 20 Hz to 20 kHz; engaged together, the two multiply
and band-limit the detection.

**A filter does nothing until its switch is on**, and both switches rest off. A fresh
instance detects with the full band — with both switches off the detection is exactly what
it was in v1.0.0 — and the knobs' resting positions, 80 Hz for the high pass and 6 kHz for
the low pass, each straight up, are corners held in readiness rather than filters in
circuit. Throwing the switch is what puts a filter in the detection path, and the throw is
eased in at the same pace as a corner move rather than stepped.

Below an engaged high-pass corner — or above an engaged low-pass one — the programme reads
quieter to the detectors than it is, and the process answers as it always does when
something reads quiet: it treats it. Raise the high-pass corner in **up** and more of the
low end takes the lift; lower the low-pass corner and the top takes it instead; in
**down** the matching expander pushes the same regions further away. The plot follows
directly, because the gains it draws are the gains the knobs are changing.

**Any setting leaves the round trip exact.** Both passes weight their detection
identically, so whatever the encoder does under one pair of corners the decoder undoes
under the same pair. The one discipline is agreement, as with the trims: an encode and its
decode must see the same filters — switch states and corners both — so for two-instance
work set them the same at both ends.

## mix

**Mix is not a panel control.** It is a parameter, reachable from your host's parameter
list or an automation lane; its knob went to the side-chain pair, and it is left off the
panel deliberately — its working setting for the round trip is its default, and what it
does beyond that is a mastering and parallel-processing job, not a knob to reach for while
calibrating. The capability is intact and unchanged: it blends the processed path against
the dry one, 0 to 200%.

- **0%** — the dry signal, delayed to match. Identical to bypass.
- **100%** — the processed signal. Default.
- **above 100%** — the difference between the two, applied again. Whatever the process is
  doing is exaggerated rather than diluted, which is the quickest way to hear what a subtle
  setting is actually up to.

**Anything but 100% breaks the round trip.** A partly encoded signal is not something a
decoder can undo, so for real two-instance work this stays at 100.

## polarity

Inverts the output. Applied after **mix**, so it flips everything the plugin sends out.

It is ramped through zero over 75 ms rather than switched, so it does not click.

## bypass

Passes the input through, delayed to the reported latency so nothing moves in time. It is
crossfaded over 75 ms rather than switched, and at full bypass the output is the input
**exactly**, sample for sample — with **input** and **output** at 0.0 dB, since those are
outside the bypass.

Use this rather than your host's bypass if you want the comparison to stay time-aligned —
and engage **auto** as well if you want it loudness-aligned.

## auto

Loudness-matching makeup on the wet path, before the blend.

Any process that changes dynamics changes loudness, and an A/B between two loudnesses is a
loudness test — the louder side wins, whatever else is true of it. **Comparisons are loudness
tests unless the levels match**, and this is the switch that makes them match. Engaged, Lime
runs the power of both paths — the processed signal and the delayed dry copy — through a
3-second integrator and applies the square root of their ratio as makeup on the wet path,
clamped to ±12 dB. The application of the correction is smoothed over 400 ms, so it settles
rather than pumps; below −70 dBFS the correction holds its last value through the silence
rather than chasing a level that is not there; and disengaging the switch returns the gain
to unity.

Because it measures against the delayed dry copy, **bypass** comparisons and **mix** blends
are level-matched: what you are judging is what the process does, not which side is louder.

## The plot

Two traces:

- **input**, filled — the spectrum arriving.
- **gain**, the line — what the process is applying, right now, at each frequency. It reads
  against the dB scale at the right: symmetric at ±28 dB, labelled every 8 dB from −24 to
  +24, around an emphasised 0 dB centre line.

Boost sits above the centre line and cut below it. In **up** the line rises where the process
is lifting; in **down** the plot shows the applied cut below the line, dipping furthest where
the expander is pushing hardest. The two directions draw as mirror images of one behaviour,
which is what they are.

The plot says so when it is not looking at anything live. With **process** at **off** it shows
a flat line at 0 dB, because that is what the process is applying. Under **bypass**, and with
**process** at **a**, it holds its last frame and greys it: A is a time-domain network with no
spectral analysis behind it, so there is nothing to publish, and a trace that kept animating
under a bypassed plugin would be claiming something untrue.

Watch the line against the fill. In **up**, where the signal is quiet the line sits high;
where something is loud it falls back toward the centre at and around it. That withdrawal is
the entire point of the system.

The curve **tapers at both extremes**, and that is deliberate rather than a limitation. Below
roughly 100 Hz and above roughly 5 kHz the process applies fixed shaping — spectral skewing
and antisaturation — which reduces how much lift survives there. Less noise reduction is
needed at the extremes, and applying it anyway would cost more in channel headroom than it
returned. What the plot shows is the *total* response of the pass in circuit, fixed shaping included,
not just what the compander is doing.

## out meter

Output peak level, 0 to −60 dB, marked every 10 dB and numbered every 20. It turns red at and
above 0 dBFS. The reading at the bottom right of the panel is the same figure as a number.

# Quick start

## To hear the round trip — one instance

1. **process b**, **mode both**.
2. Play something with quiet detail: a decaying reverb tail, a fingerboard, room sound between
   phrases, the tail of a piano chord.
3. Toggle **process** between **off** and **b**.

**off** is the modelled channel with its own noise. **b** is the same channel with the noise
reduction working across it. The difference is the point.

Then switch **mode** to **up** to hear the encoded signal on its own — the lift arrives as an
obvious brightening and level rise of the quiet parts, which is what is protecting them.

For a closer look at a subtle setting, take **mix** past 100%: above it the difference between
the processed and dry signals is applied again, so whatever the process is doing is
exaggerated rather than diluted. Put it back to 100% before doing any real work.

## To use it as an insert — one instance

1. **process b**, **mode up**, on a track or bus.
2. Set **mix** for the amount — 100% is the full lift, below it a blend.
3. Listen for what rises: tails, room, the small events the mix was sitting on.

**down** is the same insert leaning the other way. Neither needs decoding — an insert is a
use in itself, not half of a round trip.

## To pair around a channel — two instances

1. First instance, **before** the channel: **process b**, **mode up**.
2. Second instance, **after** the channel: **process b**, **mode down**.
3. Same process in both. This is not optional.
4. If anything in the path changes level, align with **Set Up** — see below.

## Across an outboard insert

1. Instance one in **up**, feeding your hardware send.
2. Instance two in **down**, on the return.
3. Engage **Set Up** on both, from your host's parameter list — it is not on the panel.
   Instance one now sends calibration noise down the send.
4. On the return instance, **auto compare** starts alternating the return against the internal
   reference every four seconds.
5. Adjust **down in** until the two halves of the cycle sound identical in level and in
   tone.
6. Disengage **Set Up** on both.

The trims fix the level. If the two still differ in *tone* rather than level, that is your
hardware's frequency response and no trim will fix it — but you now know it is there, which is
the other reason the calibration signal is noise rather than a tone.

## Bouncing an encoded file

Nothing special. Instance in **up**, bounce. To play it back, an instance in **down** with
the same process. The encoded file will sound bright and squashed on its own; that is
correct.

# How the process works

This section is not needed to use Lime, but it explains what you are hearing.

## The dual-path arrangement

The signal passes through unmodified on the main path. A side chain derives a small
*differential* component from it, and that component is **added** on the way out and
**subtracted** on the way back:

$$ \text{up:}\quad y = (1 + G(x))\,x \qquad\qquad \text{down:}\quad z = \frac{y}{1 + G(z)} $$

The decode is a feedback arrangement, and that is what makes it complementary rather than
merely approximate. Two consequences matter in use:

- At high level the side chain contributes almost nothing, so loud material passes essentially
  untouched. Distortion and overshoot stay low because the thing doing the work is idle when
  the signal is loud.
- A gain error in the channel produces a tracking error of roughly the same size in dB, about
  30 dB below peak — not a runaway. This is why the system is usable at all across a real
  channel, and why the trims are a refinement rather than a prerequisite.

## Staggering, and why it beats parallel bands

Five stages in **b**, at three thresholds, in series. Because they are in series their transfer
functions multiply, which means their decibel effects add — and it means the *spectral
discrimination* steepens with each stage rather than staying fixed.

Three of the stages sit above the 800 Hz split and two below, so the high half reaches
approximately 24 dB of low-level lift and the low half approximately 16 dB. That ratio is not
arbitrary: high-frequency noise is what a channel is worst at and what the ear finds most
objectionable, so that is where the effort goes.

## Least treatment

Each stage stops acting once the signal in its band is loud enough not to need it, and there is
a further mechanism that stops the gain moving once enough reduction has been achieved.

The upshot is that Lime is *not* doing something to your signal most of the time. It acts on
the quiet, and near a loud sound it withdraws — not just at the exact frequency of that sound
but across the region the sound is masking anyway, because reduction there would be inaudible
and would cost channel headroom.

You can watch this happen. Play a sustained note and watch the gain curve dig a valley around
it while staying high everywhere else.

## Spectral skewing and antisaturation

Two fixed networks sit in the main path, ahead of and after the stages:

- **Spectral skewing** — a 12 kHz low-pass with a shelf, and a 40 Hz high-pass with a shelf.
  These deliberately reduce the system's sensitivity at the extremes of the band, where a
  channel's own response is least predictable.
- **Antisaturation** — fixed shelving above roughly 5 kHz and below roughly 100 Hz, which
  keeps the encoded signal from asking the channel for headroom exactly where a channel has
  least to give.

Together they produce approximately 2 dB of shaping at 5 kHz, 6 dB at 10 kHz, and 10 dB at
25 Hz and 15 kHz. This is what makes the plotted curve taper, and it is why the low-level lift
you actually get is around 22 dB in the upper midrange rather than the full 24.

Both networks are implemented as exactly invertible filters, so the decode undoes them
precisely rather than approximately.

# Latency

Lime is a two-pass spectral process and it has real latency:

| Sample rate | Latency | |
|---|---|---|
| 44.1 kHz | 10240 samples | 232 ms |
| 48 kHz | 10240 samples | 213 ms |
| 88.2 kHz | 20480 samples | 232 ms |
| 96 kHz | 20480 samples | 213 ms |
| 176.4 kHz | 20480 samples | 116 ms |
| 192 kHz | 20480 samples | 107 ms |

It is **reported to the host**, so any DAW with delay compensation lines everything back up
automatically. The panel does not print it — a figure you cannot act on is not worth the space,
and the host already knows it.

**The reported figure never changes when you move a control.** **up** and **down** each need
one pass rather than two, and **bypass** and **off** need none, but every one of them is padded
to the same figure. A plugin whose latency moved when a switch was flipped would force the host
to recompensate mid-session, which hosts handle badly and which would shift your session under
you.

Two practical consequences:

- Do not use Lime on a live monitoring path. 200 ms is not a monitoring latency.
- If your host does not compensate plugin latency on a particular bus type, it will not
  compensate this one either. Check by nulling against a bypassed copy.

# Precision and smoothing

Everything internal is **64-bit double**, including the analysis transforms. If your host
offers double-precision processing the entire path is double end to end; if it offers float,
the input is widened to double, processed, and narrowed on the way out. The behaviour is
identical either way.

Every control reaches the audio through a one-pole filter: **95% of any change lands within
75 ms**. Nothing you do — by hand or from an automation lane — puts a step into the signal.
A restored session is the exception it should be: settings recalled with a session or preset
are in force from the first sample rather than fading in.

The switches cannot be interpolated: **off**, **a** and **b** are different processes rather
than three values of one thing, and so are **up** and **down**. Changing one of them fades the
process out to the dry signal, throws the switch while it is inaudible, and fades back —
**1500 ms** across the whole move, on a raised-cosine curve that starts and ends at rest rather
than at speed. You will hear the programme continuously throughout; what you will not hear is
a click, and you should be able to see the move take its time.

That is deliberately much slower than a control smoother, because a switch is not a knob. At
the speed a knob wants, the change arrives as an event rather than as a transition. **bypass**
is not affected and keeps its own 75 ms, since a bypass that takes most of a second to answer
is useless for the comparison it exists to make.

# Specifications

| | |
|---|---|
| Formats | AU, VST3 (macOS) |
| Channels | mono → mono or stereo → stereo; channels unlinked |
| Internal precision | 64-bit double throughout |
| Host precision | double native, float supported by conversion |
| Sample rates | 11.025 kHz to 192 kHz |
| Latency | reported; 107–232 ms depending on sample rate |
| Processes | A (four bands), B (five staggered stages) |
| Low-level lift, B | about 24 dB high frequencies, about 16 dB low |
| Low-level lift, A | 10 dB flat, rising to 15 dB at 15 kHz |
| Thresholds, B | about −30, −48, −62 dB relative to reference |
| Spectral split, B | fixed at 800 Hz, resolved onto 2048 log-spaced sections |
| Trims | ±10 dB, four |
| Masters | ±24 dB in and out |
| Mix | 0 – 200% |
| Auto gain | correction ±12 dB, 3 s integration, 400 ms application |
| Control smoothing | one-pole, 95% in 75 ms |
| Process and mode switching | raised-cosine, 1500 ms across the change |
| Tail | 250 ms |

## Parameters, for automation

| Parameter | Range | Default |
|---|---|---|
| Process | Off / A / B | B |
| Mode | Up / Down / Both | Both |
| Bypass | off / on | off |
| Up In | −10 … +10 dB | 0.0 |
| Up Out | −10 … +10 dB | 0.0 |
| Down In | −10 … +10 dB | 0.0 |
| Down Out | −10 … +10 dB | 0.0 |
| Set Up (host only) | off / on | off |
| Input | −24 … +24 dB | 0.0 |
| Output | −24 … +24 dB | 0.0 |
| Mix | 0 … 200% | 100.00 |
| Auto | off / on | off |
| Polarity | off / on | off |

# Troubleshooting

**It sounds bright and squashed.**
That is **up** at full strength — an upward compressor lifting everything quiet, which on
dense material at 100% **mix** is a lot. As an insert, back **mix** off until it is a blend.
In the paired duty this is the encoded signal, meant to be decoded rather than listened to;
use **both** to hear the finished round trip from one instance.

**The decoded result is dull, or pumps.**
The two instances disagree about something. Check, in this order: same **process**, and then
levels — align with **Set Up** if there is anything in the path that could have changed
level.

**Nothing seems to happen.**
Check **process** is not **off**, **bypass** is not engaged, and **mix** is not at 0%. If you
are running the pair, check that your channel actually adds noise: across a bit-transparent
path there is nothing for the round trip to remove, and in **both** mode with the modelled
channel you should hear a clear difference.

**Everything is late.**
Latency is reported but your host is not compensating for it. This is a host setting, usually
per bus type. Verify by nulling the output against a delayed dry copy.

**A click when I switch process or mode.**
There should not be one — those switches crossfade over 1500 ms, fading out to the dry signal,
changing, and fading back. If you hear a click, it is worth reporting.

**Switching process or mode takes a moment.**
It takes 1500 ms, deliberately. Off, A and B are different processes rather than three values
of one thing, so they cannot be interpolated; the switch is thrown while the wet path is
silent and the programme carries on through the dry path meanwhile. Bypass is not affected and
still answers immediately.

**The plot says "no signal".**
No audio has reached the analysis yet. It appears once the engine has published a frame.

**The plot is grey and not moving.**
Either **bypass** is engaged or **process** is on **a**. Both are states where what the trace
would show is not what you are hearing, so it holds its last frame rather than misleading you.

# Known limitation

The round trip through **up** and **down** does not yet null perfectly. It should return the
original signal to better than −100 dB; measured, it manages about **−80 dB** on steady tones
and **−51 dB** on dense broadband material near reference level (−67 dB at −26 dB).

The cause is understood and recorded in the project's design notes. In brief: the per-bin gain
this process applies is not exactly invertible through a pair of overlapping analyses — a
spectrum with abrupt bin-to-bin gain changes is not the transform of any signal, so no decoder
can undo it exactly. The remedy is to realise the gain through exactly invertible time-domain
filters instead, driven by the same analysis. That component is built and measured at −234 dB
in isolation; integrating it is outstanding work.

In practice the residual sits far below the noise of any channel the paired use is meant to
defeat.
It is stated here rather than omitted because it is a real shortfall against the design target.

---

*Aptitude Audio | Lime | v1.0.1*
