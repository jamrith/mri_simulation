# How the MPRAGE Pulseq Sequences Work

---

## 1. What Is an MPRAGE Sequence?

MPRAGE (Magnetisation Prepared Rapid Gradient Echo) is the workhorse
sequence for high-resolution T1-weighted brain imaging. It combines two
ideas:

1. **Inversion pulse** — a strong 180° RF pulse that flips all the
   magnetisation upside-down. Different tissues (grey matter, white
   matter, CSF) recover at different rates (their T1 relaxation times),
   so after a carefully chosen delay the contrast between them is
   maximised.

2. **Rapid 3D gradient-echo readout** — immediately after the inversion
   delay, a fast train of small-flip-angle excitations reads out one
   "slab" of 3D k-space. Each excitation acquires one line of k-space.

The sequence then waits for the magnetisation to recover, fires another
inversion pulse, and repeats with the next slab of k-space lines — until
the full 3D volume is filled.

### Timing parameters

| Symbol      | Meaning |
|-------------|---------|
| **TI**      | Inversion Time — delay between the 180° pulse and the centre of the readout train (where the most important k-space lines are acquired). Controls grey/white matter contrast. |
| **TRout**   | Outer TR — total time for one inversion cycle (inversion + readout train + dead time). |
| **TRinner / ESP** | Echo Spacing — time for one excitation + readout inside the train. |
| **TE**      | Echo Time — delay from excitation to the centre of the readout window. |
| **BW**      | Bandwidth per pixel — how fast the ADC samples during readout. Higher BW = shorter readout = shorter TE, but more noise. |

---

## 2. How Pulseq Builds a Sequence

Pulseq is an open vendor-neutral format for describing MRI sequences.
The MATLAB toolbox (`extern/pulseq/matlab/+mr/`) provides building
blocks:

### Building blocks

- **`mr.makeBlockPulse`** — creates a rectangular RF excitation pulse.
- **`mr.makeAdiabaticPulse`** — creates the 180° inversion pulse. Uses
  a *hyperbolic secant* waveform designed to invert magnetisation
  uniformly even if the B1 field is imperfect. The waveform is computed
  by `sigpy` (a Python library), which is why the `mrisim` conda
  environment is needed.
- **`mr.makeTrapezoid`** — creates a trapezoidal gradient waveform
  (ramp up → flat top → ramp down). Used for readout gradients,
  phase-encoding gradients, and spoilers.
- **`mr.makeAdc`** — defines when the scanner's analogue-to-digital
  converter samples the signal.

### Assembly

The script builds the sequence block-by-block inside nested loops:

```
for each phase-encoding line (outer loop):
    play inversion pulse
    wait TI delay
    for each partition-encoding line (inner loop):   ← "turbo factor" iterations
        play RF excitation
        play readout gradient + ADC
        play phase/partition encoding gradients
        play spoiler gradients
    end
    wait until TRout is reached
end
```

Each call to `seq.addBlock(...)` appends one time-block to the sequence.
A 320×256 matrix produces 320 outer × 256 inner = 81,920 readout blocks
(plus inversion/delay blocks).

### RF spoiling

After each excitation, the RF phase is incremented by a fixed amount
(84°). This prevents coherent transverse magnetisation from building up
between TRs — without it, the signal would contain unwanted stimulated
echoes that corrupt the T1 weighting. The ADC phase tracks the RF phase
so the receiver demodulates correctly.

### Gradient spoiling

After each readout, a strong gradient "spoiler" is played along the
readout axis (typically 3× the maximum k-space excursion). This
dephases any remaining transverse magnetisation across the voxel,
complementing the RF spoiling.

---

## 3. How TE Is Determined

The echo time is the interval from the centre of the excitation RF pulse
to the centre of the ADC readout window. In a symmetric-echo Pulseq
MPRAGE, it breaks down as:

```
TE = rf_duration/2 + ringdown + prewinder_duration + gradient_rise + readout_duration/2
```

| Component             | Typical value | What sets it |
|-----------------------|---------------|--------------|
| rf_duration/2         | 0.05 ms       | Hard pulse length (100 µs) |
| RF ringdown           | 0.02 ms       | Hardware dead-time after RF transmit |
| Prewinder duration    | ~0.8 ms       | Must ramp a gradient to shift k-space from 0 to −kmax before readout |
| Gradient rise time    | ~0.05 ms      | Readout gradient ramp, set by slew rate |
| readout_duration/2    | ~2.4 ms       | Half the ADC sampling window; set by BW and matrix size |

The dominant term is **readout_duration/2**. With BW = 208 Hz/pixel and
320 readout samples, the full readout lasts 4.8 ms, so the echo can
never occur earlier than ~2.4 ms after excitation in a symmetric readout.

### Why our TE (2.9 ms) differs from the HCP protocol (2.14 ms)

The HCP protocol reports TE = 2.14 ms, which is *less than half the
readout duration* (2.4 ms). This is physically impossible with a
symmetric echo — the Siemens product sequence uses **asymmetric
(partial) echo readout**:

- The readout prewinder is reduced so the gradient echo occurs *before*
  the centre of the ADC window.
- The ADC starts sampling before the echo, capturing more of the
  post-echo k-space than the pre-echo side.
- The missing k-space samples at the leading edge are reconstructed via
  conjugate symmetry (partial Fourier in readout).

Our Pulseq scripts use **symmetric echo** — the prewinder positions the
echo exactly at the ADC centre. This is simpler and doesn't require
partial Fourier reconstruction, but sets a TE floor of ~readout/2.

For simulation purposes the difference is minor: slightly more T2* decay
at 2.9 ms vs 2.14 ms, but the T1 contrast from the inversion
preparation (which is the whole point of MPRAGE) is unaffected.

---

## 4. GRAPPA Acceleration

Fully sampling all 320 phase-encoding lines takes a long time. GRAPPA
speeds this up by **skipping lines** in the outer (phase-encoding) loop:

- With acceleration factor R = 2, only every 2nd line is acquired.
- A small central region of k-space (the **ACS** / autocalibration
  lines — 32 in our case) is still fully sampled. These lines are used
  to calibrate the GRAPPA reconstruction kernel, which learns how to
  predict the missing lines from the acquired neighbours.
- The missing lines are filled in computationally during reconstruction.

The acceleration is only applied to the **outer** (phase-encoding) loop.
The inner (partition-encoding) loop is always fully sampled — every
inversion block still plays the full 256 partition lines. This means
GRAPPA R = 2 roughly halves the number of inversion blocks (320 → 176),
cutting scan time nearly in half.

### Labels

In Pulseq, each acquired line is tagged with metadata labels so the
reconstruction software knows what to do with it:

| Label   | Meaning |
|---------|---------|
| `LIN`   | Phase-encoding line number (outer loop index) |
| `PAR`   | Partition-encoding line number (inner loop index) |
| `REF`   | This line is a GRAPPA reference/calibration line |
| `IMA`   | This line contributes to the final image |
| `NOISE` | Noise prescan (no RF, just ADC — measures receiver noise floor) |

Lines in the ACS region that also fall on the regular undersampling grid
get both `REF=true` and `IMA=true`. Lines in the ACS region that would
otherwise be skipped get `REF=true, IMA=false` — used only for
calibration, not the final image.

---

## 5. Timing Raster Constraints

MRI hardware operates on fixed clock rasters. Pulseq enforces three:

| Raster               | Default  | What it constrains |
|-----------------------|----------|--------------------|
| Gradient raster       | 10 µs    | Gradient waveform sample spacing and durations |
| ADC raster            | 100 ns   | ADC dwell time (interval between samples) |
| Block duration raster | 10 µs    | Total duration of every sequence block |

These constraints interact: the ADC dwell time must be a multiple of
100 ns, **and** the total ADC duration (`dwell × N_samples`) must
contribute to a block duration that lands on the 10 µs grid.

### Example: BW = 210 Hz/pixel with 320 readout points

The natural dwell time is:

```
dwell = 1 / (210 × 320) = 14,880.95 ns   ← not on any raster
```

Rounding to the 100 ns ADC raster gives 14,900 ns, but:

```
14,900 × 320 = 4,768,000 ns = 4,768 µs
4,768 / 10 = 476.8   ← not an integer → block raster violation
```

The fix: dwell must be a multiple of **500 ns** — the smallest value
satisfying both the 100 ns ADC raster and the requirement that
`dwell × 320` is a multiple of 10,000 ns.

```
dwell = 15,000 ns   →   BW = 208.3 Hz/pixel (close enough)
15,000 × 320 = 4,800,000 ns = 4,800 µs
4,800 / 10 = 480 ✓
```

The general rule: for N readout samples, the dwell must be a multiple of
`lcm(adc_raster, block_raster / gcd(N, block_raster / adc_raster))`.
