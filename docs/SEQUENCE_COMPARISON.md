# Sequence Comparison

## 1. Demo Sequence vs HCP Sequence

| Parameter | Demo (`writeMPRAGE`) | HCP (`mprage_hcp_260321`) |
|---|---|---|
| Matrix | 192 x 240 x 256 | 256 x 320 x 320 |
| FOV | 192 x 240 x 256 mm | 179.2 x 224 x 224 mm |
| Voxel size | 1.0 x 1.0 x 1.0 mm | 0.7 x 0.7 x 0.7 mm |
| Flip angle | 7 deg | 8 deg |
| TI | 1100 ms | 1000 ms |
| TRout | 2500 ms | 2400 ms |
| BW | 200 Hz/px | 208 Hz/px |
| Readout oversampling | 1x | 1x |
| RO spoiling | 3x kmax | 3x kmax |
| RF spoiling | 84 deg | 84 deg |
| MaxGrad | 24 mT/m | 30 mT/m |
| MaxSlew | 100 T/m/s | 130 T/m/s |
| Orientation | Sagittal | Sagittal |
| GRAPPA (grappa variant) | R=2, 32 ACS | R=2, 32 ACS |

The sequences are structurally identical — same MPRAGE design, same
spoiling strategy, same k-space traversal order. The demo is a generic
lower-resolution MPRAGE; the HCP version pushes to 0.7 mm isotropic
with slightly more aggressive timing (shorter TI/TR) and higher gradient
performance.


## 2. Our HCP Sequence vs the Actual HCP Protocol

| Parameter | HCP Protocol | Our Sequence | Match? |
|---|---|---|---|
| Orientation | Sagittal | Sagittal | Yes |
| Voxel size | 0.7 x 0.7 x 0.7 mm | 0.7 x 0.7 x 0.7 mm | Yes |
| FOV read | 224 mm | 224 mm | Yes |
| FOV phase | 100% (= 224 mm) | 224 mm | Yes |
| Slices per slab | 256 | 256 | Yes |
| Base resolution | 320 | 320 | Yes |
| Phase resolution | 100% | 100% | Yes |
| TR | 2400 ms | 2400 ms | Yes |
| TI | 1000 ms | 1000 ms | Yes |
| Flip angle | 8 deg | 8 deg | Yes |
| GRAPPA R | 2 | 2 | Yes |
| ACS lines | 32 | 32 | Yes |
| RF spoiling | On | On (84 deg) | Yes |
| Turbo factor | 256 | 256 (inner loop) | Yes |
| BW | 210 Hz/px | 208.3 Hz/px | ~1% off (raster rounding) |
| ESP | 7.6 ms | 7.66 ms | ~1% off (follows from BW) |
| TE | 2.14 ms | 2.9 ms | No — see below |
| Inversion type | Non-selective | Non-selective (adiabatic hypsec) | Yes |
| Fat suppression | Water excitation (fast) | None | No — not implemented |
| Phase oversampling | ~10% (352 PE lines) | 10% (352 PE lines) | Yes |
| Partial Fourier | Off | Off | Yes |
| Scan time (GRAPPA) | 7:40 | 7:41 | Yes |


### What matches well

Geometry, resolution, contrast timing (TI, TR, flip angle), acceleration
(GRAPPA R=2, 32 ACS), turbo factor, orientation, and RF spoiling are all
faithful to the HCP protocol. These are the parameters that dominate
image contrast and spatial resolution.


### What differs

**TE (2.9 vs 2.14 ms)** — Our sequence uses symmetric echo (the
gradient echo is centred in the ADC window). The Siemens product
sequence uses asymmetric (partial) echo, where the echo occurs before
the ADC centre, allowing TE to be shorter than half the readout
duration. This adds slightly more T2* weighting to our sequence but does
not affect the T1 contrast from the inversion preparation, which is the
primary source of contrast in MPRAGE.

**Fat suppression** — The HCP protocol uses "water excitation (fast)", a
binomial or spectral-spatial RF pulse that excites only water protons
and suppresses fat signal. Our sequence uses a simple hard (block) pulse
that excites both water and fat. For brain simulation with KomaMRI's
brain phantom (which does not model fat), this has no effect.

**BW (208.3 vs 210 Hz/px)** — Forced by the timing raster constraint
(dwell must be a multiple of 500 ns for a 320-point readout to satisfy
the 10 us block duration raster). The 1% difference is negligible for
image quality and SNR.

**Phase oversampling and scan time** — The HCP protocol uses ~10% phase
oversampling, extending the phase-encoding dimension from 320 to 352
lines. With GRAPPA R=2 and 32 ACS lines this gives 192 sampled PE
steps (176 undersampled + 16 extra ACS), and 192 x 2.4 s = 460.8 s =
7:41, matching the reported scan time of 7:40. Our sequence now includes
this oversampling. The oversampled lines are discarded after
reconstruction (they exist to avoid aliasing artefacts at the FOV
edges).
