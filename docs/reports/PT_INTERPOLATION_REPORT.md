# Pt Wannier interpolation comparison (final acceptance round)

## Reference data checks

* `Pt_band.kpt`: header `593`, 593 k points, max|k| = 0.75 -> fractional
  coordinates (confirmed: interpreting them as 1/A makes the comparison
  worse by 8-25x). Path W-L-G-X-W-G.
* `Pt_band.dat`: 18 bands (num_wann), 593 path points, energy range
  7.6406375-42.004637 eV.
* At Gamma (path point with k=0):
  * reference vs model: **0.000 meV** (18 bands),
  * reference vs DFT eigenvalues (lowest 18 of 24): **3617.038 meV**,
  * model vs DFT eigenvalues: **3617.039 meV**.
  => the reference and the chk/eig data agree on the DFT mesh at Gamma;
  the discrepancy is between mesh points (interpolation scheme).
* `Pt.wout` has no `use_ws_distance` record; there is no `Pt.wsvec`.
  Wannier centres are present in the chk (used for the test below).

## Interpolation conventions vs the reference (index-aligned, 593 k)

| convention | median | max |
|---|---|---|
| old fallback: min-norm R representatives, ndegen=1 | 39.50 meV | 773.80 meV |
| raw alias R=0..N-1, ndegen=1 | 731.71 meV | 15822.2 meV |
| WS n1 (-N/2..N/2), ndegen=class (current) | 54.13 meV | 1141.8 meV |
| R=-N/2..N/2-1, ndegen=class | 82.62 meV | 1425.8 meV |


Note: the *old fallback* used min-norm representatives (roughly
R=-N/2..N/2-1) with ndegen=1. It reproduces the mesh but is
non-Hermitian off-grid (1.27 eV) and breaks PT; its smaller deviation
from Pt_band.dat is accidental, not evidence of correctness.

None of the HR-only reconstructions reproduces `Pt_band.dat`; the old alias
set is numerically closer to the reference but is non-Hermitian off-grid
(1.27 eV) and breaks PT, so it cannot be used as a validation target. The
residual 40-55 meV median is too large to be explained by numerical round-off
and points to a **different Wannier setup** in the reference run
(disentanglement window, `use_ws_distance`, or centres-dependent WS
corrections) rather than a simple k-unit/band-index error.

A centre-phase gauge test (H_nm -> e^{ik(tau_m-tau_n)}H_nm) leaves the
eigenvalues exactly invariant, so eigenvalues alone cannot identify
`use_ws_distance`; its effect enters through the R-vector assignment and would
require the original hr.dat/wsvec data.

## Official hr.dat

Still **BLOCKED** on this Windows machine: no wannier90 executable, and the
public tutorial path contains no `hr.dat`. The input set is compatible with
`restart = plot` + `write_hr = .true.` (Pt.win currently has
`#restart = plot` commented and no write_hr), so this is an environment issue,
not a data incompatibility. Exact recipe for a Linux/HPC run:

```
cp Pt.win wannier90.win
# add: restart = plot   and   write_hr = .true.
wannier90.x wannier90
# -> wannier90_hr.dat (official interpolation reference)
```

## Status of Pt

Until that official interpolation comparison passes, Pt is classified as
**software-interface / error-detection benchmark**, not a validated
quantitative SOC/Floquet material. `pt_kramers_offgrid` stays **WARN**.
