# Release notes — floquet_soc v0.1.7-rc1 (untagged working tree)

**Type:** feature/refactor. Physics core unchanged (still the v0.1.6
wsvec-aware H/spin/Peierls implementation).

## Why

The MoSSe work had put material-specific spin-texture and DFT-comparison
plotting in `materials/MoSSe/scripts/`. That is the wrong layer: the same
operation is needed for PbNV, Pt and any future material. The programs now
live in the framework; material folders keep only data, configs and results.

## New framework programs

### `floquet_soc/observables/spin_texture_map.py`
- `resolve_band(cfg, model, band)` — `vbm` / `cbm` / integer / config default.
- `full_kpts(cfg, n, kmax)` — BZ-plane mesh (`plane_kpts`, PbNV convention).
- `window_kpts(lattice, center_frac, half, n)` — square window around K/Γ in
  1/Angstrom for valley zooms.
- `compute_map(model, spin_at_k, band, kpts, mask_split=...)` — energies,
  all-band spins, selected-band spins, nearest-neighbour gap; suppresses
  gauge-dependent in-plane arrows inside a degeneracy window.
- `plot_map(...)` — one PbNV-style panel (bwr S_z, black arrows/dots,
  S_z colorbar, no title).

### CLI
```bat
floquet.py spin_texture --config CFG --mode full   --band vbm --n 41 --kmax 1.65
floquet.py spin_texture --config CFG --mode window --center 0.333 0.333 0 --half 0.05 --n 41 --mask-split 0.002
floquet.py dft-spin-check --config CFG --dft FILE.dat --band vbm --outdir OUT
```

### `floquet_soc/observables/floquet_points.py`
- `floquet.py floquet_points --config CFG --points K,Kp,Gamma --gamma 0 0.02 0.05 --pol circ+ circ- linear`
- band-resolved Ebar / S_bar / dS / w_main / w_sum / micromotion / overlap;
- K/K' valley asymmetry aggregation; K' derived from -K when the config only
  lists K; degenerate matches flagged `reliable=false`;
- writes `floquet_points.json` plus spin/valley figures.
- `--rashba`: effective Gamma->K valence-doublet slope alpha_eff(gamma);
- `--grid`: small K/K' grids with Ebar/S_z/|S_xy| per point;
- `materials/MoSSe/scripts/floquet_round2.py` removed (framework command
  replaces it).

### `floquet_soc/validation/dft_spin.py`
- `read_spin_texture_dat` (kx ky E Sx Sy Sz, `#Grid-Size`, `#Band-Index`).
- `compare` (energy offset/residual; raw spin RMS/correlation; global scale
  fit; normalized-direction RMS with degeneracy mask).
- `save_and_plot` (JSON/NPZ + two-panel texture + scatter/energy figure).

## Material layer

- `materials/MoSSe/scripts/phase1_figures.py` now calls
  `spin_texture_map.full_kpts/window_kpts/compute_map/plot_map`; no duplicate
  spin-texture plotting code remains.
- `materials/MoSSe/scripts/dft_spin_compare.py` deleted; use
  `floquet.py dft-spin-check`.
- Materials keep: data, `configs/*.json`, `phase0_probe.py`,
  `phase1_static.py`, `phase1_figures.py` (bands/Rashba), `floquet_trial_K.py`,
  `floquet_round2.py`.

### `floquet_soc/core/provenance.py` + `floquet.py provenance`
- version, config SHA256, every existing input file size/SHA256, spin source,
  numerics and drive settings;
- injected into `floquet_points.json` and `dft-spin-check` JSON.

### `docs/`
- `REQUIREMENTS_GAP_MATRIX.md`: requirement-by-requirement PASS/WARN/GAP;
- `PHASE_A_ACCEPTANCE.md`: input/spin-basis facts, checks, WARN/GAP.

### New SOC/spin workflows
- `floquet.py floquet_bands`: static E(k) vs Floquet spectral weight A0(k,E)
  (shared axes, Gamma at the boundary, high-symmetry lines).
- `floquet.py soc_scan`: standardized alpha_R(gamma), Delta_spin(gamma),
  PST metrics (collinearity/axis/rotation) in one JSON/CSV/PNG.
- `floquet_points --subspace`: degenerate-subspace averaged spin for
  PbNV/Pt; single-band labels are no longer used there.
- Benchmarks 1 and 5: `spin_degenerate`, `magnetic_soc` (both PASS).

### Phase A/B acceptance
- `docs/PHASE_A_ACCEPTANCE.md`, `docs/PHASE_B_ACCEPTANCE.md`;
- three-material validate: PbNV 27/0/2, Pt 25/0/4, MoSSe 28/0/1;
- Floquet-Rashba benchmark 8 PASS / 0 FAIL;
- MoSSe band-edge Floquet spin PASS; PbNV/Pt marked WARN (tiny gaps /
  metallic degeneracies) — single-band labels are not used as physics.

## Validation

- New `tests/test_spin_texture_map.py` (4 tests: band resolution, grids,
  map/mask, DFT read/compare/save) — pass.
- `test_version_string` against VERSION=0.1.7-rc1 — pass.
- Frozen `v0.1.6-rc1` tag/snapshot is untouched; this is an untagged working
  tree until a new release is requested.
