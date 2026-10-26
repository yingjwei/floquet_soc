# Phase 2 report: Pt physics minimal fixes

Date: 2026-10-22/23. Scope: minimal fixes only.

## Before / after key numbers

| quantity | before | after |
|---|---|---|
| Pt off-grid Hermiticity | **1.27 eV** | **1.81e-14 eV** |
| Pt off-grid Kramers residual (median/max) | 67 meV / 0.89 eV | **0.43 meV / 4.54 meV** |
| Pt on-mesh Kramers residual | 7.9e-8 eV (model) | 4.6e-14 eV |
| spin operator off-grid Hermiticity | 0 (masked by symmetrisation) | 6.7e-16 (unmasked, WS interp) |
| SOC splitting for Pt | 48.8/62/126 meV "spin-paired" | NOT_APPLICABLE; delta_PT 4.37 meV max (path) |
| Pt PST | M=(0.609,0.277,0.115) claimed | NOT_APPLICABLE: coll 0.368-0.413, axis dot 0.345, invariant spin 0.0024 |
| Pt band-vs-reference | PASS (nearest metric, hidden outliers) | WARN (aligned: median 54.1 meV, max 1141.8 meV) |
| m_sp convergence observable | dEbar=dS=0 (m_sp insensitive) | dWsum 1.34e-08 at m_sp=16, sum-rule error 1.0 (genuine m_sp dependence) |

## Pt benchmark after fix

PASS=14, WARN=3, FAIL=0. Key PASS rows: off-grid Hermiticity, on-mesh
Kramers, off-grid Kramers (median<1 meV, max<10 meV), SOC NOT_APPLICABLE.
WARN rows: hr fallback, PST gauge, band reference.

> 冻结说明（2026-10-23）：这是 phase-2 中间结果；冻结版最终 Benchmark 7 为
> PASS=13 / WARN=5 / FAIL=0（`pt_kramers_offgrid` 降级为 WARN，新增
> `pt_quantitative_status`），以 `RELEASE_VALIDATION.md` §7 与
> `outputs/benchmarks/pt/report.md` 为准。

## PbNV (kept)

* TRIM pair splittings 0.309/0.386/0.487/0.412 meV (valence), consistent with
  the reported 0.31-0.49 meV; sub-meV residual = Wannier truncation.
* PST M=(0.9570,0.0401,0.0029), collinearity 0.9570, axis 0.9999 z; gauge
  test axis dot 0.998 (stable), collinearity spread 0.038 -> kept with
  subspace caveat.
* Per-band spin statistics in `PT_PHYSICS_AUDIT.md` (all four bands
  Ising-like, |Sz|>0.8 over >96% of the mesh).

## Modified files

* `floquet_soc/core/wannier.py` — `ws_supercell()`; `from_chk_eig` uses the
  WS/ndegen reconstruction for unshifted meshes; stores `interpolation`.
* `floquet_soc/core/spin_wannier.py` — WS/ndegen spin interpolation, no
  forced Hermitisation, `hermiticity_offgrid()`, off-grid check, lattice arg.
* `floquet_soc/observables/kramers.py` — new: delta_PT, inter-doublet gaps.
* `floquet_soc/observables/pst.py` — `gauge_invariance_test()`.
* `floquet_soc/observables/soc_splitting.py` — docstring: generic gaps only.
* `floquet_soc/cli.py` — PT-aware `soc` (NOT_APPLICABLE), Pt PST
  NOT_APPLICABLE + gauge data, PbNV PST gauge check, converge dWsum output.
* `floquet_soc/workflows/converge.py` — m_sp-sensitive `dev_wsum` and
  production criterion.
* `floquet_soc/benchmarks/pt.py` — Hermiticity/Kramers/SOC/PST rows,
  index-aligned band comparison as WARN.
* `configs/pt.json`, `configs/pbnv.json` — `system` classification.
* `tests/test_pt_physics.py` — new regression tests (10 tests).
* reports: `PT_PHYSICS_AUDIT.md`, `PT_KRAMERS_REPORT.md`,
  `PT_GAUGE_INVARIANCE_REPORT.md`, `PT_HAMILTONIAN_REPORT.md`,
  `PT_PHASE2_REPORT.md`; benchmark report regenerated.

## Tests

`python -m unittest tests.test_pt_physics` -> 10 tests OK
(2 data-dependent Pt tests skipped if `benchmarks_data/Pt` is absent).
Affected existing modules (spin_wannier, stage_c) were re-run; see the final
summary in the delivery message.

## Remaining limitations (explicit)

1. Official Pt `wannier90_hr.dat`: BLOCKED (no executable;
   `restart=plot`/`write_hr` feasible but not runnable here; also not in the
   public tutorial data).
2. Pt `Pt_band.dat` reference comparison: WARN only (setting mismatch
   suspected: `use_ws_distance`).
3. Off-grid Kramers residual 0.43 meV median is Wannier truncation of the
   4x4x4 mesh; it is below the ~50 meV physical SOC scale but is not
   machine-zero, and is reported as such.
4. PbNV PST remains approximate: axis stable (0.998) but collinearity has a
   3.8% gauge spread near TRIM; quoted with the caveat.
5. No bcc Fe benchmark, no new materials, no PbNV photo-gap work in this
   phase.
