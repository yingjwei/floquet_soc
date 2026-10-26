# Pt Physics Audit — final status (phase 2)

Scope: minimal fixes for Pt Benchmark 7 only. No new materials, no parameter
scans, Floquet solver core untouched.

## What was found (phase 1) and what was fixed (phase 2)

| issue | before | after | status |
|---|---|---|---|
| off-grid Hermiticity of chk+eig H(k) | **max |H-H^dag| = 1.27 eV** | **1.81e-14 eV** | FIXED |
| Kramers residual on DFT mesh | 7.9e-8 eV (model), 2e-10 eV (DFT) | 4.6e-14 eV (model), 2e-10 eV (DFT) | PASS |
| off-grid Kramers residual (400 random k) | median **67 meV**, max **0.89 eV** | median **0.43 meV**, max **4.54 meV** | PASS (interp-truncation level) |
| spin operator mesh interpolation | 1e-14 | 1.44e-15 | PASS |
| spin operator off-grid Hermiticity | masked by forced symmetrisation | **6.7e-16** without masking | PASS |
| '48.8 / 62 / 126 meV SOC splittings' | mislabelled inter-doublet + spurious within-pair gaps | `soc` reports delta_PT and inter-doublet gaps; no SOC spin splitting for PT systems | FIXED (NOT_APPLICABLE) |
| Pt PST | M=(0.609,0.277,0.115) claimed as physical | gauge test on true DFT Kramers pairs: collinearity 0.368-0.413 (base 0.397), axis dot min 0.345; invariant subspace spin mean 0.0024 | FIXED (NOT_APPLICABLE) |
| official Wannier90 hr.dat | absent | still absent (no binary, repo has no hr.dat) | BLOCKED / documented |
| benchmark band comparison | PASS with nearest-reference metric (hidden outliers) | WARN, index-aligned: median 54.1 meV, max 1141.8 meV; reference may use different settings | FIXED (honest WARN) |

Root cause of the interpolation error: `WannierModel.from_chk_eig` used the
R=0..N-1 alias set with `ndegen=1`. That set is not symmetric under R -> -R,
so the interpolant is non-Hermitian off grid (1.27 eV) and breaks PT (67 meV
median). The fix uses the Wannier90 WS/ndegen convention (R=-N/2..N/2,
ndegen = residue-class multiplicity), which was independently validated on
PbNV against the official `wannier90_hr.dat` (R/ndegen exact, interpolation
agreement 8.9e-6 eV limited by the hr text precision).

## Answers to the six audit questions

1. 48.8/62/126 meV are **static** (no drive in `cmd_soc`).
2. Band indices are **0-based** numpy columns.
3. Formula: `adjacent_splitting` = E_{i+1}-E_i; the old
   `spin_paired_splitting` greedily paired nearest-energy bands preferring
   antiparallel spins.
4. Yes: 48.8 meV = inter-Kramers-doublet spacing; 62.0/126.1 meV = spurious
   within-pair splittings from the broken interpolation.
5. Yes: the alias interpolation broke PT off-grid (now fixed to ~0.4 meV
   median, 4.5 meV max, which is Wannier-truncation level, not physics).
6. No physical mechanism: Pt is non-magnetic and centrosymmetric; DFT mesh
   degeneracy is 2e-10 eV.

## Pt Benchmark 7 after the fix

Summary: **PASS=14, WARN=3, FAIL=0** (report:
`outputs/benchmarks/pt/report.md`). PASS rows include off-grid
Hermiticity, on-mesh Kramers, off-grid Kramers (tolerance median<1 meV,
max<10 meV; 4x4x4 truncation, target SOC scale ~50 meV), and explicit
SOC/PST NOT_APPLICABLE rows. WARN rows: fallback hr reconstruction, PST gauge,
band-vs-reference comparison.

> 冻结说明（2026-10-23）：本表是 phase-2 中间结果；冻结版最终 Benchmark 7 为
> **PASS=13 / WARN=5 / FAIL=0**（`pt_kramers_offgrid` 在最终验收中降级为 WARN，
> 另加 `pt_quantitative_status`），以 `RELEASE_VALIDATION.md` §7、
> `outputs/benchmarks/pt/report.md` 与 `STATUS_REGISTRY.md` 为准。

## PbNV (not overturned)

* TRIM valence-pair splittings (official hr model): Gamma 0.309 meV,
  X 0.386 meV, Y 0.487 meV, S 0.412 meV; conduction pair 0.11-0.32 meV.
  This matches the ~0.31-0.49 meV statement in `bands.npz`; the residual is
  Wannier-interpolation truncation at TRIM, not a physical PT breaking.
* PST (bands 110-113, 323 k): M=(0.9570,0.0401,0.0029), collinearity 0.9570,
  axis ~z (stable: gauge test axis dot min 0.998; collinearity spread 0.038).
  Keep as an approximate physical result with the subspace caveat.
* Per-band spin statistics (effective momentum region = fraction of the
  19x17 mesh where |S|>0.5):

| band | mean |S| | mean Sz | frac |S|>0.5 | frac |Sz|>0.8 |
|---|---|---|---|---|
| 110 | 0.973 | -0.031 | 100.00% | 94.12% |
| 111 | 0.973 | +0.031 | 100.00% | 94.12% |
| 112 | 0.529 | +0.014 | 53.25% | 25.08% |
| 113 | 0.629 | -0.013 | 66.87% | 36.84% |

All four bands are predominantly out-of-plane (Ising-like), consistent with
PST axis ~z.

## Acceptance criteria

* Pt Hermiticity on/off grid: PASS (1.8e-14 eV).
* Static Pt Kramers on the validated model: PASS (mesh 4.6e-14 eV; off-grid
  median 0.43 meV = interpolation truncation).
* SOC module no longer names fake splittings as spin splitting: PASS
  (NOT_APPLICABLE branch, config-based system classification).
* Degenerate-system PST does not depend on single-state gauge: PASS
  (NOT_APPLICABLE + gauge-invariance test on true DFT Kramers pairs).
* Inapplicable metrics return explicit NOT_APPLICABLE/WARN: PASS.
* Official Wannier90 band comparison: **BLOCKED** (no executable; hr.dat not
  in the public dataset) — documented, not hidden.
