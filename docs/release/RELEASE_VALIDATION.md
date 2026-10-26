# RELEASE VALIDATION — floquet_soc v0.1.0-rc1

Date: 2026-10-23. Release class:

**Benchmarked Wannier–Floquet–SOC research framework — validated numerical
core, material-specific quantitative validation required.**

## 1. Version and environment

| item | value |
|---|---|
| version | **v0.1.0-rc1** (file `VERSION`, `floquet_soc.version.__version__`) |
| CLI | `python floquet.py --version` -> `floquet_soc 0.1.0-rc1` |
| OS | Windows (10.0.26200) |
| Python | 3.12.10 |
| numpy / scipy / matplotlib | 2.5.3 / 1.18.1 / 3.11.2 |
| dependencies | numpy, scipy, matplotlib only (no compiled extensions) |
| tests | `python -m unittest discover -s tests` |

## 2. Input data provenance

| data | source | files | notes |
|---|---|---|---|
| PbNV (Benchmark 6) | user VASP + Wannier90 run, `D:\Floquet\PbNV` | win/hr/chk/eig/spn(formatted LWRITE_SPN)/centres | official `wannier90_hr.dat` present; old-run backup `C:\Users\15756\.dsh\PbNV_oldrun_backup_2026-10-22` |
| Pt (Benchmark 7) | `wannier-developers/wannier-tutorials`, `2022_05_Trieste/DAY2_WBERI/advanced/spin_hall/data_Pt` | Pt.win/chk/eig/spn/nnkp/wout/band.dat/band.kpt | **no official hr.dat**; chk+eig WS/ndegen reconstruction |
| analytic models | code (`benchmarks/rashba.py`, `rashba_dresselhaus.py`, `floquet_rashba.py`) | — | parameter-free validation |

## 3. Validated numerical capabilities (PASS)

* Wannier Hamiltonian construction: `from_files` (hr), `from_chk_eig`
  (WS/ndegen reconstruction for unshifted meshes), on-mesh exactness.
* Hermiticity: Pt H(k) off-grid `max|H-H†| = 1.81e-14 eV`;
  spin operator mesh `1.4e-15`, off-grid `6.7e-16` (no forced symmetrisation).
* Kramers diagnostics (`observables/kramers.py`): δ_PT, inter-doublet gaps.
* SOC splitting: generic band gaps + explicit system classification
  (PT non-magnetic → SOC splitting NOT_APPLICABLE); spin-paired splitting only
  for non-centrosymmetric systems.
* PST: metric + gauge-invariance test (random U(2) inside near-degenerate
  subspaces, gauge-invariant subspace spin).
* Floquet: Peierls drive, midpoint evolution, Sambe/van Vleck, spin-resolved
  spectra, micromotion spin average; convergence with photon-weight deficits.
* Stage-D convergence: n_step and m_sp with per-state deficit `D_i`.

## 4. Physics models passed (PASS)

* Analytic Rashba: Δ(k)=2α|k| and spin ⟂ k (machine precision).
* Rashba+Dresselhaus: PST eigenvalues/axis analytic (dev <1e-8).
* Floquet-Rashba: zero-drive limits, photon-weight normalization, Sambe vs
  evolution, spin-resolved spectral identity.
* PbNV: 156 Wannier, 19×17=323 k spin texture; SOC spin-paired splittings;
  PST approximate (see WARN below); TRIM T-Kramers residuals 0.309/0.386/
  0.412/0.487 meV (interp truncation); validate 26 PASS / 2 WARN / 0 FAIL.

## 5. Approximate support (WARN, explicit)

* **PbNV PST**: M=(0.9570,0.0401,0.0029), collinearity 0.9570, axis ≈ z.
  Gauge test: axis dot min 0.998, collinearity spread **0.038**; bootstrap
  0.9421–0.9730; near-Γ region 0.865 vs far region 0.964. Approximate
  Ising-like result only; not a universal PST standard.
* **Pt off-grid Kramers residual**: median 0.43 meV / max 4.54 meV =
  numerical interpolation error, `pt_kramers_offgrid = WARN`.
* Photon truncation: m_sp=16 gives max per-state deficit 1.0 (2 states);
  m_sp≥32 converged. `m_sp=32` is a **test-parameter convergence choice for the
  current Pt/PbNV settings, not a universal default**; each material must run
  the deficit check (`converge.npz: w_deficit_*`, `accepted_m_sp`).

## 6. Not quantitatively validated (do not use for production claims)

* **Pt**: software-interface / error-detection benchmark. Official
  Wannier90 `hr.dat` unavailable; `Pt_band.dat` not reproduced
  (index-aligned median 54.1 meV, max 1.14 eV); off-grid Kramers residual is
  numerical. Do not use Pt for quantitative SOC/Floquet material results.
* Any material with shifted k mesh: chk+eig alias fallback is
  grid-consistent only; off-grid physics NOT validated (spin operator check
  says `NOT_VALIDATED`).
* Photon-resolved spectra for states with `D_i > 1e-3`: UNCONVERGED for that
  m_sp window; must increase m_sp or exclude the state.

## 7. Status-token policy (CLI + metadata)

The release keeps non-PASS states visible; they are never hidden by a global
PASS count:

| token | where preserved | meaning |
|---|---|---|
| `WARN` | validate.md/json, benchmark report.md/json, CLI summary | known limitation / approximate |
| `NOT_APPLICABLE` | benchmark detail text, CLI `soc`/`pst` messages, npz `note` | quantity not defined for this symmetry |
| `UNCONVERGED` | per-state photon deficits (`w_deficit_*`, `n_deficit_gt_1e3`), CLI m_sp table | state/window not converged |
| `NEEDS_INPUT` | CLI BLOCKED messages, validate rows | missing input files / transform |
| `BLOCKED` | benchmark/validate reports | environment/data missing (e.g. official hr.dat) |

Pt Benchmark 7 final: **PASS=13, WARN=5, FAIL=0** (WARN includes
`pt_kramers_offgrid`, `pt_quantitative_status`, `pt_pst_gauge_invariance`,
`pt_hr_reconstruction`, `pt_band_interpolation_vs_wannier90`).

## 8. Automated tests (frozen regressions)

`tests/test_release_regressions.py` freezes: version, Pt fake-SOC output,
Kramers definition, non-Hermitian alias vs Hermitian WS interpolation,
bad .spn quality, TRS-violating operator, Pt PST gauge metadata, photon
deficit fields. Tolerances: Hermiticity 1e-12; Kramers mesh 1e-6 eV; gauge
spread >1e-3 is real; photon convergence |deficit| <1e-3.

Full-suite result of this release: see `RELEASE_TEST_RESULTS.txt`
(generated by the final `unittest discover` run).

## 9. Minimal reproducible example

```bash
python floquet.py --version
python floquet.py validate     --config configs/pbnv.json --outdir outputs/PbNV
python floquet.py spin_static  --config configs/pbnv.json --outdir outputs/PbNV
python floquet.py pst          --config configs/pbnv.json --outdir outputs/PbNV
python floquet.py soc          --config configs/pt.json   --outdir outputs/Pt
python floquet.py benchmark --which pt --outdir outputs/benchmarks
python -m unittest discover -s tests
```

## 10. Final project directory (top level of the deliverable)

```
floquet_soc/                 framework package (core, drives, solvers,
                             observables, workflows, benchmarks, cli.py)
floquet.py                   CLI entry point
configs/                     pbnv.json, pt.json, rashba*.json
benchmarks_data/Pt/          Pt tutorial dataset (no hr.dat)
tests/                       unittest suite incl. release regressions
outputs/PbNV/                PbNV benchmark outputs
outputs/Pt/                  Pt benchmark / interface outputs
outputs/benchmarks/pt/       Benchmark 7 report + figure
PT_*.md, STAGE_*.md          Stage / audit / acceptance reports
RELEASE_VALIDATION.md        this document
NEW_MATERIAL_GUIDE.md        input + validation checklist
VERSION                     0.1.0-rc1
```

## 11. Known unresolved limitations (not blocking the release)

1. Official Pt `wannier90_hr.dat`: missing (no wannier90 executable here; not
   in public data). Recipe in `PT_INTERPOLATION_REPORT.md`.
2. Pt reference-band comparison: WARN only (likely different Wannier settings).
3. PbNV PST: 3.8% gauge spread (approximate).
4. Shifted-mesh chk+eig interpolation: off-grid not validated.
5. m_sp convergence is material/drive dependent; always rerun the deficit test.
