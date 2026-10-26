# Final acceptance status (four levels)

## Level 1 - mathematics and software implementation: PASS

* WS/ndegen reconstruction Hermitian to 1.8e-14 eV off grid.
* Spin operator: mesh 1.4e-15, off-grid Hermiticity 6.7e-16, no forced
  symmetrisation.
* New physics tests 10/10 OK; affected existing modules 17 tests OK (2 skips);
  PbNV validate 26 PASS / 2 WARN / 0 FAIL; Pt validate 24 PASS / 4 WARN
  （2026-10-23 冻结复跑：当前代码输出 24/4；旧记录的 23/4 来自 phase-2 之前
  生成的 validate 文件，冻结快照中旧文件保留为
  `outputs/Pt/validate_pre_freeze_2026-10-23.*`）。

## Level 2 - static physics symmetry validation: PARTIAL (Pt WARN, PbNV PASS)

* PbNV: T-Kramers at TRIM only, residuals 0.31-0.49 meV (Wannier truncation),
  real SOC pair splitting away from TRIM; PST axis stable (0.998) but
  collinearity spread 0.038 -> **WARN (approximate Ising-like)**.
* Pt: DFT mesh Kramers exact (2e-10 eV), model mesh 4.6e-14 eV (PASS);
  off-grid residual median 0.43 meV / max 4.54 meV is numerical
  interpolation error -> **pt_kramers_offgrid WARN**. SOC splitting and PST
  are **NOT_APPLICABLE** for this PT-symmetric system.

## Level 3 - Floquet numerical convergence: PASS

* n_step: dEbar=1.8e-6 eV, dS=1.7e-5 at 120 steps, decreasing.
* m_sp: per-state photon deficit max 1.0 at m_sp=16 -> 2 states
  unconverged; m_sp=32/48 deficit 3.9e-14, 0 states >1e-3 ->
  accepted m_sp >= 32; production [160, 32] PASS.

## Level 4 - real-material quantitative production accuracy: NOT YET

* Pt: official Wannier90 hr.dat unavailable; `Pt_band.dat` not reproduced
  (index-aligned median 54 meV, max 1.14 eV); off-grid Kramers residual is
  numerical. Pt is a software-interface / error-detection benchmark only.
* PbNV: official hr.dat present and validated; PST approximate (WARN).
  PbNV can be used as an approximate 2D Ising-like reference, not for
  precision SOC constants.

## Pt Benchmark 7 status

PASS=13, WARN=5, FAIL=0. New WARN rows: `pt_kramers_offgrid`
(numerical residual, not physical) and `pt_quantitative_status`
(interface benchmark only). No PASS was manufactured by loosening physics
tolerances.

## Blocked / external step

Generate the official Pt hr.dat in a Linux/HPC environment with
`restart = plot` + `write_hr = .true.` (recipe in `PT_INTERPOLATION_REPORT.md`),
then re-run: Hamiltonian comparison, Kramers on/off grid, SOC classification,
PST gauge test. Only then can Pt be promoted beyond an interface benchmark.
