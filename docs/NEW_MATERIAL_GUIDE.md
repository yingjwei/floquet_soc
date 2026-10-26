# NEW MATERIAL GUIDE — floquet_soc v0.1.0-rc1

Scope: adding a new material as an independent research task; the generic
core must not be re-developed. Follow this checklist before quoting any
physical conclusion.

## 1. Required Wannier90 input files

| file | needed for | notes |
|---|---|---|
| `wannier90.win` | lattice, num_wann, mp_grid, bands, kpoints | kept with the run |
| `wannier90_hr.dat` | Hamiltonian | **preferred**; official interpolation |
| `wannier90.chk` | U matrices, Wannier centres, k mesh | required for spin rotation |
| `wannier90.eig` | DFT eigenvalues on the mesh | required with chk |
| `wannier90.spn` | spin matrices | VASP LWRITE_SPN or QE pw2wannier90; must be in the Wannier band basis or rotatable with chk |
| `wannier90.nnkp` / `.wout` | provenance / checks | optional but recommended |
| `wannier90_centres.xyz` | centres | optional |

Without `hr.dat`, `chk+eig` is used as a **grid-consistent fallback**:
* unshifted mesh -> WS/ndegen reconstruction (Hermitian, PT residual ~meV);
* shifted mesh  -> alias fallback: grid values exact, off-grid physics
  **NOT VALIDATED** (spin check token `NOT_VALIDATED`).

## 2. Static Hamiltonian checks (must pass)

1. `python floquet.py validate --config <cfg>` — 0 FAIL. WARN rows must be read.
2. On-mesh reproduction: `checkpoint_hamiltonian_consistency` < 1e-3 eV.
3. Off-grid Hermiticity: `H(k)-H(k)^dag` < 1e-10 eV (WS backends).
4. Off-grid Kramers (PT systems): report δ_PT; exact zero is expected
   physically — if meV-level residual, mark as numerical and do not claim
   degeneracy PASS.
5. Static bands vs a DFT/Wannier90 reference (if available): set an explicit
   tolerance based on the Wannier interpolation quality; do not hide outliers
   with nearest-band metrics.

## 3. Symmetry classification (required in the config)

```json
"system": {"type": "pt_nonmagnetic", "inversion": true, "time_reversal": true}
```
Allowed values: `pt_nonmagnetic` (SOC splitting NOT_APPLICABLE; report δ_PT
and inter-doublet gaps), `noncentrosymmetric_soc` (real SOC spin splitting
allowed; use spin-paired analysis with gauge test), `magnetic` (separate
exchange vs SOC; effective Zeeman under drive), `pst_candidate` (PST metrics
only if gauge-invariant).

## 4. Drive configuration

```json
"drive": {"polarization": "circular", "sign": -1, "hw": 0.4, "gamma": 0.1}
```
* `hw` in eV, `gamma` broadening; `a_ref` from lattice unless given.
* `numerics`: `n_step` and `m_sp`; run `converge` for the actual material.
* Photon truncation: require per-state deficit `|1 - sum_m w| < 1e-3` for all
  states used in spectra/PST; `m_sp=32` is only a tested choice, not a rule.

## 5. SOC analysis choice

* PT non-magnetic: run `soc`; report Kramers residual + inter-doublet gaps;
  do not call any gap a SOC spin splitting.
* Non-centrosymmetric: spin-paired splitting is meaningful only when the pair
  is non-degenerate; use the gauge test and a splitting threshold.
* Linear Rashba fit: when the gap-edge pair splitting is linear in |k| around a
  high-symmetry anchor (e.g. Γ, Y), run
  `python floquet.py rashba --config <cfg>`; the config section `rashba_fit`
  sets `anchor`, `direction` (kpath label or vector), `nk`, `kmax`, `kmin`,
  `step`, `r2_thr`, `alpha_ratio_thr`, `split_r2_thr`. It fits
  `E(k) = a·k² + E0 ± α·|k|` on both branches, requires opposite-sign α and a
  linear ΔE=2α_R|k| (split_r2_thr), and reports α_R (eV·Å) plus
  LINEAR/NOT_LINEAR. A high branch R² alone is NOT sufficient (two degenerate
  branches can be fitted with a fake same-sign α).
* Magnetic: separate exchange splitting, SOC contribution, and drive-induced
  effective Zeeman; do not call photo-induced splittings "SOC enhancement".
* PST: require the gauge test (random U(2) in near-degenerate subspaces);
  axis dot <0.99 or collinearity spread >1e-2 -> WARN/NOT_APPLICABLE.

## 6. Floquet checks before physical conclusions

1. `n_step` convergence: dEbar < tol_ebar, dS < tol_spin.
2. `m_sp` photon-weight convergence: per-state deficits < 1e-3,
   `accepted_m_sp` recorded; otherwise mark affected states UNCONVERGED.
3. Static limits: γ→0 reproduces static bands/spin texture.
4. Matching: Floquet↔static overlap purity saved (`purity`); never compare
   by quasienergy ordering.

## 7. Output metadata to keep

* Figure policy (framework-wide): **no top titles/suptitles**; one unified type
  scale in every figure (axis/colorbar label 11, ticks 10, legend 9, gamma-point
  labels 10); colorbar labels use the axis-label size. Do not override these
  per figure (`floquet_soc/plot_style.py::STYLE/FONTS`).

* npz files: kpts, energies, spins, `purity`, `system_type`, `note`,
  gauge fields, convergence fields; JSON/MD reports: WARN/NOT_APPLICABLE/
  UNCONVERGED/BLOCKED tokens.
* README for the material with file index and reproducibility commands.

## 8. Config template

See `configs/pbnv.json` (2D non-centrosymmetric) and `configs/pt.json`
(PT non-magnetic, chk+eig backend). Start from the closest one and change
only material-specific fields.
