# Benchmark 7: fcc Pt (spinor SOC, chk+eig input)

Summary: **PASS**=13, **FAIL**=0, **WARN**=5, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `pt_model_source` | **PASS** | num_wann=18, num_bands=24, kpoints=64, mp_grid=[4, 4, 4]; source=wannier90.chk+eig Fourier interpolation (+win metadata); no hr.dat |
| `pt_hr_reconstruction` | **WARN** | no official wannier90_hr.dat available: H(R) rebuilt from chk+eig with the Wannier90 WS/ndegen convention (wannier90_ws); Hermitian off-grid and validated against the official PbNV hr.dat, but Pt off-grid accuracy is limited by the 4x4x4 mesh. Backend label: grid-consistent WS/ndegen reconstruction; official-hr comparison pending. |
| `checkpoint_transformation` | **PASS** | U unitarity 1.35e-13, v isometry 1.35e-13, lattice error 9.10e-18 |
| `checkpoint_hamiltonian_consistency` | **PASS** | max|H_W-H_hr| = 4.62e-14 eV (tol 0.001 eV over 64 k) |
| `spin_wannier_hermiticity` | **PASS** | max|S-S^dag| = 0.00e+00 |
| `spin_wannier_hermiticity_offgrid` | **PASS** | WS/ndegen interpolation: max|S(k)-S(k)^dag| off-grid = 6.66e-16 |
| `spin_wannier_interpolation` | **PASS** | max|S(k_mesh)-S_W(k_mesh)| = 1.44e-15 |
| `spin_wannier_trs` | **PASS** | max|sorted eig S(-k) + reverse sorted eig S(k)| = 0.0000 (0 expected for a TRS-preserving non-magnetic system; large values indicate mismatched runs or inaccurate spin matrices) |
| `pt_hermiticity_offgrid` | **PASS** | max|H(k)-H(k)^dag| on 200 off-grid k = 1.81e-14 eV (tol 1e-10 eV; Hermiticity must be checked before Kramers) |
| `pt_kramers_mesh` | **PASS** | delta_PT=max_n|E_{2n-1}-E_{2n}| on the original DFT mesh: max 7.934e-05 meV, median 7.532e-10 meV (tol 1e-6 eV = 1e-3 meV, from DFT eigenvalue precision 2e-10 eV) |
| `pt_kramers_offgrid` | **WARN** | off-grid Kramers residual: median 0.433 meV, max 4.54 meV at k=[-0.1027, -0.3645, 0.0031]. WARN, not PASS: for a strictly PT-symmetric static crystal the exact residual is zero; these meV values are numerical/interpolation error of the 4x4x4 WS/ndegen Hamiltonian, not physical. The error must be tracked against an official Wannier90 hr.dat before any PASS. |
| `pt_quantitative_status` | **WARN** | Pt is a software-interface and error-detection benchmark only at present: the official Wannier90 hr.dat is unavailable, the Pt_band.dat reference is not reproduced (index-aligned median 54 meV / max 1.14 eV), and the off-grid Kramers residual is a numerical artefact. Do NOT use Pt for quantitative SOC/Floquet material conclusions until the official-interpolation comparison passes. |
| `pt_soc_splitting_not_applicable` | **PASS** | PT-symmetric non-magnetic static Pt: SOC spin splitting is NOT APPLICABLE (Kramers pairs must be degenerate). `soc` reports delta_PT and inter-doublet gaps, no spin-paired splitting. |
| `pt_pst_gauge_invariance` | **WARN** | NOT_APPLICABLE: single-band PST is gauge dependent for a Kramers-degenerate PT system; the random-U(2) test shows the PST axis is arbitrary (dot product down to 0.21). |
| `pt_spn_quality` | **PASS** | input .spn (QE pw2wannier90): eig range [-1.0000,+1.0000], sum-rule [0.333,3.000] (Pauli bound 1, sum 3) |
| `pt_spin_texture_finite` | **PASS** | 729 k-points x 18 bands; figure band 9 |
| `pt_band_polarization_bound` | **PASS** | max per-band |<sigma>| = 1.000 (Pauli bound 1) |
| `pt_band_interpolation_vs_wannier90` | **WARN** | index-aligned deviation over 593 path points: median 54.1 meV, max 1141.8 meV. WARN, not PASS: the official Pt hr.dat is unavailable and Pt_band.dat may have been produced with different Wannier settings (use_ws_distance); this comparison cannot validate the off-grid Hamiltonian. |
