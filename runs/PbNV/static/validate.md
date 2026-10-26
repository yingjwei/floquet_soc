# validate: PbNV

Summary: **PASS**=27, **FAIL**=0, **WARN**=2, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `model_parse` | **PASS** | kind=wannier90, num_wann=156 |
| `win_metadata` | **PASS** | num_wann=156, num_bands=192, spinors=True, kpoints=36, mp_grid=[6, 6, 1]; source=wannier90.win |
| `H_R_hermiticity` | **PASS** | max|H(-R)-H(R)^dag| = 0.00e+00 (tol 1e-10) |
| `Hk_hermiticity` | **PASS** | max|H(k)-H(k)^dag| = 8.88e-15 (tol 1e-10, 5 k probes incl. off-mesh, wsvec=no) |
| `ndegen_consistency` | **PASS** | ndegen in [1,4], mismatched R/-R pairs: 0 |
| `wannier_centres` | **PASS** | tau source = wannier90_centres.xyz |
| `static_gap` | **PASS** | indirect gap on chosen k-path = 1.487060 eV (num_electrons=112) |
| `spn_format` | **PASS** | wannier90.spn: formatted-spn, header='Written by VASP:    unknown system', num_bands=192, num_kpts=36 |
| `spn_hermiticity` | **PASS** | max|S-S^dag| = 4.44e-16 |
| `spn_real_diagonal` | **PASS** | max|Im diag| = 2.22e-16 |
| `spn_pauli_eigenvalue_bound` | **PASS** | eigenvalues in [-1.000,+1.000]; Pauli bound [-1,1] |
| `spn_spin_sum_rule` | **WARN** | eig(sum_i S_i^2) in [0.719,3.000] (3 for a spin-complete basis; <3 means spin weight outside the Wannier band window) |
| `spn_band_polarization_bound` | **PASS** | max per-band |<sigma>| = 1.000 |
| `spin_hermiticity` | **PASS** | max|S-S^dag| = 4.44e-16 (tolerance 1e-10) |
| `spin_real_diagonal` | **PASS** | max|Im diag| = 2.22e-16 (tolerance 1e-10) |
| `spin_eigenvalue_bound` | **PASS** | eigenvalues in [-1.000,+1.000], Pauli bound [-1,1] |
| `spin_diagonal_polarization_bound` | **PASS** | max |<sigma>_band| = 1.000 (Pauli bound 1) |
| `spin_wannier_interpolation_set` | **PASS** | R set source = wannier90_hr.dat (49 R vectors, ndegen in [1,4]) |
| `spin_wannier_hermiticity` | **PASS** | max|S-S^dag| = 0.00e+00 |
| `spin_wannier_trs` | **PASS** | max|sorted eig S(-k) + reverse sorted eig S(k)| = 0.0017 (0 expected for a TRS-preserving non-magnetic system; large values indicate mismatched runs or inaccurate spin matrices) |
| `spin_wannier_interpolation` | **PASS** | max|S(k_mesh)-S_W(k_mesh)| = 9.67e-12 |
| `spin_wannier_hermiticity_offgrid` | **PASS** | WS/ndegen interpolation: max|S(k)-S(k)^dag| off-grid = 6.66e-16 |
| `spin_wannier_eigenvalue_bound` | **PASS** | eig range [-1.000,+1.000] (Pauli bound [-1,1]) |
| `spin_wannier_sum_rule` | **WARN** | eig(sum_i S_i^2) in [2.786,3.000] (3 for a spin-complete basis; <3 means spin weight outside the Wannier window) |
| `checkpoint_hamiltonian_consistency` | **PASS** | max|H_W-H_hr| = 9.73e-06 eV (tol 0.001 eV over 36 k) |
| `checkpoint_transformation` | **PASS** | U unitarity 8.90e-14, v isometry 8.79e-14, lattice error 2.22e-16 |
| `spin_operator_wannier_basis` | **PASS** | chk rotation applied (wannier90.chk); max|H_W-H_hr| = 9.73e-06 eV |
| `spin_texture_wannier` | **PASS** | gauge-consistent spin operators available; static spin texture can be computed (see `spin_static`) |
| `drive_config` | **PASS** | drive parameters parsed |
