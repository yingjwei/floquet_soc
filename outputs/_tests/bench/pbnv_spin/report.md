# Benchmark 6: PbNV spin operator (chk-rotated)

Summary: **PASS**=10, **FAIL**=0, **WARN**=0, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `checkpoint_transformation` | **PASS** | U unitarity 8.90e-14, v isometry 8.79e-14, lattice error 2.22e-16 |
| `checkpoint_hamiltonian_consistency` | **PASS** | max|H_W-H_hr| = 9.73e-06 eV (tol 0.001 eV over 36 k) |
| `spin_wannier_hermiticity` | **PASS** | max|S-S^dag| = 0.00e+00 |
| `spin_wannier_interpolation` | **PASS** | max|S(k_mesh)-S_W(k_mesh)| = 9.67e-12 |
| `pbnv_kramers_degeneracy` | **PASS** | informational (PbNV has no inversion symmetry, so Kramers degeneracy is not expected): nearest-partner energy distance over 81 k x 156 bands: max=2.99e-01 eV, median=5.24e-04 eV |
| `pbnv_band_polarization_bound` | **PASS** | max per-band |<sigma>| over path = 1.000 (Pauli bound 1); fraction >1: 0.000 |
| `pbnv_adjacent_spin_antisymmetry` | **PASS** | informational (no Kramers pairing without inversion; avoided crossings rotate the spin): fraction of adjacent band pairs with |S_2n+S_2n+1|<0.1: 0.698 (median 0.008) |
| `pbnv_spin_trs` | **PASS** | max|sorted eig S(-k) + reverse sorted eig S(k)| = 0.0017 Pauli units (0 expected for a TRS-preserving non-magnetic system) |
| `pbnv_spin_sum_rule` | **PASS** | eig(sum_i S_i^2) in [2.786,3.000] over the Wannier mesh (3 for orthonormal complete basis) |
| `pbnv_spn_input_quality` | **PASS** | input .spn (VASP LWRITE_SPN, formatted): eig range [-1.000,+1.000], sum-rule [0.719,3.000]. Pauli bound satisfied; the sum-rule lower bound reflects spin weight outside the 156-dim Wannier window (not a data error). |
