# Stage-B benchmark: driven Rashba model

Summary: **PASS**=8, **FAIL**=0, **WARN**=0, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `B1_zero_drive_quasienergy` | **PASS** | max|eps - E_static mod hw| = 6.66e-16 eV (tol 1e-10) |
| `B1_zero_drive_spin_texture` | **PASS** | max|S_bar - S_static| = 1.15e-13 (tol 1e-10) |
| `B2_photon_weight_normalization` | **PASS** | max|sum_m w - 1| = 3.75e-14 (m_sp=12) |
| `B3_sambe_vs_evolution` | **PASS** | max quasienergy mismatch Sambe vs evolution = 8.18e-12 eV (n_harm=3, tol 1e-5) |
| `B4_effective_convergence` | **PASS** | max|eps - eig(H_eff)| for w=1,2,4,8 eV: [1.27e-08, 3.16e-09, 7.91e-10, 1.98e-10]; decreasing=True, ratio(8/1)=1.56e-02 (tol: final <1e-4 eV) |
| `B5_spin_weight_identity` | **PASS** | max|sum_m s_am,i - S_bar| = 6.11e-16 |
| `B5_spin_spectrum_integrates_to_spin` | **PASS** | max|∫A^S dE - sum_alpha S_bar| = 6.20e-16 (sigma=2 meV, grid +-8 eV) |
| `B6_micromotion_reported` | **PASS** | max|S_bar - S_t0| = 8.948e-03 (difference between the correct periodic average and the t=0 approximation) |
