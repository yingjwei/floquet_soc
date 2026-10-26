# Benchmark 2: Rashba SOC model

Summary: **PASS**=6, **FAIL**=0, **WARN**=0, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `rashba_splitting_analytic` | **PASS** | max |DeltaE - 2 alpha k| = 6.94e-18 eV (tol 1e-10) |
| `rashba_spin_momentum_locking` | **PASS** | max |S - S_analytic| = 4.37e-16 (tol 1e-10) |
| `rashba_Sz_zero` | **PASS** | max |S_z| = 4.16e-16 |
| `rashba_spin_norm_half` | **PASS** | max ||S| - 1/2| = 2.78e-16 (hbar units) |
| `rashba_pst_metric_inplane` | **PASS** | M eigenvalues = [ 0.5  0.5 -0. ], expected (0.5,0.5,0) (max dev 2.2e-16) |
| `rashba_hermiticity` | **PASS** | max|H-H^dag| = 0.0e+00 |
