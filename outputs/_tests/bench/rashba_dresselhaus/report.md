# Benchmark 3: Rashba+Dresselhaus, PST

Summary: **PASS**=5, **FAIL**=0, **WARN**=0, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `rd_splitting_analytic` | **PASS** | max |DeltaE - 2|d(k)|| = 2.78e-17 eV |
| `rd_alpha_eq_beta_collinear` | **PASS** | max |S - +-(1/2)(1,-1,0)/sqrt2| = 1.54e-15 |
| `rd_pst_collinearity_one` | **PASS** | M eigenvalues = [1. 0. 0.] -> (1,0,0) |
| `rd_pst_axis` | **PASS** | PST axis = [ 0.70710678 -0.70710678  0.        ], angle to (1,-1,0)/sqrt2 = 0.00e+00 deg |
| `rd_alpha_ne_beta_noncollinear` | **PASS** | collinearity(alpha!=beta) = 0.6500 < 1 |
