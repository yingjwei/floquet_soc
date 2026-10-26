# Benchmark 5: Rashba ferromagnet (exchange vs SOC)

Summary: **PASS**=4, **FAIL**=0, **WARN**=0, **BLOCKED**=0, **SKIP**=0

| check | status | detail |
|---|---|---|
| `exchange_splitting_at_Gamma` | **PASS** | Delta E(Gamma) = 0.400000 eV, 2M = 0.400000 eV (tol 1e-12) |
| `Gamma_spin_is_out_of_plane` | **PASS** | |<S_z>| = 0.500000, max|S_xy| = 0.00e+00 |
| `finite_k_splitting_has_soc` | **PASS** | Delta E(k=0.1) = 0.400500 eV differs from 2M by 5.00e-04 eV (SOC contribution) |
| `gamma0_static_limit` | **PASS** | max folded |eps_F - E_static| at k=0.1 = 2.22e-16 eV (tol 1e-10, modulo hw=0.4) |
