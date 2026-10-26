# PbNV PST credibility statement (WARN, kept as approximate)

Data: `outputs/PbNV/pst.npz`, bands 110-113, 19x17=323 k.

| quantity | value |
|---|---|
| M eigenvalues | (0.9570, 0.0401, 0.0029) |
| collinearity | **0.9570** |
| PST axis | (0.009, 0.0104, 0.9999) ~ z |
| mean polarization | 0.797 |
| gauge test: axis dot min | 0.998 |
| gauge test: collinearity spread | **0.038** |
| bootstrap over k: collinearity | 0.9561 (range 0.9421-0.9730) |

Per-band collinearity (band alone): {110: 0.95, 111: 0.95, 112: 0.965, 113: 0.964}.
Region dependence: {'far |k|>=0.15': 0.964, 'near |k|<0.15': 0.865, 'splitting>=5meV': 0.994, 'splitting<5meV': 0.842}.

Effective k region (|S|>0.5 fraction, |Sz|>0.8 fraction):

| band | mean |S| | frac |S|>0.5 | frac |Sz|>0.8 |
|---|---|---|---|
| 110 | 0.973 | 100.0% | 94.1% |
| 111 | 0.973 | 100.0% | 94.1% |
| 112 | 0.529 | 53.3% | 25.1% |
| 113 | 0.629 | 66.9% | 36.8% |

Interpretation: PbNV is non-centrosymmetric; TRS alone only enforces
Kramers degeneracy at TRIM points (Gamma, X, S, Y), where the official hr
model gives 110-111 residuals 0.309/0.386/0.412/0.487 meV (and 0.11-0.32 meV
for 112-113). Away from TRIM the pair splitting is real (up to ~30 meV on the
path, ~61 meV on the 19x17 mesh). The PST axis is stable (dot 0.998) but the
collinearity has a 3.8% gauge spread, so the result is **WARN (approximate
Ising-like)**, not a precision PST constant. Effective region: the
out-of-plane polarization is strong across the whole mesh (table above); the
gauge ambiguity is concentrated near TRIM (small splitting).

Contrast with Pt: Pt has P+T, so Kramers degeneracy holds at **every k**
(PT-Kramers), making single-band spins/PST gauge-dependent everywhere; PbNV
has T only, so the degeneracy is restricted to TRIM and away from TRIM the
single-band splitting is physical. They must not be treated with the same
PST rules.
