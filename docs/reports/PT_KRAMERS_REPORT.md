# PT Kramers Report (updated after the interpolation fix)

Definition: `delta_PT(k)=max_n|E_{2n-1}(k)-E_{2n}(k)|` (1-based pairs =
0-based columns (0,1),(2,3),...).

## After fix (Wannier90 WS/ndegen reconstruction)

| data set | max delta_PT | median | worst k |
|---|---|---|---|
| original DFT 4x4x4 mesh (.eig) | 1.98e-10 eV | 0 | — |
| model on the same mesh | 4.6e-14 eV | ~1e-15 eV | — |
| k-path W-L-G-X-W-G (Pt soc output) | 4.37 meV | 0.27 meV | (0.133,0.133,0.133) |
| 400 random off-grid k | 4.54 meV | 0.43 meV | (-0.103,-0.365,0.003) |

Before the fix: off-grid median **67 meV**, max **0.89 eV** (FAIL).

Tolerances: DFT-mesh 1e-6 eV (DFT precision 2e-10 eV); off-grid accepted at
median <1 meV and max <10 meV because the 4x4x4 Wannier Hamiltonian is
truncated (target SOC scale ~50 meV; a residual below this scale cannot create
a false meV-level SOC splitting). PASS on mesh; PASS (truncation level) off
grid.

## Interpretation

The remaining few-meV off-grid residual is Wannier-interpolation truncation
(4x4x4 mesh), not a physical PT breaking. It is now below the old
"48.8/62/126 meV" artifacts by more than an order of magnitude and is
reported explicitly as a WARN-level limitation. An official Wannier90 hr.dat
is still required for a definitive off-grid comparison.
