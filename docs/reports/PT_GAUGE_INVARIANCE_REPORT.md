# PT Gauge-Invariance Report (updated)

## Test on the true DFT Kramers subspaces

Old test rotated eigenvectors inside *spurious* non-degenerate pairs created
by the broken interpolation; that was not a genuine gauge transformation.
The new test uses the original DFT mesh, where the Kramers pairs are exactly
degenerate (delta_PT max 1e-10 eV): 64 k x 9 pairs = 576 pair rotations, 5
random U(2) seeds each.

| quantity | base | random-U(2) seeds |
|---|---|---|
| PST collinearity | 0.3965 | 0.3678 - 0.4131 |
| PST axis dot with base | 1 | min **0.345** |
| gauge-invariant subspace spin `0.5 Tr(P_pair sigma)` | — | mean norm **0.0024**, max 0.335 |

The single-band PST axis is not gauge invariant: rotating inside the exact
Kramers subspaces changes the axis completely. The gauge-invariant subspace
spin is essentially zero on average (0.0024), so there is no persistent spin
texture in the degenerate subspace either.

Verdict: **FAIL for single-band PST / NOT_APPLICABLE for Pt as a PST
material.** Pt remains a good benchmark for spin operators, Hermiticity and
degenerate-subspace handling. The PST command now returns an explicit
NOT_APPLICABLE status and stores the gauge-test data
(`outputs/Pt/pst.npz: gauge_collinearity, gauge_axis, gauge_invariant_spin`).

## PbNV (non-centrosymmetric)

Near-degenerate subspaces (splitting <5 meV, mainly near TRIM) were rotated
randomly: collinearity spread 0.038, PST axis dot min 0.998, invariant
subspace spin mean 0.035. The PbNV PST axis is stable (dot 0.998) and the
collinearity remains 0.957; therefore the PbNV result is kept as an
approximate physical result, with the near-degeneracy caveat.
