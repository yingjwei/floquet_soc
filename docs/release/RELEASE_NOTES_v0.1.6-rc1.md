# Release notes — floquet_soc v0.1.6-rc1

**Type:** numerical-core patch (Wannier90 interpolation convention). No API
break, no changes to the physics formulas.
**Parent:** v0.1.5-rc1. **Trigger:** MoSSe wsvec interpolation audit.

## Why (audit result: v0.1.5 and earlier were incomplete — FAIL)

MoSSe's Wannier90 run has `use_ws_distance=.true.` and **no**
`write_ndegen_applied` token in `wannier90_wsvec.dat`. Therefore the R list in
`wannier90_hr.dat` is *folded*, and the pair-dependent supercell shifts listed
in `wannier90_wsvec.dat` must be applied before the plain Fourier transform
(Wannier90 `ws_apply_ndegen` + plain Fourier on the expanded list). The
previous interpolation used the folded R list only.

Consequences measured with an independent implementation
(`materials/MoSSe/scripts/wsvec_audit.py`, see `MOSSE_WSVEC_AUDIT.md`):

| quantity | plain (v0.1.5) | Wannier90 ws-aware | difference |
|---|---|---|---|
| on-grid H / S | — | — | 2.3e−13 / 5.6e−15 (identical) |
| off-grid max \\|ΔE\\| near Γ→K/K/K′ | — | — | up to **0.66 meV** |
| Γ→K Rashba α_R | 64.84 meV·Å | **44.28 meV·Å** | **−31.7%** |
| Γ→M Rashba α_R | 64.06 | **53.73** | −16% |
| K-point Floquet γ=0.02 \\|ΔS̄\\| | 5.11e−5 | 1.54e−4 | scheme difference 1.61e−4 (3.1× the signal) |

wsvec statistics: 260100 blocks = 225 R × 34² pairs; 11,626 blocks with
non-zero supercell shift; 506 blocks with 2–3 equidistant shifts; expanded
R list 225 → **261**.

## Changes

### `core/wannier.py`
- `read_wsvec`: now returns the full block structure (`Rblocks`, `ijs`,
  `n_wsvec`, `shift_start`, `shifts_flat`) plus `header`,
  `use_ws_distance`, `write_ndegen_applied`.
- `WannierModel(..., wsvec_path=...)` / `from_files` / `from_chk`: build the
  expanded representation (`Rlist_full`, `HR_full`) via Wannier90
  `ws_apply_ndegen` weights `1/(ndegen·ws_ndeg)`; `write_ndegen_applied=true`
  is honored (no double counting); `Hk/Hk_batch` use the expanded list.
- `displacement_tensor` uses the expanded R list (same shortest Bravais
  vectors as H(k)); `__setstate__` keeps old pickled models safe.

### `core/spin_wannier.py`
- `WannierSpinOperator` applies the same expansion to the Wannier-basis spin
  operator; `R_source="wannier90_wsvec.dat"`; off-grid Hermiticity 1.1e−15.

### `drives/peierls.py`
- `peierls_Hk` and `peierls_fourier_components` use the expanded R list and
  `d = (R+shift)_cart + τ_j − τ_i`, so the Peierls phase and the
  interpolation see the same hopping vectors. A=0 still reproduces `Hk`
  exactly (0.0).

### `cli.py`
- `build_model` passes `files["wsvec"]` into the model.

### `core/cache.py`
- Cache key includes the framework version and the analytic-model block
  (prevents stale unpickled models across patches).

### tests
- New `tests/test_wsvec.py` (5 tests): flags/statistics, expanded model,
  on-grid equivalence + off-grid difference, Peierls vector consistency,
  spin-operator expansion, Rashba-convention regression (44.3 vs 64.8).

## Validation

- patched framework vs independent ws-aware implementation:
  H(k) diff **0.0**, spin **1.4e−15**, Peierls H(k,A) **0.0**;
- MoSSe `validate`: **28 PASS / 0 FAIL / 1 WARN** (spin R source
  `wannier90_wsvec.dat`, 261 R);
- full suite: see `RELEASE_TEST_RESULTS_v0.1.6-rc1.txt`;
- PbNV / Pt have no `_wsvec.dat` → `HR_full=None`, behavior unchanged
  (regression tests pass).

## Compatibility

- On-grid results (all DFT-mesh quantities, K/K′ splittings, k-mesh
  comparisons) are unchanged to machine precision.
- Off-grid quantities **do change** where they should: Rashba α_R, off-grid
  band energies/spin texture, Peierls-driven Hamiltonians and Floquet
  observables. Any v0.1.5-or-earlier off-grid number must be recomputed.
- `from_chk_eig` still uses a reconstructed folded Fourier representation
  (no wsvec input); documented approximation, not used for MoSSe.
