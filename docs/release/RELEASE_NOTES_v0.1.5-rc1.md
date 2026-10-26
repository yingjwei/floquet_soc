# Release notes — floquet_soc v0.1.5-rc1

**Type:** numerical-core patch (no API break; no physics-model change).
**Parent:** v0.1.4-rc1. **Patch target:** the Wannier90 WS interpolation set and
its use by the Wannier-basis spin operator.

> Naming note: `v0.1.1-rc1` … `v0.1.4-rc1` were already frozen earlier in the
> project, so this patch is numbered **v0.1.5-rc1** ("the WS/interpolation
> patch") rather than re-using v0.1.1.

## Why

The MoSSe SOC Wannier90 run uses an **odd 15×15×1** k mesh.  Two related
defects appeared when validating it:

1. `core/wannier.py::ws_supercell` built the R set with
   `np.arange(-N // 2, N // 2 + 1)`.  For odd N, Python floor division gives
   `-N//2 == -(N+1)//2`, so N=15 produced **16** values (−8…7) instead of 15
   (−7…7), and N=1 produced 2 values instead of 1.  Result: for MoSSe the
   helper returned 512 R vectors with `ndegen ∈ {2,4,8}` (16×16×2) instead of
   the correct 225.  `H(k)` was unaffected (it uses the R list from
   `wannier90_hr.dat`), but the Wannier-basis spin operator used the wrong
   set and its **off-grid interpolation became non-Hermitian**
   (`max|S−S†| = 0.999`), which `validate` correctly flagged as FAIL.

2. Even with the R-count fixed, a reconstructed box set is not the Wannier90
   WS set for non-orthogonal/hexagonal meshes: MoSSe's actual WS set
   (`wannier90_hr.dat`, 225 vectors) reaches |R| ≈ 10, i.e. it is the
   hexagonal Wigner–Seitz cell, not the box −7…7.  The spin operator must
   interpolate with **the same R/ndegen set as H(k)**.

## Changes

### `floquet_soc/core/wannier.py`
- **Fix** `ws_supercell`: ranges are now `-(N//2) … N//2` (15 values for
  odd N=15, 17 for even N=16); `ndegen` folding is unchanged and remains
  correct for both parities.
- **New** `read_wsvec(path)`: parses `wannier90_wsvec.dat`
  (`use_ws_distance=.true.`): unique R set, block count, `n_wsvec`, and
  optionally the listed supercell vectors.  Used for consistency tests.

### `floquet_soc/core/spin_wannier.py`
- `WannierSpinOperator(..., interpolation_set=(R, ndegen))`: the operator can
  now use the Wannier90 R set directly.  Validation accepts even-N WS sets,
  where `len(R) > nk` (e.g. PbNV 6×6 mesh → 49 R with `ndegen = 1/2/4`), with
  the conditions `R` covers the mesh residue classes one-to-one and
  `Σ_R 1/ndegen = nk`.
- `build_wannier_spin_operator` passes `model.Rlist/model.ndegen` when the
  model carries an `hr.dat`-based R set, so **spin and H use identical
  interpolation data**.  `R_source` and a new check
  `spin_wannier_interpolation_set` expose this in `validate`.
- Fallback to the reconstructed WS box (unshifted mesh) / alias grid
  (shifted mesh) is retained when no Wannier90 R set is available.

### `floquet_soc/core/cache.py`
- **Fix** the parsed-model disk cache key: analytic-model configs have no data
  files, so hashing only config-dir + file stats gave the same cache key to
  `rashba.json` and `rashba_dresselhaus.json` in the same directory; a stale
  model could be loaded (caught by `test_pst_analytic`, which got
  Rashba instead of Rashba–Dresselhaus and returned collinearity 0.5 instead
  of 1.0).  The key now also includes `kind`, `name` and the serialized
  `model` block.

### `tests/test_ws_grid.py` (new, 8 tests)
- odd/even/single-point `ws_supercell` shapes, symmetry under R→−R, ndegen;
- MoSSe: `wannier90_wsvec.dat` R set == `hr.dat` R set, block count =
  `nrpts × nw²`, all wsvec vectors are supercell multiples of `mp_grid`;
- MoSSe: spin operator uses `hr.dat` R set; **off-grid Hermiticity
  `1.11e-15`** and mesh round-trip `2.42e-11`;
- MoSSe: `from_chk_eig` reproduces `chk+hr` on the DFT mesh (`<1e-3 eV`),
  regression for the odd-grid `ws_supercell` fix;
- PbNV (even 6×6 mesh): `len(R) = 49 > nk = 36`, `ndegen` 1/2/4 accepted,
  off-grid Hermiticity machine precision (guards the new validation rule).

## Validation

`validate` with unchanged data:

| material | before (v0.1.4 core) | after |
|---|---|---|
| MoSSe (15×15×1) | — (`off-grid spin` FAIL, 25 PASS/1 FAIL/2 WARN with v0.1.0 core) | **28 PASS / 0 FAIL / 1 WARN** |
| PbNV | 26 PASS / 2 WARN | **27 PASS / 0 FAIL / 2 WARN** |
| Pt | 24 PASS / 4 WARN | **25 PASS / 0 FAIL / 4 WARN** |

MoSSe off-grid spin Hermiticity: `9.99e-1 → 1.11e-15`; spin R set source now
`wannier90_hr.dat` (225 R, ndegen all 1).  PbNV/Pt `H(k)` and on-grid spin
values are unchanged; only the off-grid spin interpolation set changed (and,
for even grids, it is now the same WS set as H, which slightly changes
off-grid interpolation but not grid values).

Full test suite (85 tests): see `RELEASE_TEST_RESULTS_v0.1.5-rc1.txt`;
`test_version_string` now compares `__version__` with `VERSION` instead
of a hard-coded old release string.

## Compatibility / limitations

- `H(k)` path is unchanged; on-grid spin operators are unchanged.
- `from_chk_eig` still uses a reconstructed WS box set (no wsvec input in its
  API); for odd meshes this is now sized correctly, and the MoSSe test shows
  it reproduces `chk+hr` on the mesh.  Off-grid `from_chk_eig` interpolation
  for hexagonal meshes remains an approximation (documented in its docstring).
- Deep disentangled Wannier bands are still not physical; the MoSSe report
  keeps them as WARN and uses only the frozen gap-edge manifold.
