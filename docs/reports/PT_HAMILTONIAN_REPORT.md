# PT Hamiltonian Report (updated after the WS/ndegen fix)

## Fixed: interpolation backend

Before (alias R=0..N-1, ndegen=1):
* off-grid `max|H(k)-H(k)^dag|` = **1.27 eV** (non-Hermitian),
* off-grid Kramers residual median 67 meV / max 0.89 eV,
* `spin_wannier_hermiticity` was masked by a forced 0.5(S+S^dag).

After (Wannier90 WS/ndegen: R=-N/2..N/2, ndegen = residue multiplicity):
* off-grid `max|H(k)-H(k)^dag|` = **1.81e-14 eV**,
* on-mesh reproduction 4.6e-14 eV,
* off-grid Kramers median 0.43 meV / max 4.54 meV,
* spin operator: mesh 1.44e-15, off-grid Hermiticity 6.7e-16 (no forced
  symmetrisation), R/ndegen convention identical to the Hamiltonian.

The WS/ndegen rule was validated on PbNV against the official
`wannier90_hr.dat`: R set 49/49 exact, ndegen 0 mismatches, interpolation
agreement 8.9e-6 eV (limited by the hr file's 1e-6 printing precision).

## Official Wannier90 hr.dat for Pt: still BLOCKED

* No `wannier90`/`wannier90.exe` on PATH or in searched install dirs.
* The tutorial repository path contains no hr.dat (404 for `Pt_hr.dat` and
  `wannier90_hr.dat` via jsDelivr/gh-proxy).
* `Pt.win` has `#restart = plot` commented and no `write_hr`;
  `Pt.nnkp/.chk/.eig/.win` are present, so `restart=plot`+`write_hr=.true.`
  is feasible given an executable (the Windows absence is not a data
  incompatibility).

Consequence: the official-hr comparison and the `Pt_band.dat` reference
comparison remain pending. Index-aligned deviation vs `Pt_band.dat` is
median 54.1 meV / max 1141.8 meV; this is reported as WARN, not PASS,
because the reference may have been produced with different Wannier settings
(e.g. `use_ws_distance`). Backend label:

**Grid-consistent WS/ndegen reconstruction (Hermitian); official-hr
comparison pending.**

## Spin-operator TRS check

`spin_wannier_trs` on the DFT mesh: residual 0.0000 (PASS). Off-grid TRS of
the *operator* in the Wannier gauge is not directly comparable without the
gauge transform, so no claim is made beyond the mesh check.
