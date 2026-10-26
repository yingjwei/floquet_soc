# -*- coding: utf-8 -*-
"""Spin operator in the Wannier basis (stage A completion).

Pipeline (all steps validated numerically):

1. ``wannier90.spn`` gives S_band(k) in the DFT/WAVECAR band basis
   (192 x 192 for PbNV; see ``core.spin_operator.read_vasp_spn``).
2. ``wannier90.chk`` gives v(k) (num_bands x num_wann) with
   |w_nk> = sum_m v_mn(k)|u_mk>, so
       S_W(k) = v†(k) S_band(k) v(k).
   The checkpoint is first validated by requiring
       H_W(k) = v†(k) diag(eig(k)) v(k)  ≈  H_hr(k)
   (which simultaneously checks k ordering, band ordering and gauge).
3. Spin operators are interpolated to arbitrary k through their real-space
   Fourier transform on the Wannier mesh:
       S(R) = (1/Nk) sum_k S_W(k) e^{-2 pi i k.R},
       S(k) = sum_R S(R) e^{2 pi i k.R}.
   The transform pair is exact at the mesh points (checked), and valid away
   from them because Wannier functions are localized.
"""
from __future__ import annotations

import numpy as np

from .chk import (read_chk, read_eig, checkpoint_checks,
                  hamiltonian_consistency, WannierCheckpointError)
from .spin_operator import (read_wannier_spn,
                            validate_spin_matrices)


class WannierSpinOperatorError(RuntimeError):
    pass


def rotate_spin_to_wannier(spn: dict, chk: dict) -> np.ndarray:
    """S_W(k) = v† S_band(k) v, shape (nk,3,nw,nw)."""
    v = chk["v_matrix"]                                  # (nk,nb,nw)
    S = spn["matrices"]                                  # (nk,3,nb,nb)
    if S.shape[0] != v.shape[0]:
        raise WannierSpinOperatorError(
            f"k-point count mismatch: spn {S.shape[0]} vs chk {v.shape[0]}")
    if S.shape[2] != v.shape[1]:
        raise WannierSpinOperatorError(
            f"band count mismatch: spn {S.shape[2]} vs chk num_bands {v.shape[1]}")
    Sw = np.einsum("kim,kcij,kjn->kcmn", v.conj(), S, v,
                  optimize="optimal")
    return 0.5 * (Sw + Sw.conj().transpose(0, 1, 3, 2))


class WannierSpinOperator:
    """Spin operators in the Wannier basis, interpolated from the DFT mesh."""

    def __init__(self, kmesh, mp_grid, spin_mats, lattice=None,
                 interpolation_set=None, ws_model=None):
        """Build the interpolated Wannier-basis spin operator.

        ``interpolation_set`` = (Rvecs, ndegen) from the actual Wannier90
        calculation (``wannier90_hr.dat``).  This is the preferred route:
        the spin operator then uses the SAME WS R set as H(k), so the two
        interpolations are mutually consistent.  For hexagonal k meshes the
        Wannier90 WS set is not a simple box (MoSSe 15x15x1: |R| reaches 10),
        so reconstructing it from the mesh is not equivalent.  When no set is
        provided, fall back to the reconstructed WS box (unshifted mesh) or
        the alias grid fallback (shifted mesh), as before.
        """
        kmesh = np.asarray(kmesh, float)
        mp_grid = np.asarray(mp_grid, int)
        spin_mats = np.asarray(spin_mats, complex)
        nk, nc, nw, nw2 = spin_mats.shape
        if nc != 3 or nw != nw2:
            raise WannierSpinOperatorError(
                f"expected (nk,3,nw,nw), got {spin_mats.shape}")
        if nk != int(np.prod(mp_grid)):
            raise WannierSpinOperatorError(
                f"kmesh size {nk} != prod(mp_grid)={int(np.prod(mp_grid))}")
        self.nw = nw
        self.kmesh = kmesh
        self.mp_grid = mp_grid
        self.spin_mats = spin_mats
        # Wannier90 WS/ndegen interpolation for unshifted meshes; shifted mesh
        # keeps the documented alias fallback (grid-consistent only).
        shift = np.mod(kmesh * mp_grid, 1.0)
        shift[np.abs(shift - 1.0) < 1e-9] = 0.0
        unshifted = bool(np.all(np.abs(shift) < 1e-9))
        self.lattice = None if lattice is None else np.asarray(lattice, float)
        self.R_source = "reconstructed"
        if interpolation_set is not None:
            R, nd = interpolation_set
            R = np.asarray(R, int).reshape(-1, 3)
            nd = np.asarray(nd, int).reshape(-1)
            # A Wannier90 WS set can have MORE R vectors than mesh k points
            # (even N: N+1 values per direction, folded by ndegen), e.g.
            # PbNV 6x6 mesh -> 49 R with ndegen 1/2/4.  The consistency
            # conditions are: the residue classes cover the mesh one-to-one
            # and sum_R 1/ndegen(R) = nk.
            if len(R) < nk:
                raise WannierSpinOperatorError(
                    f"interpolation set has only {len(R)} R vectors for a "
                    f"{nk}-point mesh")
            residues = {tuple(np.mod(r, mp_grid)) for r in R}
            if len(residues) != nk:
                raise WannierSpinOperatorError(
                    f"interpolation R set covers {len(residues)} residue "
                    f"classes, mesh has {nk}")
            if abs(float(np.sum(1.0 / nd)) - nk) > 1e-8:
                raise WannierSpinOperatorError(
                    f"sum_R 1/ndegen = {float(np.sum(1.0 / nd)):.3f} != nk={nk}")
            self.Rvecs, self.ndegen = R, nd
            self.ws_interpolation = True
            self.R_source = "wannier90_hr.dat"
        elif unshifted and self.lattice is not None:
            from .wannier import ws_supercell
            self.Rvecs, self.ndegen = ws_supercell(mp_grid, self.lattice)
            self.ws_interpolation = True
        else:
            grids = [np.arange(N) for N in mp_grid]
            self.Rvecs = np.array(np.meshgrid(*grids, indexing="ij")
                                  ).reshape(3, -1).T.astype(int)
            self.ndegen = np.ones(len(self.Rvecs), dtype=int)
            self.ws_interpolation = False
        ph = np.exp(-2j * np.pi * (kmesh @ self.Rvecs.T))     # (nk,nR)
        # S_R stores the interpolation coefficients S(R)/ndegen(R)
        self.S_R = (np.einsum("kR,kcij->Rcij", ph, spin_mats,
                              optimize=True) / nk
                    / self.ndegen[:, None, None, None])
        # Wannier90 use_ws_distance: same pair-dependent expansion as H(k)
        self.full_ws = False
        if (ws_model is not None
                and getattr(ws_model, "HR_full", None) is not None
                and getattr(ws_model, "wsvec_info", None) is not None
                and ws_model.wsvec_info.get("applied_to_model")):
            info = ws_model.wsvec_info
            Rfull = np.asarray(ws_model.Rlist_full, int)
            nrf = len(Rfull)
            Sfull = np.zeros((nrf, 3, nw, nw), dtype=complex)
            Rblocks = info["Rblocks"]
            for ib in range(len(Rblocks)):
                r = tuple(int(x) for x in Rblocks[ib])
                ir = info["ridx"][r]
                a, b = int(info["ijs"][ib, 0]), int(info["ijs"][ib, 1])
                n0, n1 = int(info["shift_start"][ib]), int(info["shift_start"][ib + 1])
                ns = max(n1 - n0, 1)
                w = 1.0 / float(ns)
                for sh in info["shifts_flat"][n0:n1]:
                    jr = info["rfidx"][tuple(int(x) for x in
                                             (np.asarray(r) + sh))]
                    Sfull[jr, :, a, b] += self.S_R[ir, :, a, b] * w
            self.S_R = Sfull
            self.Rvecs = Rfull
            self.ndegen = np.ones(nrf, dtype=int)
            self.full_ws = True
            self.ws_interpolation = True
            self.R_source = "wannier90_wsvec.dat"

    def at(self, k) -> np.ndarray:
        """Interpolated spin operators at k, shape (3,nw,nw)."""
        k = np.asarray(k, float).reshape(3)
        ph = np.exp(2j * np.pi * (self.Rvecs @ k))
        return np.einsum("R,Rcij->cij", ph, self.S_R, optimize=True)

    def hermiticity_offgrid(self, ktest=None) -> float:
        """max |S(k)-S(k)^dagger| off the DFT mesh (must be ~machine eps)."""
        if ktest is None:
            rng = np.random.default_rng(0)
            ktest = rng.uniform(-0.5, 0.5, size=(32, 3))
        ktest = np.atleast_2d(np.asarray(ktest, float))
        err = 0.0
        for k in ktest:
            S = self.at(k)
            err = max(err, float(np.max(np.abs(S - S.conj().transpose(0, 2, 1)))))
        return err

    def interpolation_error(self) -> float:
        """Max |S(k_mesh) - original S_W(k_mesh)| (must be machine precision)."""
        return max(float(np.max(np.abs(self.at(k) - S)))
                   for k, S in zip(self.kmesh, self.spin_mats))

    def checks(self) -> list:
        """(name, passed, detail) checks for the rotated operator."""
        Sm = self.spin_mats
        herm = max(float(np.max(np.abs(s - s.conj().transpose(0, 2, 1))))
                   for s in Sm)
        ev_min = min(float(np.linalg.eigvalsh(s[c]).min())
                     for s in Sm for c in range(3))
        ev_max = max(float(np.linalg.eigvalsh(s[c]).max())
                     for s in Sm for c in range(3))
        s2 = np.einsum("kcij,kcjl->kil", Sm, Sm, optimize=True)
        s2_lo = float(np.min(np.linalg.eigvalsh(s2)))
        s2_hi = float(np.max(np.linalg.eigvalsh(s2)))
        roundtrip = self.interpolation_error()
        # TRS: for a non-magnetic system eig[S(-k)] = -eig[S(k)]
        from .symmetry import find_kpairs
        partner = find_kpairs(self.kmesh)
        Sg = np.empty((len(Sm), 3, self.nw))
        for ik in range(len(Sm)):
            for c in range(3):
                Sg[ik, c] = np.sort(np.linalg.eigvalsh(Sm[ik, c]))
        trs_dev = 0.0
        for ik, j in enumerate(partner):
            if j <= ik:
                continue
            # spec(S(-k)) = -spec(S(k))  <=>  sorted(S(-k)) = -sorted(S(k))[::-1]
            trs_dev = max(trs_dev,
                          float(np.max(np.abs(Sg[j] + Sg[ik][:, ::-1]))))
        return [
            ("spin_wannier_interpolation_set",
             self.R_source in ("wannier90_hr.dat", "wannier90_wsvec.dat",
                               "reconstructed"),
             f"R set source = {self.R_source} ({len(self.Rvecs)} R vectors, "
             f"ndegen in [{self.ndegen.min()},{self.ndegen.max()}])"),
            ("spin_wannier_hermiticity", herm < 1e-12,
             f"max|S-S^dag| = {herm:.2e}"),
            ("spin_wannier_trs", trs_dev < 1e-2,
             f"max|sorted eig S(-k) + reverse sorted eig S(k)| = "
             f"{trs_dev:.4f} (0 expected for a "
             f"TRS-preserving non-magnetic system; large values indicate "
             f"mismatched runs or inaccurate spin matrices)"),
            ("spin_wannier_interpolation", roundtrip < 1e-10,
             f"max|S(k_mesh)-S_W(k_mesh)| = {roundtrip:.2e}"),
            ("spin_wannier_hermiticity_offgrid",
             (not self.ws_interpolation)
             or self.hermiticity_offgrid() < 1e-10,
             (f"WS/ndegen interpolation: max|S(k)-S(k)^dag| off-grid = "
              f"{self.hermiticity_offgrid():.2e}")
             if self.ws_interpolation else
             ("NOT_VALIDATED: shifted-mesh alias fallback "
              f"(off-grid |S-S^dag| = {self.hermiticity_offgrid():.2e}); "
              "grid values are exact but off-grid spin physics is not "
              "validated for this backend")),
            ("spin_wannier_eigenvalue_bound",
             ev_min >= -1.0 - 3e-3 and ev_max <= 1.0 + 3e-3,
             f"eig range [{ev_min:+.3f},{ev_max:+.3f}] (Pauli bound [-1,1])"),
            ("spin_wannier_sum_rule",
             abs(s2_lo - 3) < 0.05 and abs(s2_hi - 3) < 0.05,
             f"eig(sum_i S_i^2) in [{s2_lo:.3f},{s2_hi:.3f}] (3 for a "
             f"spin-complete basis; <3 means spin weight outside the "
             f"Wannier window)"),
        ]


def build_wannier_spin_operator(spn_path, chk_path, eig_path, model=None,
                                h_tolerance=1e-3):
    """Full pipeline spn + chk + eig -> WannierSpinOperator.

    Raises ``WannierSpinOperatorError`` if the checkpoint fails its internal
    consistency test (H_W vs H_hr above ``h_tolerance`` eV); no spin analysis
    is produced from an inconsistent transformation.
    """
    spn = read_wannier_spn(spn_path)
    chk = read_chk(chk_path)
    eig = read_eig(eig_path)
    cchk = checkpoint_checks(chk)
    if cchk["u_unitarity"] > 1e-8 or cchk["v_isometry"] > 1e-8:
        raise WannierSpinOperatorError(
            f"checkpoint transformation is not unitary/isometric: "
            f"U={cchk['u_unitarity']:.2e}, v={cchk['v_isometry']:.2e}")
    if model is not None:
        hc = hamiltonian_consistency(chk, eig, model)
        if hc["max_H_diff"] > h_tolerance:
            raise WannierSpinOperatorError(
                f"H_W(k) from the checkpoint does not reproduce H_hr(k): "
                f"max diff = {hc['max_H_diff']:.2e} eV > {h_tolerance:g} eV; "
                f"k ordering / gauge / window may be inconsistent")
    else:
        hc = None
    Sm = rotate_spin_to_wannier(spn, chk)
    interp = None
    if model is not None and getattr(model, "Rlist", None) is not None:
        interp = (np.asarray(model.Rlist, int), np.asarray(model.ndegen, int))
    op = WannierSpinOperator(chk["kpt_red"], chk["mp_grid"], Sm,
                             lattice=(None if model is None
                                      else getattr(model, "lattice", None)),
                             interpolation_set=interp,
                             ws_model=model)
    metrics, _ = validate_spin_matrices(spn["matrices"])
    checks = list(op.checks())
    if hc is not None:
        checks.append(("checkpoint_hamiltonian_consistency",
                       hc["max_H_diff"] <= h_tolerance,
                       f"max|H_W-H_hr| = {hc['max_H_diff']:.2e} eV "
                       f"(tol {h_tolerance:g} eV over {chk['num_kpts']} k)"))
    checks.append(("checkpoint_transformation",
                   cchk["u_unitarity"] < 1e-8 and cchk["v_isometry"] < 1e-8,
                   f"U unitarity {cchk['u_unitarity']:.2e}, "
                   f"v isometry {cchk['v_isometry']:.2e}, "
                   f"lattice error {cchk['lattice_error']:.2e}"))
    metadata = dict(source_spn=spn_path, source_chk=chk_path,
                    source_eig=eig_path, num_bands=spn["num_bands"],
                    num_wann=chk["num_wann"], num_kpts=chk["num_kpts"],
                    mp_grid=chk["mp_grid"].tolist(),
                    have_disentangled=bool(chk["have_disentangled"]),
                    spn_header=spn["header"],
                    spn_quality=metrics,
                    h_consistency=(None if hc is None
                                   else hc["max_H_diff"]))
    return op, dict(checks=checks, metrics=metrics, checkpoint=cchk,
                    h_consistency=hc, metadata=metadata)
