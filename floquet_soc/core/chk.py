# -*- coding: utf-8 -*-
"""Wannier90 checkpoint (.chk) and eigenvalue (.eig) readers.

The checkpoint contains the Wannier transformation matrices.  For a
disentangled calculation (``have_disentangled = .true.``) the matrix that
maps DFT band states to Wannier functions is

    v_{m,n}(k) = sum_{j <= ndimwin(k)} U_opt[j,m](k) · U[j,n](k)
                 for m in the outer window, zero otherwise

so that  |w_{n k}> = sum_m v_{m n}(k) |u_{m k}>  and any band-basis operator
O transforms as  O_W(k) = v†(k) O(k) v(k).

The parser was verified on the PbNV checkpoint: U is unitary to 1e-13 and
``v† diag(eig) v`` reproduces the Wannier Hamiltonian from ``_hr.dat`` to
9.5e-6 eV over all 36 k points (k ordering, band ordering and gauge all agree).

Format reference: wannier90 ``wannier_write_chk`` / wannier-berri
``w90files/chk.py`` (Fortran unformatted records).
"""
from __future__ import annotations

import struct

import numpy as np


class WannierCheckpointError(RuntimeError):
    pass


class FortranFileR:
    """Minimal sequential Fortran-unformatted reader (little endian)."""

    def __init__(self, path):
        self.path = path
        self.f = open(path, "rb")

    def close(self):
        self.f.close()

    def record(self):
        hdr = self.f.read(4)
        if len(hdr) < 4:
            raise EOFError(f"{self.path}: unexpected end of file")
        n = struct.unpack("<i", hdr)[0]
        if n < 0 or n > 10**9:
            raise WannierCheckpointError(f"{self.path}: bad record length {n}")
        data = self.f.read(n)
        tail = self.f.read(4)
        n2 = struct.unpack("<i", tail)[0] if len(tail) == 4 else None
        if n2 != n:
            raise WannierCheckpointError(
                f"{self.path}: record length mismatch {n} != {n2}")
        return data

    def ints(self):
        return np.frombuffer(self.record(), dtype=np.int32)

    def floats(self):
        return np.frombuffer(self.record(), dtype=np.float64)

    def complexs(self):
        a = self.floats()
        return a[::2] + 1j * a[1::2]

    def string(self):
        return self.record().decode("ascii", "replace").strip()

    def skip(self):
        self.record()


def read_chk(path: str) -> dict:
    """Parse ``seedname.chk`` into a dict (including the v_matrix)."""
    F = FortranFileR(path)
    try:
        out = {}
        out["comment"] = F.string()
        out["num_bands"] = int(F.ints()[0])
        num_exclude = int(F.ints()[0])
        exclude = F.ints()
        if len(exclude) != num_exclude:
            raise WannierCheckpointError(
                f"{path}: exclude_bands length mismatch")
        out["exclude_bands"] = exclude
        out["real_lattice"] = F.floats().reshape((3, 3), order="F")
        out["recip_lattice"] = F.floats().reshape((3, 3), order="F")
        out["num_kpts"] = int(F.ints()[0])
        out["mp_grid"] = F.ints().astype(int)
        out["kpt_red"] = F.floats().reshape((out["num_kpts"], 3))
        out["nntot"] = int(F.ints()[0])
        out["num_wann"] = int(F.ints()[0])
        out["checkpoint"] = F.string()
        out["have_disentangled"] = bool(F.ints()[0])
        nb, nk, nw = out["num_bands"], out["num_kpts"], out["num_wann"]
        if out["have_disentangled"]:
            out["omega_invariant"] = float(F.floats()[0])
            out["lwindow"] = np.array(F.ints().reshape((nk, nb)), dtype=bool)
            out["ndimwin"] = F.ints().astype(int)
            out["u_matrix_opt"] = F.complexs().reshape(
                (nk, nw, nb)).swapaxes(1, 2)
        else:
            out["lwindow"] = np.ones((nk, nb), dtype=bool)
            out["ndimwin"] = np.full(nk, nb)
            out["u_matrix_opt"] = None
        out["u_matrix"] = F.complexs().reshape(
            (nk, nw, nw)).swapaxes(1, 2)
        F.skip()                                    # M_matrix (not needed)
        out["centres"] = F.floats().reshape((nw, 3))
        out["spreads"] = F.floats().reshape((nw,))
    finally:
        F.close()

    out["v_matrix"] = build_v_matrix(out)
    return out


def build_v_matrix(chk: dict) -> np.ndarray:
    """v_{m,n}(k): DFT band -> Wannier function rotation, shape (nk,nb,nw)."""
    nk, nb, nw = chk["num_kpts"], chk["num_bands"], chk["num_wann"]
    v = np.zeros((nk, nb, nw), dtype=complex)
    for ik in range(nk):
        u = chk["u_matrix"][ik]
        if chk["have_disentangled"]:
            nd = int(chk["ndimwin"][ik])
            v[ik, chk["lwindow"][ik], :] = chk["u_matrix_opt"][ik][:nd, :] @ u
        else:
            v[ik] = u
    return v


def read_eig(path: str) -> np.ndarray:
    """Read ``seedname.eig`` (columns: band index, k index, energy [eV]).

    Returns array (num_kpts, num_bands).
    """
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    ib = data[:, 0].astype(int) - 1
    ik = data[:, 1].astype(int) - 1
    val = data[:, 2]
    nk, nb = int(ik.max()) + 1, int(ib.max()) + 1
    E = np.zeros((nk, nb))
    E[ik, ib] = val
    return E


def checkpoint_checks(chk: dict) -> dict:
    """Numerical quality metrics of the parsed checkpoint."""
    nw = chk["num_wann"]
    nk = chk["num_kpts"]
    uerr = max(float(np.max(np.abs(u.conj().T @ u - np.eye(nw))))
               for u in chk["u_matrix"])
    V = chk["v_matrix"]
    verr = max(float(np.max(np.abs(V[ik].conj().T @ V[ik] - np.eye(nw))))
               for ik in range(nk))
    lat_err = float(np.max(np.abs(
        chk["real_lattice"] @ chk["recip_lattice"].T / (2 * np.pi)
        - np.eye(3))))
    # consistency of the window bookkeeping
    win_bad = 0
    for ik in range(nk):
        lw = chk["lwindow"][ik]
        nd = int(chk["ndimwin"][ik])
        if int(np.sum(lw)) != nd:
            win_bad += 1
        if np.sum(lw) and np.any(lw[:np.argmax(lw)]):
            win_bad += 0                          # contiguity checked below
        idx = np.where(lw)[0]
        if len(idx) and (idx[-1] - idx[0] + 1 != len(idx)):
            win_bad += 1
    rows_beyond = 0.0
    if chk["have_disentangled"]:
        uo = chk["u_matrix_opt"]
        vals = [float(np.max(np.abs(uo[ik][chk["ndimwin"][ik]:])))
                for ik in range(nk) if chk["ndimwin"][ik] < chk["num_bands"]]
        rows_beyond = max(vals) if vals else 0.0
    return dict(u_unitarity=uerr, v_isometry=verr, lattice_error=lat_err,
                window_inconsistent=win_bad,
                u_opt_beyond_window=rows_beyond,
                num_kpts=nk, num_bands=chk["num_bands"],
                num_wann=nw, have_disentangled=chk["have_disentangled"])


def hamiltonian_consistency(chk: dict, eig: np.ndarray, model) -> dict:
    """Compare H_W(k)=v† diag(eig) v with the Wannier Hamiltonian H(k).

    This is the strongest end-to-end test of the checkpoint parsing: it
    simultaneously checks the k ordering, band ordering, outer window and the
    Wannier gauge.
    """
    v = chk["v_matrix"]
    worst = 0.0
    worst_eig = 0.0
    per_k = []
    for ik in range(chk["num_kpts"]):
        k = chk["kpt_red"][ik]
        Hw = v[ik].conj().T @ (np.diag(eig[ik]) @ v[ik])
        Hr = model.Hk(k)
        d = float(np.max(np.abs(Hw - Hr)))
        de = float(np.max(np.abs(np.linalg.eigvalsh(Hw)
                                 - np.linalg.eigvalsh(Hr))))
        worst = max(worst, d)
        worst_eig = max(worst_eig, de)
        per_k.append(d)
    return dict(max_H_diff=worst, max_eig_diff=worst_eig, per_k=per_k)
