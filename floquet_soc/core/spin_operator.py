# -*- coding: utf-8 -*-
"""Spin operators for Wannier/analytic models.

Two sources are supported:

1. ``pauli`` — exact spin-1/2 operators ``S_i = sigma_i / 2`` for analytic
   2-band models (Rashba, Dresselhaus, Kane-Mele, ...).  Benchmark-verified.

2. ``spn`` — spin matrices read from a VASP-produced ``wannier90.spn`` file.
   The VASP binary format is decoded here.  IMPORTANT: the file produced by
   ``vaspspn`` (header "Created from WAVECAR at ...") contains matrices in the
   *WAVECAR band basis* (num_bands x num_bands, typically > num_wann), NOT in
   the Wannier basis.  Rotating them into the Wannier basis requires the
   Wannier transformation matrices (wannier90.chk / .u.mat / .u_dis.mat).
   If those files are missing the framework refuses to guess (see
   ``MissingTransformError``).

Gauge-consistency rule
----------------------
A spin operator may only be combined with H(k) if both are expressed in the
same basis and the same k-point gauge.  ``load_spin_for_model`` enforces this
and records the decision in the returned metadata.
"""
from __future__ import annotations

import os
import struct

import numpy as np

SIGMA = np.array([[[0, 1], [1, 0]],
                  [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)
#: spin-1/2 operators S_i = sigma_i/2 (units of hbar)
SPIN_HALF = 0.5 * SIGMA


class SpinInputError(RuntimeError):
    pass


class MissingTransformError(RuntimeError):
    """Raised when spin matrices cannot be brought to the Wannier basis."""

    def __init__(self, message, missing_files=None):
        super().__init__(message)
        self.missing_files = list(missing_files or [])


# ------------------------------------------------------------------ helpers
def _is_hermitian(M, tol=1e-12):
    return float(np.max(np.abs(M - M.conj().T))) <= tol


def validate_spin_matrices(mats, label="spin"):
    """Numerical quality checks for a set of (nk,3,nb,nb) spin matrices.

    Returns a dict of metrics plus a list of (name, passed, detail) checks.
    Checks that must hold for *any* representation in an orthonormal basis:
      * Hermiticity;
      * real diagonal;
      * eigenvalues within [-1, 1] when matrices are Pauli-like (S=sigma/2
        would be [-0.5,0.5], but VASP .spn convention is Pauli = sigma).
    """
    mats = np.asarray(mats)
    if mats.ndim != 4 or mats.shape[1] != 3 or mats.shape[2] != mats.shape[3]:
        raise ValueError(f"{label}: expected (nk,3,n,n), got {mats.shape}")
    nk, _, nb, _ = mats.shape
    herm = max(float(np.max(np.abs(m - m.conj().T)))
               for m in mats.reshape(-1, nb, nb))
    diag_im = max(float(np.max(np.abs(np.imag(np.diag(m)))))
                  for m in mats.reshape(-1, nb, nb))
    eig_min, eig_max = np.inf, -np.inf
    trace_max = 0.0
    s2_lo, s2_hi = np.inf, -np.inf
    diag_norm_max = 0.0
    for ik in range(nk):
        S = mats[ik]
        for c in range(3):
            ev = np.linalg.eigvalsh(S[c])
            eig_min = min(eig_min, float(ev.min()))
            eig_max = max(eig_max, float(ev.max()))
            trace_max = max(trace_max, abs(float(np.trace(S[c]).real)))
        S2 = S[0] @ S[0] + S[1] @ S[1] + S[2] @ S[2]
        ev2 = np.linalg.eigvalsh(S2)
        s2_lo, s2_hi = min(s2_lo, float(ev2.min())), max(s2_hi, float(ev2.max()))
        pol = np.sqrt(sum(np.diag(S[c].real) ** 2 for c in range(3)))
        diag_norm_max = max(diag_norm_max, float(np.max(pol)))

    checks = [
        ("spin_hermiticity", herm < 1e-10,
         f"max|S-S^dag| = {herm:.2e} (tolerance 1e-10)"),
        ("spin_real_diagonal", diag_im < 1e-10,
         f"max|Im diag| = {diag_im:.2e} (tolerance 1e-10)"),
        ("spin_eigenvalue_bound", (eig_min >= -1.0 - 1e-3)
         and (eig_max <= 1.0 + 1e-3),
         f"eigenvalues in [{eig_min:+.3f},{eig_max:+.3f}], Pauli bound [-1,1]"),
        ("spin_diagonal_polarization_bound", diag_norm_max <= 1.0 + 1e-3,
         f"max |<sigma>_band| = {diag_norm_max:.3f} (Pauli bound 1)"),
    ]
    metrics = dict(hermiticity=herm, diag_imag=diag_im,
                   eig_min=float(eig_min), eig_max=float(eig_max),
                   trace_max=trace_max, S2_min=float(s2_lo),
                   S2_max=float(s2_hi), diag_pol_max=float(diag_norm_max),
                   num_bands=int(nb), num_kpts=int(nk))
    return metrics, checks


# ------------------------------------------------------------------ VASP spn
def read_vasp_spn(path: str) -> dict:
    """Decode the VASP/WAVECAR-derived binary ``wannier90.spn`` file.

    Layout (verified against the writer in
    wannier-berri/utils/vaspspn.py):

    * Fortran-unformatted record 0: 60-byte ASCII header;
    * record 1: int32 [num_bands, num_kpts];
    * one record per k point: 3*N*(N+1)/2 complex128 values, ordered as
      for m in range(N):
          for n in range(m+1):
              SX[n,m], SY[n,m], SZ[n,m]
      i.e. upper triangle (column major) with the 3 spin components
      interleaved.  Matrices are reconstructed as Hermitian.
    """
    with open(path, "rb") as f:

        def rec():
            raw = f.read(4)
            if len(raw) < 4:
                raise SpinInputError(f"{path}: truncated record header")
            n = struct.unpack("<i", raw)[0]
            if n < 0 or n > 10**9:
                raise SpinInputError(f"{path}: bad record length {n}")
            data = f.read(n)
            tail = f.read(4)
            if len(tail) != 4 or struct.unpack("<i", tail)[0] != n:
                raise SpinInputError(f"{path}: record length mismatch")
            return data

        header = rec().decode("ascii", "replace").strip()
        dims = rec()
        if len(dims) != 8:
            raise SpinInputError(f"{path}: expected 2 int32 in dims record")
        num_bands, num_kpts = struct.unpack("<ii", dims)
        tri = num_bands * (num_bands + 1) // 2
        expected = 16 * 3 * tri
        mats = np.empty((num_kpts, 3, num_bands, num_bands), dtype=complex)
        idx = [(i, j) for j in range(num_bands) for i in range(j + 1)]
        for ik in range(num_kpts):
            data = rec()
            if len(data) != expected:
                raise SpinInputError(
                    f"{path}: k-record {ik} has {len(data)} bytes, expected "
                    f"{expected} (=3*{num_bands}*{num_bands+1}/2*16)")
            vec = np.frombuffer(data, dtype=np.complex128).reshape(-1, 3)
            for c in range(3):
                M = np.zeros((num_bands, num_bands), dtype=complex)
                for k, (i, j) in enumerate(idx):
                    M[i, j] = vec[k, c]
                diag = np.diagonal(M).copy()
                M = M + M.conj().T - np.diag(np.diag(M).real)
                np.fill_diagonal(M, diag)     # Wannier90 keeps diag as written
                mats[ik, c] = M
    return dict(header=header, num_bands=num_bands, num_kpts=num_kpts,
                matrices=mats, path=path, format="vasp-binary-spn")


def read_formatted_spn(path: str) -> dict:
    """Read a *formatted* Wannier90 ``.spn`` file (``spn_formatted = .true.``).

    Format (Wannier90 ``conv_read_spn_fmt``): a header line, then
    ``num_bands num_kpts``, then for each k point, for m in 1..num_bands and
    n in 1..m, three complex numbers (Sx, Sy, Sz) for element (n,m).  The
    upper triangle is stored; the lower triangle is filled by conjugation.
    """
    with open(path, "r", encoding="ascii", errors="replace") as f:
        header = f.readline().rstrip()
        dims = f.readline().split()
        if len(dims) < 2:
            raise SpinInputError(f"{path}: cannot read 'num_bands num_kpts'")
        num_bands, num_kpts = int(dims[0]), int(dims[1])
        rest = f.read()
    vals = np.fromstring(rest, sep=" ")
    n_complex = len(vals) // 2
    tri = num_bands * (num_bands + 1) // 2
    if n_complex != 3 * tri * num_kpts:
        raise SpinInputError(
            f"{path}: expected {3 * tri * num_kpts} complex entries, found "
            f"{n_complex}; is this really a formatted spn file?")
    cplx = vals[: 2 * n_complex].reshape(-1, 2)
    cplx = cplx[:, 0] + 1j * cplx[:, 1]
    idx = [(i, j) for j in range(num_bands) for i in range(j + 1)]
    mats = np.zeros((num_kpts, 3, num_bands, num_bands), dtype=complex)
    k = 0
    for ik in range(num_kpts):
        for (i, j) in idx:
            for c in range(3):
                mats[ik, c, i, j] = cplx[k]
                k += 1
        for c in range(3):
            M = mats[ik, c]
            diag = np.diagonal(M).copy()
            M = M + M.conj().T - np.diag(np.diag(M).real)
            np.fill_diagonal(M, diag)         # keep diag exactly as written
            mats[ik, c] = M
    return dict(header=header, num_bands=num_bands, num_kpts=num_kpts,
                matrices=mats, path=path, format="formatted-spn")


def read_wannier_spn(path: str) -> dict:
    """Auto-detect and read either a formatted or a binary (VASP) spn file."""
    with open(path, "rb") as f:
        head = f.read(8)
    is_text = bool(head) and all(32 <= b < 127 or b in (9, 10, 13)
                                 for b in head)
    if is_text:
        return read_formatted_spn(path)
    return read_vasp_spn(path)


#: files that would allow rotating band-basis spin matrices into Wannier basis
TRANSFORM_CANDIDATES = ("wannier90.chk", "wannier90.u.mat",
                        "wannier90.u_dis.mat", "wannier90.amn",
                        "wannier90.mmn", "wannier90.eig", "wannier90.nnkp")


def find_transform_files(spn_path: str) -> dict:
    d = os.path.dirname(os.path.abspath(spn_path)) or "."
    return {name: os.path.join(d, name) for name in TRANSFORM_CANDIDATES
            if os.path.exists(os.path.join(d, name))}


def _pauli_spin_at_k(_k=None):
    """``spin_at_k(k) -> S = sigma/2`` for a 2-band analytic model.

    The Pauli source has no k dependence, but every consumer
    (``run_points``, ``gauge_invariance_test``, ``compute_map``,
    ``dft_spin.compare``, ``soc_scan``) expects a *callable*.  Without this
    adapter those commands reported BLOCKED on analytic models while
    ``spin_static``/``pst``/``spin_floquet`` worked, and
    ``pst.gauge_invariance_test`` crashed calling an ndarray.
    """
    return SPIN_HALF


def load_spin_for_model(config: dict, model=None) -> dict:
    """Load spin operators for a model according to the material config.

    Returns dict:
        ops       (3,n,n) Pauli-spin operators (analytic) or None
        spn       parsed spn data (or None)
        basis     "pauli" | "wannier" | "vasp-bands"
        status    "ok" | "blocked"
        metadata  provenance/quality information
        checks    list of (name, passed, detail)
    Raises MissingTransformError when a rotation is required but the
    transformation files are absent.
    """
    spin_cfg = config.get("spin", {})
    source = spin_cfg.get("source", "pauli")
    checks = []
    if source == "pauli":
        nw_model = getattr(model, "nw", None)
        if nw_model is not None and nw_model != 2:
            raise SpinInputError(
                f"Pauli spin operators S=sigma/2 are only valid for 2-band "
                f"models; this model has nw={nw_model}. Provide wannier90.spn "
                f"+ .chk/.eig for a material spin operator.")
        return dict(ops=SPIN_HALF, spin_at_k=_pauli_spin_at_k, spn=None,
                    basis="pauli", status="ok",
                    metadata=dict(source="analytic Pauli matrices S=sigma/2",
                                  units="hbar"),
                    checks=[("spin_operator_source_pauli", True,
                             "analytic spin-1/2 operators S=sigma/2")])
    if source != "spn":
        raise SpinInputError(f"unknown spin source {source!r}")

    spn_path = spin_cfg.get("file")
    if not spn_path or not os.path.exists(spn_path):
        raise SpinInputError(f"spin source 'spn' but file not found: {spn_path}")
    spn = read_wannier_spn(spn_path)
    metrics, spin_checks = validate_spin_matrices(spn["matrices"],
                                                  label=os.path.basename(spn_path))
    checks.extend(spin_checks)
    nw = getattr(model, "nw", config.get("num_wann"))
    basis = "vasp-bands" if (nw is None or spn["num_bands"] != nw) else "wannier"
    metadata = dict(source=os.path.basename(spn_path), format=spn["format"],
                    header=spn["header"], num_bands=spn["num_bands"],
                    num_kpts=spn["num_kpts"], num_wann=nw, basis=basis,
                    metrics=metrics, units="Pauli (eigenvalue bound 1)")

    if basis == "wannier":
        # unlikely for VASP files (num_bands==num_wann); accept with warning
        checks.append(("spin_basis_matches_wannier", True,
                       "num_bands == num_wann; still verify gauge manually"))
        return dict(ops=None, spn=spn, basis="wannier", status="ok",
                    metadata=metadata, checks=checks)

    # ---- rotation to the Wannier basis needs the checkpoint + eigenvalues
    d = os.path.dirname(os.path.abspath(spn_path))
    chk_path = spin_cfg.get("chk") or os.path.join(d, "wannier90.chk")
    eig_path = spin_cfg.get("eig") or os.path.join(d, "wannier90.eig")
    have_chk = os.path.exists(chk_path)
    have_eig = os.path.exists(eig_path)
    if not (have_chk and have_eig):
        missing = [p for p, ok in ((chk_path, have_chk), (eig_path, have_eig))
                   if not ok]
        candidates = [p for p in (
            os.path.join(d, "wannier90.u.mat"),
            os.path.join(d, "wannier90.u_dis.mat"),
            os.path.join(d, "wannier90_u.mat"),
            os.path.join(d, "wannier90_u_dis.mat")) if os.path.exists(p)]
        extra = ("; found U-matrix files "
                 + ", ".join(os.path.basename(p) for p in candidates)
                 + " but .u.mat rotation is not implemented yet; provide "
                   "wannier90.chk + wannier90.eig"
                 ) if candidates else ""
        msg = (
            f"{os.path.basename(spn_path)} has {spn['num_bands']} bands but the "
            f"Wannier model has num_wann={nw}: the spin matrices are in the "
            f"WAVECAR band basis, not the Wannier basis. Rotation requires "
            f"wannier90.chk (U matrices) and wannier90.eig; missing: "
            f"{', '.join(missing)}{extra}")
        raise MissingTransformError(msg, missing_files=missing)
    from .spin_wannier import (build_wannier_spin_operator,
                               WannierSpinOperatorError)
    try:
        op, info = build_wannier_spin_operator(spn_path, chk_path, eig_path,
                                               model=model)
    except WannierSpinOperatorError as exc:
        raise SpinInputError(f"spin-operator rotation failed: {exc}") from exc
    checks.extend(info["checks"])
    metadata = dict(source=os.path.basename(spn_path),
                    format=spn["format"], header=spn["header"],
                    num_bands=spn["num_bands"], num_kpts=spn["num_kpts"],
                    num_wann=nw, basis="wannier",
                    transform=os.path.basename(chk_path),
                    eig=os.path.basename(eig_path),
                    h_consistency_ev=(None if info["h_consistency"] is None
                                      else info["h_consistency"]["max_H_diff"]),
                    metrics=info["metrics"],
                    spn_quality=info["metrics"],
                    units="Pauli (eigenvalue bound 1)")
    return dict(ops=None, spin_at_k=op.at, spin_operator=op, spn=spn,
                basis="wannier", status="ok", metadata=metadata,
                checks=checks)


def spin_texture_from_spn(spn: dict, k_indices=None) -> np.ndarray:
    """Diagonal spin polarization per band from spn matrices, (nk,nb,3).

    Valid only if the spn matrices are in the same basis as the band
    eigenstates used for interpretation; here only meaningful for the VASP
    band basis (it is the spin expectation of each WAVECAR band).
    """
    mats = spn["matrices"]
    nk, _, nb, _ = mats.shape
    if k_indices is not None:
        mats = mats[np.asarray(k_indices)]
        nk = mats.shape[0]
    tex = np.empty((nk, nb, 3))
    for ik in range(nk):
        for c in range(3):
            tex[ik, :, c] = np.real(np.diag(mats[ik, c]))
    return tex
