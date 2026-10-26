# -*- coding: utf-8 -*-
"""Generic Wannier90 tight-binding reader.

No hard coding of num_wann / num_electrons / orbital types: everything is read
from the input files (``.win``, ``_hr.dat``, ``_centres.xyz``).

Conventions
-----------
* R vectors are integer lattice translations (in units of the real-space
  lattice vectors given in ``.win``).
* k points are fractional (crystal) coordinates; the Bloch phase is
  ``exp(2*pi*i * k . R)`` (Wannier90 convention).
* ``H(k) = sum_R H(R) exp(2*pi*i k.R) / ndegen(R)``.
"""
from __future__ import annotations

import os
import re

import numpy as np


def ws_supercell(ngrid, lattice):
    """Wannier90 WS interpolation set R and degeneracy factors ndegen.

    R_i runs over ``-(N_i//2) .. N_i//2``: N+1 values for even N and N values
    for odd N (e.g. N=15 -> -7..7, N=16 -> -8..8), which is the Wannier90
    convention and makes the set symmetric under R -> -R for both parities.
    ndegen(R) = number of R' in this set with R' == R (mod N).
    For an unshifted k mesh the interpolant
    H(k)=sum_R H(R)exp(2*pi*i k.R)/ndegen(R) reproduces H(k) on the mesh and
    is Hermitian for every k (the set is symmetric under R -> -R and ndegen
    is symmetric).

    Notes
    -----
    v0.1.5-rc1 fix: the previous code used ``np.arange(-N // 2, ...)``; for odd
    N Python floor division gives ``-N//2 == -(N+1)//2``, which produced an
    asymmetric, oversized R set (N=15 -> 16 values -8..7; N=1 -> 2 values
    -1..0).  For the MoSSe 15x15x1 mesh that made the spin-operator
    interpolation use 512 R points with ndegen 2/4/8 instead of the 225
    points with ndegen 1 written by Wannier90, so the off-grid spin operator
    was non-Hermitian (max|S-S^dag| ~ 1).  The H(k) reader was unaffected
    (it uses the R list from ``_hr.dat``), but the spin operator and
    ``from_chk_eig`` use this helper.
    """
    from collections import Counter
    import numpy as np
    ngrid = [int(x) for x in ngrid]
    ranges = [np.arange(-(N // 2), N // 2 + 1) for N in ngrid]
    R = np.array(np.meshgrid(*ranges, indexing="ij")).reshape(3, -1).T
    R = R.astype(int)
    key = np.mod(R, np.asarray(ngrid, int))
    cnt = Counter(tuple(int(v) for v in row) for row in key)
    ndegen = np.array([cnt[tuple(int(v) for v in row)] for row in key], int)
    return R, ndegen


class WannierInputError(RuntimeError):
    pass


class HermiticityError(WannierInputError):
    """H(k) is not Hermitian beyond the accepted tolerance.

    Raised instead of silently averaging, because a non-Hermitian H(k) means
    the real-space data (``_hr.dat`` R/-R pairs, or ``wsvec.dat`` coverage) is
    wrong; ``np.linalg.eigvalsh`` would otherwise read one triangle and return
    plausible-looking but wrong eigenvalues.
    """


def assert_hermitian(H, k=None, rtol: float = 1e-10, what: str = "H(k)") -> float:
    """Return max|H - H^dag|, raising :class:`HermiticityError` if too large.

    The threshold is *relative* (``rtol`` times the matrix scale, with an
    absolute floor of 1) because the residual of a legitimate Fourier sum
    grows with |H|: a 156x156, |H| ~ 57 eV model sits at ~1e-14, twenty
    orders below the default threshold.

    Once this passes, ``eigvalsh`` (one triangle, ~2x faster than ``eigh``)
    is a legitimate choice rather than a hidden assumption.
    """
    H = np.asarray(H)
    res = float(np.max(np.abs(H - np.conj(np.swapaxes(H, -1, -2)))))
    scale = max(1.0, float(np.max(np.abs(H))))
    tol = rtol * scale
    if res > tol:
        where = "" if k is None else f" at k={np.asarray(k, float).tolist()}"
        raise HermiticityError(
            f"{what} not Hermitian{where}: max|H - H^dag| = {res:.3e} "
            f"> tol {tol:.3e} (scale {scale:.3e}). Check the _hr.dat R/-R "
            f"pairs and the wannier90_wsvec.dat shift coverage; a one-sided "
            f"wsvec block set produces exactly this.")
    return res


# ------------------------------------------------------------------ .win
def _fortran_float(txt):
    """Parse a Fortran-style real (1.0d0, 1.0D-2, .5) as float."""
    return float(str(txt).strip().replace("D", "E").replace("d", "e"))


def read_win(path: str) -> dict:
    """Parse the subset of ``.win`` needed by the framework."""
    with open(path, encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    info = dict(path=path, num_wann=None, num_bands=None, spinors=False,
                lattice=None, atoms=[], kpoints=None, mp_grid=None,
                projections=[], use_bloch_phases=False,
                dis_win_min=None, dis_win_max=None,
                dis_froz_min=None, dis_froz_max=None)
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        s = raw.split("!")[0].strip()
        low = s.lower()
        if low.startswith("num_wann"):
            info["num_wann"] = int(s.split("=")[1].split()[0])
        elif low.startswith("num_bands"):
            info["num_bands"] = int(s.split("=")[1].split()[0])
        elif low.startswith("spinors"):
            info["spinors"] = "t" in s.split("=")[1].strip().lower()[:4]
        elif low.startswith("use_bloch_phases"):
            info["use_bloch_phases"] = "t" in s.split("=")[1].strip().lower()[:4]
        elif low.startswith("dis_win_min"):
            info["dis_win_min"] = _fortran_float(s.split("=")[1].split()[0])
        elif low.startswith("dis_win_max"):
            info["dis_win_max"] = _fortran_float(s.split("=")[1].split()[0])
        elif low.startswith("dis_froz_min"):
            info["dis_froz_min"] = _fortran_float(s.split("=")[1].split()[0])
        elif low.startswith("dis_froz_max"):
            info["dis_froz_max"] = _fortran_float(s.split("=")[1].split()[0])
        elif low.startswith("mp_grid"):
            info["mp_grid"] = [int(x) for x in s.split("=")[1].split()]
        elif low.startswith("begin unit_cell_cart"):
            i += 1
            lat = []
            while i < n and not lines[i].strip().lower().startswith("end"):
                p = lines[i].split()
                if len(p) >= 3:
                    lat.append([_fortran_float(x) for x in p[:3]])
                i += 1
            info["lattice"] = np.array(lat, float)
        elif low.startswith("begin atoms_cart"):
            i += 1
            atoms = []
            while i < n and not lines[i].strip().lower().startswith("end"):
                p = lines[i].split()
                if len(p) >= 4:
                    atoms.append((p[0], np.array([_fortran_float(x) for x in p[1:4]])))
                i += 1
            info["atoms"] = atoms
        elif low.startswith("begin kpoints"):
            i += 1
            kk = []
            while i < n and not lines[i].strip().lower().startswith("end"):
                p = lines[i].split()
                if len(p) >= 3:
                    kk.append([_fortran_float(x) for x in p[:3]])
                i += 1
            info["kpoints"] = np.array(kk, float)
        elif low.startswith("begin projections"):
            i += 1
            proj = []
            while i < n and not lines[i].strip().lower().startswith("end"):
                s2 = lines[i].strip()
                if s2:
                    proj.append(s2)
                i += 1
            info["projections"] = proj
        i += 1

    if info["num_wann"] is None:
        raise WannierInputError(f"{path}: num_wann not found")
    if info["lattice"] is None:
        raise WannierInputError(f"{path}: unit_cell_cart not found")
    return info


# ------------------------------------------------------------------ _hr.dat
def read_hr(path: str) -> dict:
    """Parse ``seedname_hr.dat`` (Wannier90 format)."""
    with open(path, encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if not lines:
        raise WannierInputError(f"{path}: empty")
    i = 0
    # optional header comment line
    def is_int(tok):
        return re.fullmatch(r"[+-]?\d+", tok) is not None

    if not is_int(lines[0]):
        i = 1
    nw = int(lines[i]); i += 1
    nrpts = int(lines[i]); i += 1
    # degeneracy lines: nrpts integers (15 per line)
    ndegen = []
    while len(ndegen) < nrpts:
        ndegen.extend(int(x) for x in lines[i].split())
        i += 1
    ndegen = np.array(ndegen[:nrpts], int)
    Rlist = np.zeros((nrpts, 3), int)
    H = np.zeros((nrpts, nw, nw), dtype=complex)

    # Two hr.dat dialects exist in the wild:
    #  (a) a separate "R1 R2 R3" line followed by nw*nw lines "m n Re Im";
    #  (b) every line carries "R1 R2 R3 m n Re Im" (VASP's writer).
    #         We support both, and reject anything inconsistent.
    ir = -1
    cur_R = None
    count = 0
    total = 0
    n_total = nrpts * nw * nw

    def _check_block_complete(block_R, block_count):
        # A short block silently shifted every following matrix element; that
        # is the worst failure mode (wrong numbers, no error), so it must be
        # rejected here rather than left to the leftover-lines check (which
        # only catches junk at the very end of the file).
        if block_R is not None and block_count != nw * nw:
            raise WannierInputError(
                f"{path}: R-block (R={tuple(block_R)}) has {block_count} "
                f"elements, expected {nw*nw} — hr.dat truncated or format "
                f"corrupted")

    while i < len(lines) and total < n_total:
        p = lines[i].split()
        if len(p) >= 7:
            r = (int(p[0]), int(p[1]), int(p[2]))
            m, n = int(p[3]) - 1, int(p[4]) - 1
            val = complex(float(p[5]), float(p[6]))
            if r != cur_R:
                _check_block_complete(cur_R, count)
                ir += 1
                if ir >= nrpts:
                    raise WannierInputError(
                        f"{path}: more R blocks than nrpts={nrpts}")
                cur_R = r
                Rlist[ir] = r
                count = 0
            if not (0 <= m < nw and 0 <= n < nw):
                raise WannierInputError(
                    f"{path}: band index out of range ({m},{n})")
            H[ir, m, n] = val
            count += 1
            total += 1
        elif len(p) == 3:
            _check_block_complete(cur_R, count)
            ir += 1
            if ir >= nrpts:
                raise WannierInputError(
                    f"{path}: more R blocks than nrpts={nrpts}")
            cur_R = (int(p[0]), int(p[1]), int(p[2]))
            Rlist[ir] = cur_R
            count = 0
        elif len(p) == 4 and cur_R is not None:
            m, n = int(p[0]) - 1, int(p[1]) - 1
            val = complex(float(p[2]), float(p[3]))
            if not (0 <= m < nw and 0 <= n < nw):
                raise WannierInputError(
                    f"{path}: band index out of range ({m},{n})")
            H[ir, m, n] = val
            count += 1
            total += 1
        else:
            raise WannierInputError(
                f"{path}: cannot parse line {i}: {lines[i][:80]!r}")
        i += 1
    _check_block_complete(cur_R, count)
    if ir != nrpts - 1:
        raise WannierInputError(
            f"{path}: parsed only {ir + 1} of {nrpts} R blocks")
    if total != n_total:
        raise WannierInputError(
            f"{path}: parsed {total} matrix elements, expected {n_total} "
            f"(nrpts={nrpts}, nw={nw})")
    if i != len(lines):
        raise WannierInputError(
            f"{path}: {len(lines)-i} leftover lines after parsing")
    return dict(num_wann=nw, nrpts=nrpts, ndegen=ndegen, R=Rlist, H=H,
                path=path)


def write_hr(path: str, num_wann: int, Rlist, ndegen, H) -> str:
    """Write a Wannier90 ``_hr.dat`` (dialect a: R header + nw*nw lines).

    Used to export derived real-space Hamiltonians (photon-resolved ``H_p(R)``
    or the light-dressed effective Hamiltonian) in a format any Wannier90 /
    tight-binding downstream code can read directly.

    ``H`` has shape ``(nrpts, num_wann, num_wann)``; ``ndegen`` has length
    ``nrpts`` (pass ones for an already-weighted/expanded representation).
    """
    Rlist = np.asarray(Rlist, int).reshape(-1, 3)
    H = np.asarray(H, complex)
    ndegen = np.asarray(ndegen, int).reshape(-1)
    nrpts = len(Rlist)
    if H.shape != (nrpts, num_wann, num_wann):
        raise WannierInputError(
            f"write_hr: H shape {H.shape} != ({nrpts},{num_wann},{num_wann})")
    if len(ndegen) != nrpts:
        raise WannierInputError(
            f"write_hr: ndegen length {len(ndegen)} != nrpts {nrpts}")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{num_wann}\n{nrpts}\n")
        for i in range(0, nrpts, 15):
            f.write(" ".join(str(int(x)) for x in ndegen[i:i + 15]) + "\n")
        for ir in range(nrpts):
            f.write(f"  {Rlist[ir, 0]:4d}{Rlist[ir, 1]:4d}{Rlist[ir, 2]:4d}\n")
            for m in range(num_wann):
                for n in range(num_wann):
                    v = H[ir, m, n]
                    f.write(f"    {m + 1:4d}{n + 1:4d}"
                            f"  {v.real:18.12f}  {v.imag:18.12f}\n")
    return path


# ------------------------------------------------------------------ wsvec
def header_flag(header: str, key: str):
    """Parse ``key = .true./.false.`` out of a Wannier90 comment header.

    Returns ``True``/``False`` when the token is present and ``None`` when it
    is absent, so the caller can apply a documented inference instead of
    silently assuming a value.  A real ``wannier90_wsvec.dat`` header often
    contains no tokens at all (`` use_ws_distance=.true. `` is optional), which
    is exactly the case that used to fall back to the *folded* convention and
    produce wrong off-grid H(k).
    """
    txt = header.replace(" ", "").replace("\t", "").lower()
    i = txt.find(key.lower())
    if i < 0:
        return None
    rest = txt[i + len(key):].lstrip("=")
    for token, value in ((".true.", True), (".false.", False),
                         ("true", True), ("false", False),
                         ("t", True), ("f", False)):
        if rest.startswith(token):
            return value
    return None


def read_wsvec(path: str, keep_vectors: bool = False):
    """Parse ``wannier90_wsvec.dat``.

    Format: one comment header, then blocks

        R1 R2 R3  i j
        n_wsvec
        v1 v2 v3
        ... (n_wsvec lines)

    The listed vectors are *supercell* translations (multiples of ``mp_grid``)
    that must be added to ``R`` to put the second Wannier function in the
    Wigner-Seitz cell of the first one (Wannier90 ``use_ws_distance``).

    Returns (all compact):
        header            : the raw comment line
        use_ws_distance   : True/False from the header, or None when the
                           header carries no token (callers must then
                           apply the documented inference)
        write_ndegen_applied : token in the header; if True the shifts are
                           already folded in the real-space files and must NOT
                           be applied again
        R                 : unique R set (nrpts,3)
        Rblocks           : R of every block (nblocks,3)
        ijs               : (nblocks,2) Wannier indices (0-based)
        n_wsvec           : (nblocks,) number of shifts per block
        shift_start       : (nblocks+1,) offsets into shifts_flat
        shifts_flat       : (ntot,3) all listed shifts (int32)
        vectors           : list of (n,3) arrays, only if keep_vectors
    """
    Rs = []
    ijs = []
    ns = []
    vec_list = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    header = lines[0] if lines else ""
    i = 1                                   # skip the comment header
    nlines = len(lines)
    while i < nlines:
        p = lines[i].split()
        if not p:
            i += 1
            continue
        if len(p) == 5:
            r = (int(p[0]), int(p[1]), int(p[2]))
            ij = (int(p[3]) - 1, int(p[4]) - 1)
            n = int(lines[i + 1].split()[0])
            vs = []
            for k in range(n):
                q = lines[i + 2 + k].split()
                vs.append((int(q[0]), int(q[1]), int(q[2])))
            Rs.append(r)
            ijs.append(ij)
            ns.append(n)
            vec_list.append(vs)
            i += 2 + n
        else:
            i += 1
    shift_start = np.zeros(len(vec_list) + 1, dtype=np.int64)
    shift_start[1:] = np.cumsum(ns)
    shifts_flat = np.array([v for vs in vec_list for v in vs],
                           dtype=np.int32).reshape(-1, 3)
    out = dict(
        header=header,
        use_ws_distance=header_flag(header, "use_ws_distance"),
        write_ndegen_applied=header_flag(header, "write_ndegen_applied") is True,
        header_use_ws_distance=header_flag(header, "use_ws_distance"),
        R=np.array(sorted(set(Rs)), dtype=int),
        Rblocks=np.asarray(Rs, dtype=int),
        ijs=np.asarray(ijs, dtype=int),
        n_wsvec=np.asarray(ns, dtype=np.int64),
        shift_start=shift_start,
        shifts_flat=shifts_flat,
        blocks=len(Rs),
    )
    if keep_vectors:
        out["vectors"] = vec_list
    return out


# ------------------------------------------------------------------ centres
def read_centres(path: str, num_wann: int | None = None) -> np.ndarray:
    """Read Wannier centres (lines whose first token is ``X``)."""
    tau = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.split()
            if p and p[0] == "X" and len(p) >= 4:
                tau.append([float(x) for x in p[1:4]])
    tau = np.asarray(tau, float)
    if num_wann is not None and len(tau) != num_wann:
        raise WannierInputError(
            f"{path}: found {len(tau)} 'X' centres, expected num_wann="
            f"{num_wann}")
    return tau


def atoms_as_tau(win: dict) -> np.ndarray:
    """Fallback Wannier centres: atomic positions repeated per spinor orbital.

    Only used when the *_centres.xyz file is absent; result is flagged as an
    approximation in the model metadata (never silently equal to true centres).
    """
    orbitals = {"N": 8, "V": 10, "Pb": 10}
    tau = []
    for elem, pos in win["atoms"]:
        if elem not in orbitals:
            raise WannierInputError(
                f"unknown element {elem!r} for fallback centres; provide "
                f"*_centres.xyz")
        tau += [pos] * orbitals[elem]
    tau = np.asarray(tau, float)
    if win["num_wann"] and len(tau) != win["num_wann"]:
        raise WannierInputError(
            f"fallback centres give {len(tau)} orbitals != num_wann="
            f"{win['num_wann']}")
    return tau


# ------------------------------------------------------------------ model
class WannierModel:
    """Generic Wannier90 tight-binding model H(k)."""

    kind = "wannier90"

    def __init__(self, win: dict, hr: dict, tau: np.ndarray | None = None,
                 tau_source: str = "centres", wsvec_path: str | None = None,
                 wsvec_use_ws_distance: bool | None = None):
        if win["num_wann"] != hr["num_wann"]:
            raise WannierInputError(
                f"num_wann mismatch: win={win['num_wann']} hr={hr['num_wann']}")
        self.win = win
        self.hr = hr
        self.nw = hr["num_wann"]
        self.nrpts = hr["nrpts"]
        self.ndegen = hr["ndegen"]
        self.Rlist = hr["R"]
        self.HR = hr["H"]
        self.lattice = np.asarray(win["lattice"], float)
        self.a_ref = float(np.linalg.norm(self.lattice[0]))
        if tau is not None:
            self.tau = np.asarray(tau, float)
            self.tau_source = tau_source
        else:
            self.tau = atoms_as_tau(win)
            self.tau_source = "atomic-positions (fallback, approximate)"
        # Wannier90 use_ws_distance expansion (v0.1.6-rc1)
        self.wsvec_path = None
        self.wsvec_info = None
        self.Rlist_full = None
        self.HR_full = None
        self._init_wsvec(wsvec_path, wsvec_use_ws_distance)

    def _init_wsvec(self, wsvec_path: str | None,
                    override: bool | None = None):
        """Build the expanded (Wannier90 ws_apply_ndegen) representation.

        With ``use_ws_distance=.true.`` and ``write_ndegen_applied`` absent,
        the R list in ``wannier90_hr.dat`` is *folded*; the pair-dependent
        supercell shifts in ``wannier90_wsvec.dat`` must be applied before the
        plain Fourier transform:

            H_ij(k) = sum_R (1/ndegen) sum_deg (1/ws_ndeg)
                      H_ij(R) exp(2*pi*i k.(R + shift_deg(R,i,j)))

        This is Wannier90's ws_apply_ndegen + plain Fourier on the expanded
        list.  With ``write_ndegen_applied=.true.`` the real-space files are
        already expanded and no shift is applied here.
        """
        if not wsvec_path or not os.path.exists(wsvec_path):
            return
        wv = read_wsvec(wsvec_path)
        header_flag = wv["header_use_ws_distance"]
        if override is not None:
            flag, source = bool(override), "config override"
        elif header_flag is not None:
            flag, source = header_flag, "header"
        else:
            # wannier90 writes wannier90_wsvec.dat only when
            # use_ws_distance=.true., so a token-free header means "true".
            # This is an inference, not a fallback: it is recorded below.
            flag, source = True, "inferred (wsvec.dat exists)"
        has_shifts = bool(np.any(np.asarray(wv["shifts_flat"]) != 0))
        if override is None and header_flag is False and has_shifts:
            raise WannierInputError(
                f"{wsvec_path}: header says use_ws_distance=.false. but the "
                f"file lists {int(wv['blocks'])} shift blocks; the input is "
                f"self-contradictory. Set files.wsvec_use_ws_distance to "
                f"resolve it explicitly.")
        info = dict(path=wsvec_path, header=wv["header"],
                    use_ws_distance=flag,
                    use_ws_distance_source=source,
                    header_use_ws_distance=header_flag,
                    config_override=override,
                    write_ndegen_applied=bool(wv["write_ndegen_applied"]),
                    n_blocks=int(wv["blocks"]),
                    n_nonzero_shift_blocks=int(np.sum(
                        np.any(wv["shifts_flat"] != 0, axis=1))),
                    n_multi_shift_blocks=int(np.sum(wv["n_wsvec"] > 1)))
        self.wsvec_path = wsvec_path
        if not flag or wv["write_ndegen_applied"]:
            info.update(applied_to_model=False,
                        reason=("use_ws_distance=false (%s)" % source if not
                                flag else
                                "write_ndegen_applied=true (already folded "
                                "into hr.dat)"))
            self.wsvec_info = info
            return
        ridx = {tuple(int(x) for x in r): i
                for i, r in enumerate(self.Rlist)}
        Rblocks = wv["Rblocks"]
        nblocks = wv["blocks"]
        # unique expanded vectors
        vecs = set()
        for ib in range(nblocks):
            r = Rblocks[ib]
            n0 = wv["shift_start"][ib]
            n1 = wv["shift_start"][ib + 1]
            for sh in wv["shifts_flat"][n0:n1]:
                vecs.add(tuple(int(x) for x in (r + sh)))
        Rfull = np.array(sorted(vecs), dtype=int)
        rfidx = {tuple(int(x) for x in r): j for j, r in enumerate(Rfull)}
        nrf = len(Rfull)
        Hfull = np.zeros((nrf, self.nw, self.nw), dtype=complex)
        for ib in range(nblocks):
            r = Rblocks[ib]
            ir = ridx.get(tuple(int(x) for x in r))
            if ir is None:
                raise WannierInputError(
                    f"{wsvec_path}: R {tuple(r)} not present in hr.dat")
            a, b = int(wv["ijs"][ib, 0]), int(wv["ijs"][ib, 1])
            n0 = wv["shift_start"][ib]
            n1 = wv["shift_start"][ib + 1]
            ns = int(n1 - n0)
            w = 1.0 / (float(self.ndegen[ir]) * float(ns))
            for sh in wv["shifts_flat"][n0:n1]:
                jr = rfidx[tuple(int(x) for x in (r + sh))]
                Hfull[jr, a, b] += self.HR[ir, a, b] * w
        info.update(applied_to_model=True, Rlist_full=Rfull,
                    n_expanded_R=nrf, blocks=nblocks,
                    ijs=np.asarray(wv["ijs"], int),
                    shift_start=np.asarray(wv["shift_start"], int),
                    shifts_flat=np.asarray(wv["shifts_flat"], int),
                    Rblocks=np.asarray(Rblocks, int),
                    ridx=ridx, rfidx=rfidx)
        self.wsvec_info = info
        self.Rlist_full = Rfull
        self.HR_full = Hfull

    @property
    def use_ws_distance(self) -> bool:
        return getattr(self, "HR_full", None) is not None

    def __setstate__(self, state):
        """Unpickle old cached models from before the wsvec patch safely."""
        self.__dict__.update(state)
        for k in ("wsvec_path", "wsvec_info", "Rlist_full", "HR_full"):
            self.__dict__.setdefault(k, None)

    # ---- basic
    @classmethod
    def from_files(cls, win_path: str, hr_path: str,
                   centres_path: str | None = None,
                   wsvec_path: str | None = None,
                   wsvec_use_ws_distance: bool | None = None):
        win = read_win(win_path)
        hr = read_hr(hr_path)
        tau = None
        src = "atomic-positions (fallback, approximate)"
        if centres_path and os.path.exists(centres_path):
            tau = read_centres(centres_path, hr["num_wann"])
            src = "wannier90_centres.xyz"
        return cls(win, hr, tau, src, wsvec_path=wsvec_path,
                   wsvec_use_ws_distance=wsvec_use_ws_distance)

    @classmethod
    def from_chk(cls, chk_path: str, hr_path: str,
                 centres_path: str | None = None,
                 wsvec_path: str | None = None,
                 wsvec_use_ws_distance: bool | None = None):
        """Build a model when wannier90.win is unavailable.

        The missing win metadata (lattice, k points, num_wann/num_bands,
        mp_grid) is taken from ``wannier90.chk``; consistency with the hr
        file is enforced (num_wann) and the k mesh is validated separately by
        the checkpoint H-consistency check.  The reconstructed ``win`` dict is
        flagged in ``model.win["source"]`` so validation can report it.
        """
        from .chk import read_chk
        chk = read_chk(chk_path)
        hr = read_hr(hr_path)
        if chk["num_wann"] != hr["num_wann"]:
            raise WannierInputError(
                f"num_wann mismatch: chk={chk['num_wann']} "
                f"hr={hr['num_wann']}")
        win = dict(
            num_wann=int(chk["num_wann"]), num_bands=int(chk["num_bands"]),
            lattice=np.asarray(chk["real_lattice"], float),
            kpoints=np.asarray(chk["kpt_red"], float),
            mp_grid=np.asarray(chk["mp_grid"], int),
            kpoint_weights=np.full(len(chk["kpt_red"]), 1.0
                                   / len(chk["kpt_red"])),
            atoms=[], spinors=None,
            source="reconstructed from wannier90.chk (win file missing)")
        tau = None
        src = "atomic-positions (fallback, approximate)"
        if centres_path and os.path.exists(centres_path):
            tau = read_centres(centres_path, hr["num_wann"])
            src = "wannier90_centres.xyz"
        return cls(win, hr, tau, src, wsvec_path=wsvec_path,
                   wsvec_use_ws_distance=wsvec_use_ws_distance)

    @classmethod
    def from_chk_eig(cls, chk_path: str, eig_path: str,
                     win_path: str | None = None):
        """Build a model from ``.chk`` + ``.eig`` only (no ``_hr.dat``).

        H_W(k)=v^dagger diag(eig) v is evaluated on the DFT k-mesh and Fourier
        transformed to H(R).  For a shifted mesh k_i=(n_i+s_i)/N_i, R and
        R+N_i are *not* equivalent; only shifts g=q_i*N_i with q_i*s_i integer
        leave exp(2*pi*i k.R) invariant.  The reconstruction therefore uses one
        residue R_i=0..N_i-1 per direction and only allowed image shifts,
        which reproduces H_W exactly on the mesh.  Between mesh points it is a
        Fourier interpolation of the Wannier Hamiltonian; for shifted meshes
        the coefficients are alias (Dirichlet-kernel) coefficients rather than
        physical Wannier H(R), which is stated in ``win["source"]``.
        """
        from .chk import read_chk, read_eig
        chk = read_chk(chk_path)
        eig = read_eig(eig_path)
        nk = int(chk["num_kpts"])
        nw = int(chk["num_wann"])
        nb = int(chk["num_bands"])
        if eig.shape != (nk, nb):
            raise WannierInputError(
                f"eig shape {eig.shape} != (nk,num_bands)=({nk},{nb}); the "
                f"chk and eig files are not from the same run")
        V = np.asarray(chk["v_matrix"], complex)
        lattice = np.asarray(chk["real_lattice"], float)
        kpts = np.asarray(chk["kpt_red"], float)
        ngrid = tuple(int(x) for x in chk["mp_grid"])

        Hw = np.empty((nk, nw, nw), dtype=complex)
        for ik in range(nk):
            Vk = V[ik]
            Hw[ik] = np.conj(Vk.T) @ (eig[ik].astype(complex)[:, None] * Vk)
            Hw[ik] = 0.5 * (Hw[ik] + np.conj(Hw[ik].T))

        def mesh_shift(N, vals):
            """s = k*N mod 1 for this direction (0 if N==1)."""
            if N == 1:
                return 0.0
            sv = np.mod(np.asarray(vals, float) * N, 1.0)
            sv = np.round(sv, 10)
            sv[sv > 1 - 1e-9] = 0.0
            return float(sv[0])

        shift = np.mod(kpts * np.asarray(ngrid, float), 1.0)
        shift[np.abs(shift - 1.0) < 1e-9] = 0.0
        unshifted = bool(np.all(np.abs(shift) < 1e-9))
        if unshifted:
            # Wannier90 WS/ndegen reconstruction: Hermitian for every k.
            Rlist, ndegen = ws_supercell(ngrid, lattice)
            interpolation = "wannier90_ws"
        else:
            allowed = []
            for d in range(3):
                N = ngrid[d]
                s_d = mesh_shift(N, kpts[:, d])
                q = 1
                while q < 24 and abs(round(s_d * q) - s_d * q) > 1e-7:
                    q += 1
                allowed.append(q * N)          # allowed R image shift
            reps = []
            for i in range(ngrid[0]):
                for j in range(ngrid[1]):
                    for kz in range(ngrid[2]):
                        cands = [(i + a, j + b, kz + c)
                                 for a in (-allowed[0], 0, allowed[0])
                                 for b in (-allowed[1], 0, allowed[1])
                                 for c in (-allowed[2], 0, allowed[2])]
                        reps.append(min(cands, key=lambda R: float(
                            np.linalg.norm(np.asarray(R, float) @ lattice))))
            Rlist = np.array(sorted(set(reps)), dtype=int)
            ndegen = np.ones(len(Rlist), dtype=int)
            interpolation = "shifted_alias"
        nr = len(Rlist)
        HR = np.zeros((nr, nw, nw), dtype=complex)
        for ir, R in enumerate(Rlist):
            phase = np.exp(-2j * np.pi * (kpts @ R))
            # NB: H(R) itself is *not* Hermitian (H(-R)=H(R)^dagger); only the
            # full H(k) is.  Do not symmetrise individual R blocks.
            HR[ir] = np.tensordot(phase, Hw, axes=(0, 0)) / nk

        hr = dict(num_wann=nw, nrpts=nr, R=Rlist, H=HR,
                  ndegen=np.asarray(ndegen, dtype=int),
                  header="reconstructed from chk+eig")
        win = dict(num_wann=nw, num_bands=nb, lattice=lattice,
                   kpoints=kpts, mp_grid=np.asarray(ngrid),
                   kpoint_weights=np.ones(nk) / nk, atoms=[], spinors=None,
                   source=("reconstructed from wannier90.chk+eig "
                           f"({interpolation})"))
        if win_path and os.path.exists(win_path):
            w = read_win(win_path)
            for key in ("num_bands", "lattice", "atoms", "spinors",
                        "mp_grid", "kpoints", "kpoint_weights"):
                if key in w and w[key] is not None:
                    win[key] = w[key]
            win["source"] = ("wannier90.chk+eig Fourier interpolation (+win "
                             "metadata); no hr.dat")
        tau = np.asarray(chk["centres"], float)
        obj = cls(win, hr, tau, "wannier90.chk centres")
        obj.reconstructed = True          # HR is a Fourier reconstruction
        obj.interpolation = interpolation
        return obj

    @property
    def R_cart(self) -> np.ndarray:
        """Cartesian R vectors (nrpts,3), row vectors R @ lattice."""
        if not hasattr(self, "_R_cart"):
            self._R_cart = self.Rlist @ self.lattice
        return self._R_cart

    def _Hk_raw(self, k_frac) -> np.ndarray:
        k = np.asarray(k_frac, float).reshape(3)
        if getattr(self, "HR_full", None) is not None:
            phase = np.exp(2j * np.pi * (self.Rlist_full @ k))
            return np.einsum("r,rij->ij", phase, self.HR_full)
        phase = np.exp(2j * np.pi * (self.Rlist @ k))
        return np.einsum("r,rij->ij", phase / self.ndegen, self.HR)

    def Hk(self, k_frac, check: bool = True) -> np.ndarray:
        """H(k) for a single k point (fractional coordinates).

        With ``check=True`` (default) a non-Hermitian result raises
        :class:`HermiticityError` instead of being silently averaged or fed to
        ``eigvalsh``.  ``check=False`` exists for diagnostics that need to
        *measure* the residual.
        """
        H = self._Hk_raw(k_frac)
        if check:
            assert_hermitian(H, k_frac)
        return H

    def _Hk_batch_raw(self, kpts) -> np.ndarray:
        kpts = np.atleast_2d(np.asarray(kpts, float))
        if getattr(self, "HR_full", None) is not None:
            phases = np.exp(2j * np.pi * (kpts @ self.Rlist_full.T))
            return np.einsum("kr,rij->kij", phases, self.HR_full)
        phases = np.exp(2j * np.pi * (kpts @ self.Rlist.T))   # (nk,nrpts)
        return np.einsum("kr,rij->kij", phases / self.ndegen, self.HR)

    def Hk_batch(self, kpts, check: bool = True) -> np.ndarray:
        """H(k) on many k points; ``check`` as in :meth:`Hk`."""
        Hk = self._Hk_batch_raw(kpts)
        if check:
            assert_hermitian(Hk, None, what="H(k) batch")
        return Hk


    def static_bands(self, kpts) -> np.ndarray:
        return np.linalg.eigvalsh(self.Hk_batch(kpts))

    # ---- geometry for Peierls substitution
    def displacement_tensor(self) -> np.ndarray:
        """d_ij(R) = R_cart + tau_j - tau_i, shape (nrpts,nw,nw,3).

        With use_ws_distance the expanded R list is used, so the Peierls
        phase sees the same shortest Bravais vector as the interpolation.
        """
        if getattr(self, "HR_full", None) is not None:
            if not hasattr(self, "_d_tensor_full"):
                Rc = self.Rlist_full @ self.lattice
                self._d_tensor_full = (Rc[:, None, None, :]
                                       + self.tau[None, None, :, :]
                                       - self.tau[None, :, None, :])
            return self._d_tensor_full
        if not hasattr(self, "_d_tensor"):
            Rc = self.R_cart                                 # (nrpts,3)
            self._d_tensor = (Rc[:, None, None, :]
                              + self.tau[None, None, :, :]
                              - self.tau[None, :, None, :])
        return self._d_tensor

    # ---- checks
    def hermiticity_errors(self, ktest=None) -> dict:
        """Hermiticity residuals, split by representation.

        ``HR_hermiticity_folded``  : max|H(-R) - H(R)^dag| on the *folded* R
            set from ``_hr.dat``.  ``None`` when that set is not closed under
            R -> -R, in which case the residual is *undefined* (reporting 0.0
            there used to be a false PASS).
        ``Hk_hermiticity_expanded`` : max|H(k) - H(k)^dag|, evaluated through
            :meth:`_Hk_batch_raw` so it covers ``HR_full`` (wsvec-expanded)
            models as well.  The probe mesh includes **off-mesh** k points:
            wsvec shifts are integer supercell vectors, so at any k with
            k*mp_grid integral the shift phases are all 1 and the expansion
            is invisible -- an on-mesh-only probe cannot detect a one-sided
            wsvec block set.
        """
        Ridx = {tuple(r): i for i, r in enumerate(self.Rlist)}
        err_hr = 0.0
        unpaired = []
        for i, r in enumerate(self.Rlist):
            j = Ridx.get(tuple(-r))
            if j is None:
                unpaired.append(tuple(int(x) for x in r))
                continue
            err_hr = max(err_hr, float(np.max(np.abs(self.HR[i]
                                                     - self.HR[j].conj().T))))
        folded = None if unpaired else err_hr
        if ktest is None:
            ktest = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                              [0.25, 0.25, 0.0],
                              # off-mesh probes: make the wsvec shifts visible
                              [0.137, 0.083, 0.021],
                              [0.3211, 0.1129, 0.0473]])
        ktest = np.atleast_2d(np.asarray(ktest, float))
        per_k = [float(np.max(np.abs(H - H.conj().T)))
                 for H in self._Hk_batch_raw(ktest)]
        err_k = max(per_k)
        return dict(
            HR_hermiticity_folded=folded,
            HR_hermiticity=folded,          # legacy alias, same semantics
            n_unpaired_R=len(unpaired),
            unpaired_R=unpaired[:8],
            Hk_hermiticity=err_k,           # legacy alias
            Hk_hermiticity_expanded=err_k,
            Hk_per_k=np.asarray(per_k),
            ktest=ktest,
            wsvec_expanded=getattr(self, "HR_full", None) is not None,
        )


    def ndegen_consistency(self) -> dict:
        """Check that R and -R pairs share the same degeneracy factor."""
        Ridx = {tuple(r): i for i, r in enumerate(self.Rlist)}
        bad = 0
        for i, r in enumerate(self.Rlist):
            j = Ridx.get(tuple(-r))
            if j is not None and self.ndegen[i] != self.ndegen[j]:
                bad += 1
        return dict(ndegen_mismatch=bad, ndegen_min=int(self.ndegen.min()),
                    ndegen_max=int(self.ndegen.max()))
