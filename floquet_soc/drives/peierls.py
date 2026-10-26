# -*- coding: utf-8 -*-
"""Periodic driving interfaces.

Stage A implements the general vector-potential interface and a Peierls
substitution driver; the Floquet solvers that consume ``H(k,t)`` live in
``floquet_soc.solvers`` (stage B).

Vector potential (length gauge, spatial uniform)::

    A(t) = A1 cos(omega t) e1 + A2 cos(omega t + phi) e2

with amplitude related to the dimensionless drive strength by
``gamma = e A0 a_ref / hbar`` (a_ref = |a1|, the first lattice constant).
The electric field is ``E(t) = -dA/dt``.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.special import jv

EV = 1.0          # energies and hbar*omega are in eV
HBAR = 1.0        # hbar = 1 in eV units (omega is an angular frequency in eV)


class DriveConfigError(RuntimeError):
    pass


class TruncationWarning(UserWarning):
    """Raised when the Bessel series for the Peierls phase is cut too short."""


# De-duplication: a k-scan calls this once per k point, and the warning text
# embeds the offending value; warn only when a new (larger) g_max is reached.
_warned_gmax = 0.0


class Drive:
    """Configurable monochromatic drive.

    Parameters
    ----------
    hw : float
        Photon energy hbar*omega in eV.
    gamma : float
        Dimensionless drive amplitude gamma = e A0 a_ref / hbar.
    a_ref : float
        Reference length (Angstrom), normally |a1| of the lattice.
    polarization : {"linear", "circular", "elliptical"}
    e1, e2 : array-like
        Unit vectors spanning the polarization plane (Cartesian, i.e. as the
        Cartesian displacement vectors d_ij).  Defaults: x and y.
    phi : float
        Phase difference between the two components (radians).
    ellipticity : float
        A2/A1 ratio (0 = linear, 1 = circular with phi=+-pi/2).
    sign : +1 or -1
        Handedness for circular drive (Phi = sign*pi/2).
    """

    def __init__(self, hw: float, gamma: float, a_ref: float = 1.0,
                 polarization: str = "circular", e1=(1.0, 0.0, 0.0),
                 e2=(0.0, 1.0, 0.0), phi: float = 0.5 * np.pi,
                 ellipticity: float | None = None, sign: int = +1,
                 amplitude: float | None = None):
        if hw <= 0:
            raise DriveConfigError("photon energy hbar*omega must be positive")
        self.hw = float(hw)
        self.omega = self.hw / HBAR
        self.gamma = float(gamma)
        self.a_ref = float(a_ref)
        if amplitude is not None:
            self.A0 = float(amplitude)
        else:
            # gamma = e A0 a_ref / hbar ; here e=hbar=1 (length gauge in eV)
            self.A0 = self.gamma / self.a_ref
        e1 = np.asarray(e1, float)
        e2 = np.asarray(e2, float)
        self.e1 = e1 / np.linalg.norm(e1)
        e2p = e2 - np.dot(e2, self.e1) * self.e1
        if np.linalg.norm(e2p) < 1e-12:
            raise DriveConfigError("e1 and e2 are parallel")
        self.e2 = e2p / np.linalg.norm(e2p)
        pol = polarization.lower()
        if pol == "linear":
            ratio, phi = 0.0, 0.0
        elif pol == "circular":
            ratio, phi = 1.0, sign * 0.5 * np.pi
        elif pol == "elliptical":
            ratio = 1.0 if ellipticity is None else float(ellipticity)
            phi = float(phi)
        else:
            raise DriveConfigError(
                f"unknown polarization {polarization!r} (linear/circular/"
                f"elliptical)")
        self.polarization = pol
        self.ratio = ratio
        self.phi = float(phi)
        self.sign = int(sign)

    @property
    def period(self) -> float:
        return 2 * np.pi / self.omega

    def A(self, t):
        t = np.asarray(t, float)
        c1 = np.cos(self.omega * t)
        c2 = np.cos(self.omega * t + self.phi)
        shape = t.shape + (3,)
        return (self.A0 * (self.ratio * c2[..., None] * self.e2[None, :]
                           + c1[..., None] * self.e1[None, :])).reshape(shape)

    def E_field(self, t):
        """E(t) = -dA/dt (eV/Angstrom if A in 1/Angstrom)."""
        A = self.A(t)
        return -np.gradient(A, t, axis=0) if np.ndim(t) else np.zeros(3)

    def as_dict(self):
        return dict(hw=self.hw, gamma=self.gamma, A0=self.A0,
                    a_ref=self.a_ref, polarization=self.polarization,
                    ratio=self.ratio, phi=self.phi, sign=self.sign,
                    e1=self.e1.tolist(), e2=self.e2.tolist(),
                    period=self.period)

    @classmethod
    def from_config(cls, drv: dict, a_ref: float = 1.0):
        return cls(hw=drv["hw"], gamma=drv.get("gamma", 0.0), a_ref=a_ref,
                   polarization=drv.get("polarization", "circular"),
                   e1=drv.get("e1", (1, 0, 0)),
                   e2=drv.get("e2", (0, 1, 0)),
                   phi=drv.get("phi", 0.5 * np.pi),
                   ellipticity=drv.get("ellipticity"),
                   sign=drv.get("sign", +1),
                   amplitude=drv.get("amplitude"))


# ------------------------------------------------------------------ Peierls
def _hr_arrays(model):
    """Real-space arrays consistent with the interpolation convention.

    With Wannier90 use_ws_distance (and write_ndegen_applied absent) the
    model exposes the expanded R list / weights in ``Rlist_full`` /
    ``HR_full``; the Peierls phase then uses the same shortest Bravais
    vectors as H(k).  Otherwise the folded convention is used.
    """
    if getattr(model, "HR_full", None) is not None:
        return (np.asarray(model.Rlist_full, int),
                np.asarray(model.HR_full, complex), None,
                model.displacement_tensor())
    return (np.asarray(model.Rlist, int), np.asarray(model.HR, complex),
            np.asarray(model.ndegen, int), model.displacement_tensor())


def peierls_Hk(model, k_frac, A_vec) -> np.ndarray:
    """H(k,t) = sum_R H(R)/ndegen * exp(i[k.R + A.d(R)]) (Wannier model).

    ``A_vec`` is the vector potential in 1/Angstrom (Peierls phase A.d_ij).
    For A=0 this reduces exactly to model.Hk(k).  With Wannier90
    use_ws_distance the expanded R list of ``wannier90_wsvec.dat`` is used
    (same as ``model.Hk``), so the drive sees the same hopping vectors.
    """
    k = np.asarray(k_frac, float).reshape(3)
    A = np.asarray(A_vec, float).reshape(3)
    Rlist, HR, ndegen, d = _hr_arrays(model)
    phase_R = np.exp(2j * np.pi * (Rlist @ k))
    phase_A = np.exp(1j * np.einsum("rnmc,c->rnm", d, A))
    pref = phase_R if ndegen is None else phase_R / ndegen
    return np.einsum("r,rnm,rnm->nm", pref, phase_A, HR)


# ------------------------------------------------------------------ dispatch
def hkt(model, drive, k, t):
    """Time-dependent Hamiltonian H(k,t) for the generic Floquet solvers.

    * Wannier model: Peierls two-centre substitution (exact same H(k,t) as the
      validated PbNV engine);
    * analytic continuum models: minimal coupling H0(k + A(t)) in the same
      length gauge (k and A both in 1/Angstrom).  For lattice models always
      use the Wannier branch.
    """
    A = np.asarray(drive.A(np.array([float(t)]))).reshape(3)
    if hasattr(model, "HR"):
        return peierls_Hk(model, k, A)
    return model.Hk(np.asarray(k, float).reshape(3) + A)


# ------------------------------------------------------------------ Bessel form
def peierls_hr_components(model, drive, pmax=3):
    """Photon-resolved **real-space** Hamiltonian ``H_p(R)``.

    With ``A(t) = cos(wt) P + sin(wt) Q`` and ``A(t).d = g cos(wt - theta)``,

        exp(i A(t).d) = sum_p i^p J_p(g) exp(-i p theta) exp(i p w t)

    so each hopping ``H(R)`` acquires a photon ladder

        H_p(R) = H(R) i^p J_p(g_R) exp(-i p theta_R),

    and the k-space components are ``H_p(k) = sum_R H_p(R) e^{ikR}/ndegen``.
    This is the object a downstream code needs: writing ``H_p(R)`` in
    Wannier90 ``_hr.dat`` format gives an exact light-dressed real-space
    Hamiltonian (photon layers) without any in-memory Floquet diagonalization.

    Returns dict with ``Rlist``, ``ndegen`` (None for the already-expanded
    wsvec representation), ``Hp_R`` {p: (nrpts,nw,nw)}, ``gmax`` and
    ``truncation_residual``.
    """
    Rlist, HR, ndegen, d = _hr_arrays(model)
    A0, ratio, phi = drive.A0, drive.ratio, drive.phi
    P = A0 * (drive.e1 + ratio * np.cos(phi) * drive.e2)
    Q = -A0 * ratio * np.sin(phi) * drive.e2
    dd = np.einsum("rnmc,c->rnm", d, P)
    bb = np.einsum("rnmc,c->rnm", d, Q)
    g = np.hypot(dd, bb)
    theta = np.arctan2(bb, dd)
    # Truncation guard: exp(i A.d) = sum_p i^p J_p(g) e^{i p omega t}, cut at
    # |p|=pmax.  A Hermiticity check is blind to this (H(t) stays Hermitian),
    # so estimate the largest discarded Bessel weight explicitly.  The tail
    # |J_p(g)| decays like (e g / 2p)^p, so a short probe past pmax suffices.
    gmax = float(np.max(np.abs(g))) if g.size else 0.0
    tail = 0.0
    if gmax > 1e-6:
        tail = max(abs(jv(p, gmax)) for p in range(pmax + 1, pmax + 9))
        global _warned_gmax
        if tail > 1e-3 and gmax > _warned_gmax:
            _warned_gmax = gmax
            warnings.warn(
                f"Peierls Bessel truncation at pmax={pmax} drops "
                f"|J_p(g)| up to {tail:.3e} (g_max={gmax:.3f}): the driven "
                f"H(k,t) is unreliable. Raise the photon order (m_sp / pmax) "
                f"or lower the drive strength.", TruncationWarning)
    Hp_R = {}
    for p in range(-pmax, pmax + 1):
        coef = (1j ** p) * jv(p, g) * np.exp(-1j * p * theta)
        Hp_R[p] = np.asarray(HR) * coef
    return dict(Rlist=np.asarray(Rlist, int), ndegen=ndegen, Hp_R=Hp_R,
                gmax=gmax, truncation_residual=float(tail), pmax=int(pmax))


def peierls_fourier_components(model, drive, k, pmax=3):
    """Fourier components H_p(k) of the Peierls H(k,t).

    A(t) = cos(omega t) P + sin(omega t) Q  with
        P = A0 (e1 + ratio cos(phi) e2),   Q = -A0 ratio sin(phi) e2,
    so A(t).d = g cos(omega t - theta) with g=hypot(P.d, Q.d),
    theta=atan2(Q.d, P.d), and

        exp(i A(t).d) = sum_p i^p J_p(g) exp(-i p theta) exp(i p omega t).

    Returns dict p -> H_p (nw,nw), p = -pmax..pmax.
    """
    hr = peierls_hr_components(model, drive, pmax=pmax)
    phase_R = np.exp(2j * np.pi * (hr["Rlist"] @ np.asarray(k, float)))
    if hr["ndegen"] is not None:
        phase_R = phase_R / hr["ndegen"]
    return {p: np.einsum("r,rij->ij", phase_R, H)
            for p, H in hr["Hp_R"].items()}


def peierls_sample(model, drive, k, t, Hp=None, pmax=3):
    """Fast H(k,t) from Bessel Fourier components."""
    if Hp is None:
        Hp = peierls_fourier_components(model, drive, k, pmax=pmax)
    H = np.zeros((model.nw, model.nw), dtype=complex)
    for p, H_p in Hp.items():
        H += H_p * np.exp(1j * p * drive.omega * t)
    return H
