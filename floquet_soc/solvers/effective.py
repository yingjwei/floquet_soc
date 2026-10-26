# -*- coding: utf-8 -*-
"""High-frequency effective Hamiltonian (stage B).

Van Vleck / Schrieffer-Wolff expansion in the photon energy:

    H_eff(ω) = H_0 + (1/ħω) Σ_{n>0} [H_n, H_{-n}]/n + O(1/ω²)

with H_n the Fourier components of H(k,t), H_n = (1/T)∫ H(t) e^{-inωt} dt.
The formula is only valid when ħω is large compared with the bandwidth and
the drive amplitude (no multiphoton resonances).  ``HighFrequencyNotValid`` is
raised when the user asks for the expansion in a parameter regime where the
criterion is clearly violated; the convergence benchmark in
``benchmarks/floquet_rashba.py`` quantifies the error against the exact
evolution solver.
"""
from __future__ import annotations

import numpy as np

from .sambe import fourier_components


class HighFrequencyNotValid(RuntimeError):
    pass


def effective_hamiltonian(model, drive, k, order=2, n_harm=3, n_step=400,
                          check_validity=True):
    """Return H_eff (matrix), order ∈ {1 (H_0 only), 2 (first 1/ω correction)}."""
    if order > 2:
        raise ValueError(
            f"effective_hamiltonian: order={order} is not supported; "
            f"only the 0th (order=1) and first 1/omega (order=2) corrections "
            f"are implemented")
    Hp = fourier_components(model, drive, k, n_harm=n_harm, n_step=n_step)
    H0 = Hp[0]
    if order <= 1:
        return H0
    corr = np.zeros_like(H0)
    for n in range(1, n_harm + 1):
        if n in Hp and -n in Hp:
            corr += (Hp[n] @ Hp[-n] - Hp[-n] @ Hp[n]) / (n * drive.hw)
    H_eff = H0 + corr
    if check_validity:
        # The 1/ω expansion breaks at *resonances*: when ħω is not much larger
        # than a level spacing the drive can bridge.  The old scale used
        # max|eig(H0)| -- the full bandwidth, dominated by deep valence bands,
        # which made this check a permanent no-op.  Use the smallest real
        # spacing as the resonance scale; on exact degeneracies fall back to
        # the bandwidth.
        eig = np.linalg.eigvalsh(H0)
        spacings = np.diff(np.sort(eig))
        positive = spacings[spacings > 1e-9]
        if positive.size:
            scale = float(np.min(positive))
        else:
            scale = float(np.max(np.abs(eig)))
        if np.max(np.abs(corr)) > 0.5 * max(scale, 1e-12):
            raise HighFrequencyNotValid(
                f"1/omega correction ({np.max(np.abs(corr)):.3g} eV) is not "
                f"small compared with the level-spacing scale ({scale:.3g} eV); "
                f"the high-frequency expansion is not reliable here")
    return H_eff


def effective_bands(model, drive, k, order=2, **kw):
    return np.linalg.eigvalsh(effective_hamiltonian(model, drive, k,
                                                    order=order, **kw))


def effective_split_error(model, drive, k, exact_eps, **kw):
    """Compare exact quasienergies with H_eff eigenvalues (same branch).

    Returns dict with max|exact - eff| assuming both sets lie in the first
    Floquet zone (true for the benchmark parameters used here).
    """
    e = np.sort(np.asarray(exact_eps, float))
    h = np.sort(effective_bands(model, drive, k, **kw))
    n = min(len(e), len(h))
    return dict(max_abs=float(np.max(np.abs(e[:n] - h[:n]))),
                exact=e, effective=h)
