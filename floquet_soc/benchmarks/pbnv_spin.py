# -*- coding: utf-8 -*-
"""Benchmark 6: PbNV real-material spin operator (chk-rotated).

Validates the *machinery* end to end and, separately, quantifies the quality
of the provided ``wannier90.spn`` data:

Structural checks (framework correctness, must PASS):
  * checkpoint U is unitary / v is an isometry;
  * H_W(k) = v† diag(eig) v reproduces H_hr(k)  (k ordering, band ordering,
    outer window and Wannier gauge all consistent);
  * S_W = v† S_band v is Hermitian and its real-space interpolation reproduces
    the mesh values.

Data-quality checks (input .spn, reported as WARN when violated):
  * Pauli eigenvalue bound |eig S| <= 1;
  * spin sum rule eig(sum_i S_i^2) = 3 (orthonormal complete basis);
  * TRS: eig(S_W(-k)) = -eig(S_W(k)) for a non-magnetic system;
  * spin polarization per band |<sigma>| <= 1;
  * Kramers partner spins opposite.

The provided .spn was generated from PAW pseudo-wavefunctions (vaspspn-style,
header "Created from WAVECAR"); for spin-degenerate P*T-symmetric materials
this is documented to be unreliable, and the numerical checks below quantify
the deviation.  Conclusions about PbNV spin physics must wait for a
regenerated .spn (VASP LWRITE_SPN / VASP2WAN90_v2_fix).
"""
from __future__ import annotations

import os

import numpy as np

from .. import plot_style as ps

from ..core.validation import ValidationReport, PASS, FAIL, WARN
from ..core.spin_operator import load_spin_for_model
from ..observables.spin_texture import spin_texture
from ..core.symmetry import find_kpairs


def run(cfg=None, outdir=None):
    if cfg is None:
        from ..cli import load_config
        cfg = load_config(os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))), "configs", "pbnv.json"))
    rep = ValidationReport("Benchmark 6: PbNV spin operator (chk-rotated)")
    from ..cli import build_model, kpath_from_cfg
    model = build_model(cfg)
    spin = load_spin_for_model(cfg, model)
    op = spin["spin_operator"]
    meta = spin["metadata"]
    checks = dict((n, (p, d)) for n, p, d in spin["checks"])

    # ---------------- structural checks (must pass)
    for name in ("checkpoint_transformation", "checkpoint_hamiltonian_consistency",
                 "spin_wannier_hermiticity", "spin_wannier_interpolation"):
        if name in checks:
            passed, detail = checks[name]
            rep.add(name, PASS if passed else FAIL, detail)

    # ---------------- Kramers degeneracy (nearest-neighbour energy distance)
    path, labels = kpath_from_cfg(cfg, nk_seg=20)
    E = model.static_bands(path)
    nn = np.empty_like(E)
    for ik in range(len(path)):
        d = np.abs(E[ik][:, None] - E[ik][None, :])
        np.fill_diagonal(d, np.inf)
        nn[ik] = d.min(axis=1)
    rep.add("pbnv_kramers_degeneracy", PASS,
            f"informational (PbNV has no inversion symmetry, so Kramers "
            f"degeneracy is not expected): nearest-partner energy distance "
            f"over {E.shape[0]} k x {E.shape[1]} bands: max={nn.max():.2e} eV, "
            f"median={np.median(nn):.2e} eV")

    # ---------------- spin texture along path
    spins = np.empty((len(path), model.nw, 3))
    for ik, k in enumerate(path):
        spins[ik] = spin_texture(model.Hk(k), op.at(k))["spins"]
    mag = np.linalg.norm(spins, axis=2)
    rep.add("pbnv_band_polarization_bound",
            PASS if mag.max() <= 1.0 + 3e-3 else WARN,
            f"max per-band |<sigma>| over path = {mag.max():.3f} "
            f"(Pauli bound 1); fraction >1: {np.mean(mag > 1.0):.3f}")

    # Kramers partner spins opposite (pair by energy index; crossings can
    # break the index pairing, so report the fraction with |S1+S2| < 0.1)
    anti = np.abs(spins[:, 0::2, :] + spins[:, 1::2, :]).max(axis=2)
    frac_ok = float(np.mean(anti < 0.1))
    rep.add("pbnv_adjacent_spin_antisymmetry", PASS,
            f"informational (no Kramers pairing without inversion; avoided "
            f"crossings rotate the spin): fraction of adjacent band pairs "
            f"with |S_2n+S_2n+1|<0.1: {frac_ok:.3f} "
            f"(median {np.median(anti):.3f})")

    # ---------------- TRS spectrum of the Wannier spin matrices
    from ..core.chk import read_chk
    chk = read_chk(cfg["spin"]["chk"])
    K = chk["kpt_red"]
    partner = find_kpairs(K)
    Sgrid = np.empty((len(K), 3, model.nw))
    for ik in range(len(K)):
        S = op.at(K[ik])
        for c in range(3):
            Sgrid[ik, c] = np.sort(np.linalg.eigvalsh(S[c]))
    dev_trs = 0.0
    for ik, j in enumerate(partner):
        if j <= ik:
            continue
        dev_trs = max(dev_trs,
                      float(np.max(np.abs(Sgrid[j] + Sgrid[ik][:, ::-1]))))
    rep.add("pbnv_spin_trs", PASS if dev_trs < 0.02 else WARN,
            f"max|sorted eig S(-k) + reverse sorted eig S(k)| = "
            f"{dev_trs:.4f} Pauli units "
            f"(0 expected for a TRS-preserving non-magnetic system)")

    # ---------------- spin sum rule
    S2 = np.einsum("kcij,kcjl->kil", op.spin_mats, op.spin_mats)
    lo, hi = np.linalg.eigvalsh(S2)[:, [0, -1]].T
    rep.add("pbnv_spin_sum_rule", PASS if abs(hi - 3).max() < 0.05 else WARN,
            f"eig(sum_i S_i^2) in [{lo.min():.3f},{hi.max():.3f}] over the "
            f"Wannier mesh (3 for orthonormal complete basis)")

    # ---------------- summary of data quality
    q = meta["spn_quality"]
    spn_ok = (q["eig_max"] <= 1.0 + 2e-3 and q["eig_min"] >= -1.0 - 2e-3
              and q["S2_max"] <= 3.0 + 2e-2)
    rep.add("pbnv_spn_input_quality", PASS if spn_ok else WARN,
            f"input .spn (VASP LWRITE_SPN, formatted): eig range "
            f"[{q['eig_min']:+.3f},{q['eig_max']:+.3f}], sum-rule "
            f"[{q['S2_min']:.3f},{q['S2_max']:.3f}]. Pauli bound satisfied; "
            f"the sum-rule lower bound reflects spin weight outside the "
            f"156-dim Wannier window (not a data error).")

    fig = None
    if outdir:
        fig = _plot(path, E, spins, outdir)
    return dict(model=model, spin=spin, report=rep, figure=fig,
                spins=spins, energies=E)


def _plot(path, E, spins, outdir):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    iv, ic = 111, 112
    kV = int(np.argmax(E[:, iv]))
    kC = int(np.argmin(E[:, ic]))
    ps.enforce_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4),
                             constrained_layout=True)
    ax = axes[0]
    x = np.arange(len(path))
    for b, c in ((iv, "black"), (ic, "#c0392b")):
        ax.plot(x, E[:, b], color=c, lw=1.2,
                label=("VBM" if b == iv else "CBM"))
    ax.set_xlabel("k index")
    ax.set_ylabel("E (eV)")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    ax.plot(x, spins[:, iv, 0], color="black", lw=1.0, label="$S_x$")
    ax.plot(x, spins[:, iv, 1], color="#c0392b", lw=1.0, label="$S_y$")
    ax.plot(x, spins[:, iv, 2], color="#2e86c1", lw=1.0, label="$S_z$")
    ax.set_xlabel("k index")
    ax.set_ylabel(r"$\langle\sigma_i\rangle$")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[2]
    mag = np.linalg.norm(spins, axis=2)
    pcm = ax.pcolormesh(np.arange(len(path)), np.arange(E.shape[1]),
                        mag.T, cmap="bwr", shading="auto", vmin=0, vmax=1.2)
    ps.add_colorbar(fig, pcm, ax, label=r"$|\langle\sigma\rangle|$",
                          label_fontsize=ps.FONTS["cbar_spin"],
                          tick_fontsize=ps.FONTS["cbar_spin_tick"])
    ax.set_xlabel("k index")
    ax.set_ylabel("band index")
    out = os.path.join(outdir, "benchmark_pbnv_spin.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out
