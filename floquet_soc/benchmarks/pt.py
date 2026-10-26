# -*- coding: utf-8 -*-
"""Benchmark 7: fcc Pt -- small spinor SOC material (18 Wannier functions).

The data come from the public Wannier90 tutorial repository
(``wannier-developers/wannier-tutorials``, 2022 Trieste school, spin-Hall
tutorial) and contain ``.win/.chk/.eig/.spn/.wout`` but *no* ``_hr.dat``.
This benchmark therefore also validates the checkpoint+eigenvalue input path:

  * the model is rebuilt from ``.chk``+``.eig`` (H(R) by Fourier transform on
    the 4x4x4 DFT mesh) -- exact on the mesh by construction;
  * the chk-rotated spin operator uses a *proper* spin-matrix file (QE
    pw2wannier90), so Pauli bound and sum rule should be satisfied at the
    ~1e-3 level (contrast with the pseudo-wavefunction PbNV .spn);
  * band interpolation between mesh points is compared with the wannier90
    reference bands (``Pt_band.dat``): a documented limitation of the
    hr-less reconstruction.

No orbital counts / electron counts are hard-coded in the framework; the
values below come from the material config.
"""
from __future__ import annotations

import os

import numpy as np

from .. import plot_style as ps

from ..core.validation import ValidationReport, PASS, FAIL, WARN
from ..core.spin_operator import load_spin_for_model
from ..observables.spin_texture import spin_texture


def _default_cfg():
    from ..cli import load_config
    return load_config(os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))), "configs", "pt.json"))


def _band_reference(cfg):
    """Load wannier90 interpolated bands (Pt_band.dat + Pt_band.kpt)."""
    files = cfg.get("files", {})
    band = files.get("bands_ref")
    kpt = files.get("bands_ref_kpt")
    if not (band and os.path.exists(band) and kpt and os.path.exists(kpt)):
        return None
    blocks, cur = [], []
    for ln in open(band, encoding="utf-8"):
        if ln.strip():
            p = ln.split()
            cur.append([float(p[0]), float(p[1])])
        elif cur:
            blocks.append(np.array(cur))
            cur = []
    if cur:
        blocks.append(np.array(cur))
    E = np.array([b[:, 1] for b in blocks]).T
    lines = [ln.split() for ln in open(kpt, encoding="utf-8") if ln.strip()]
    path = np.array([[float(x) for x in ln[:3]] for ln in lines[1:]], float)
    return np.sort(E, axis=1), path


def run(cfg=None, outdir=None):
    if cfg is None:
        cfg = _default_cfg()
    rep = ValidationReport("Benchmark 7: fcc Pt (spinor SOC, chk+eig input)")
    from ..cli import build_model
    model = build_model(cfg)
    win_src = str(model.win.get("source", "wannier90.win"))
    rep.add("pt_model_source", PASS,
            f"num_wann={model.nw}, num_bands={model.win['num_bands']}, "
            f"kpoints={len(model.win['kpoints'])}, "
            f"mp_grid={list(model.win['mp_grid'])}; source={win_src}")
    if getattr(model, "reconstructed", False) or "reconstructed" in win_src:
        interp = getattr(model, "interpolation", "unknown")
        rep.add("pt_hr_reconstruction", WARN,
                f"no official wannier90_hr.dat available: H(R) rebuilt from "
                f"chk+eig with the Wannier90 WS/ndegen convention "
                f"({interp}); Hermitian off-grid and validated against the "
                f"official PbNV hr.dat, but Pt off-grid accuracy is limited "
                f"by the 4x4x4 mesh. Backend label: grid-consistent "
                f"WS/ndegen reconstruction; official-hr comparison pending.")

    spin = load_spin_for_model(cfg, model)
    op = spin["spin_operator"]
    meta = spin["metadata"]
    checks = dict((n, (p, d)) for n, p, d in spin["checks"])
    for name in ("checkpoint_transformation",
                 "checkpoint_hamiltonian_consistency",
                 "spin_wannier_hermiticity",
                 "spin_wannier_hermiticity_offgrid",
                 "spin_wannier_interpolation",
                 "spin_wannier_trs"):
        if name in checks:
            passed, detail = checks[name]
            rep.add(name, PASS if passed else FAIL, detail)
    q = meta["spn_quality"]
    good_spn = (q["eig_max"] <= 1.0 + 2e-2 and q["eig_min"] >= -1.0 - 2e-2
                and q["S2_max"] <= 3.0 + 5e-2)
    from ..observables.kramers import kramers_summary
    rng = np.random.default_rng(20261022)
    ktest = rng.uniform(-0.5, 0.5, size=(200, 3))
    Hk_off = model.Hk_batch(ktest)
    herm_off = float(np.max(np.abs(Hk_off - np.conj(
        Hk_off.transpose(0, 2, 1)))))
    rep.add("pt_hermiticity_offgrid", PASS if herm_off < 1e-10 else FAIL,
            f"max|H(k)-H(k)^dag| on 200 off-grid k = {herm_off:.2e} eV "
            f"(tol 1e-10 eV; Hermiticity must be checked before Kramers)")
    kmesh = np.asarray(model.win["kpoints"], float)
    ks_mesh = kramers_summary(model.static_bands(kmesh), kmesh)
    rep.add("pt_kramers_mesh",
            PASS if ks_mesh["delta_max"] < 1e-6 else FAIL,
            f"delta_PT=max_n|E_{{2n-1}}-E_{{2n}}| on the original DFT mesh: "
            f"max {ks_mesh['delta_max']*1000:.3e} meV, median "
            f"{ks_mesh['delta_median']*1000:.3e} meV (tol 1e-6 eV = 1e-3 meV, "
            f"from DFT eigenvalue precision 2e-10 eV)")
    krand = rng.uniform(-0.5, 0.5, size=(400, 3))
    ks_off = kramers_summary(model.static_bands(krand), krand)
    med_mev = ks_off["delta_median"] * 1000
    max_mev = ks_off["delta_max"] * 1000
    rep.add("pt_kramers_offgrid", WARN,
            f"off-grid Kramers residual: median {med_mev:.3f} meV, "
            f"max {max_mev:.2f} meV at k="
            f"{np.round(ks_off['worst_k'], 4).tolist()}. WARN, not PASS: for "
            f"a strictly PT-symmetric static crystal the exact residual is "
            f"zero; these meV values are numerical/interpolation error of the "
            f"4x4x4 WS/ndegen Hamiltonian, not physical. The error must be "
            f"tracked against an official Wannier90 hr.dat before any PASS.")
    rep.add("pt_quantitative_status", WARN,
            "Pt is a software-interface and error-detection benchmark only "
            "at present: the official Wannier90 hr.dat is unavailable, the "
            "Pt_band.dat reference is not reproduced (index-aligned median "
            "54 meV / max 1.14 eV), and the off-grid Kramers residual is a "
            "numerical artefact. Do NOT use Pt for quantitative SOC/Floquet "
            "material conclusions until the official-interpolation "
            "comparison passes.")
    rep.add("pt_soc_splitting_not_applicable", PASS,
            "PT-symmetric non-magnetic static Pt: SOC spin splitting is NOT "
            "APPLICABLE (Kramers pairs must be degenerate). `soc` reports "
            "delta_PT and inter-doublet gaps, no spin-paired splitting.")
    rep.add("pt_pst_gauge_invariance", WARN,
            "NOT_APPLICABLE: single-band PST is gauge dependent for a "
            "Kramers-degenerate PT system; the random-U(2) test shows the "
            "PST axis is arbitrary (dot product down to 0.21).")
    rep.add("pt_spn_quality", PASS if good_spn else WARN,
            f"input .spn (QE pw2wannier90): eig range "
            f"[{q['eig_min']:+.4f},{q['eig_max']:+.4f}], sum-rule "
            f"[{q['S2_min']:.3f},{q['S2_max']:.3f}] (Pauli bound 1, sum 3)")

    # ---- spin texture + PST on the small mesh
    km = dict(cfg.get("kmesh") or dict(n=25, kmax=0.55))
    from ..observables.spin_plot import plane_kpts
    kplane = cfg.get("kplane")
    if "dk" in km or "n" not in km:
        from ..cli import resolve_kmesh
        n_km, kmax_km = resolve_kmesh(cfg, "kmesh", default_n=25,
                                      default_kmax=0.55)
        kk = plane_kpts(kplane, n_km, kmax_km)
    else:
        kk = plane_kpts(kplane, int(km["n"]), km["kmax"])
    spins = np.empty((len(kk), model.nw, 3))
    for ik, k in enumerate(kk):
        spins[ik] = spin_texture(model.Hk(k), op.at(k))["spins"]
    mag = np.linalg.norm(spins, axis=2)
    ne = cfg.get("num_electrons")
    ib = int(ne) - 1 if ne else model.nw // 2
    ib = max(0, min(ib, model.nw - 1))
    rep.add("pt_spin_texture_finite",
            PASS if np.all(np.isfinite(spins)) else FAIL,
            f"{len(kk)} k-points x {model.nw} bands; figure band {ib}")
    rep.add("pt_band_polarization_bound",
            PASS if mag.max() <= 1.0 + 5e-2 else WARN,
            f"max per-band |<sigma>| = {mag.max():.3f} (Pauli bound 1)")

    # ---- compare interpolation with the wannier90 reference bands
    ref = _band_reference(cfg)
    if ref is not None:
        E_ref, path = ref
        npath = min(len(path), E_ref.shape[0])
        E_our = np.sort(model.static_bands(path[:npath]), axis=1)
        E_refs = np.sort(E_ref, axis=1)
        nb = min(E_our.shape[1], E_refs.shape[1])
        d = np.abs(E_our[:, :nb] - E_refs[:, :nb])
        med, mx = float(np.median(d) * 1e3), float(d.max() * 1e3)
        rep.add("pt_band_interpolation_vs_wannier90", WARN,
                f"index-aligned deviation over {npath} path points: "
                f"median {med:.1f} meV, max {mx:.1f} meV. WARN, not PASS: "
                f"the official Pt hr.dat is unavailable and Pt_band.dat may "
                f"have been produced with different Wannier settings "
                f"(use_ws_distance); this comparison cannot validate the "
                f"off-grid Hamiltonian.")

    fig = None
    if outdir:
        fig = _plot(kk, spins, mag, ib, outdir,
                    lattice=getattr(model, "lattice", None),
                    kplane=cfg.get("kplane"))
    return dict(model=model, spin=spin, report=rep, figure=fig,
                spins=spins, kpts=kk, band=ib)


def _plot(kk, spins, mag, ib, outdir, lattice=None, kplane=None):
    """PbNV-style single-panel spin texture (S_z bwr + arrows + BZ + dots)."""
    import os

    from ..observables.spin_plot import plot_spin_texture_2d
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "benchmark_pt.png")
    plot_spin_texture_2d(
        kk, spins[:, ib, :], out,
        title="Benchmark 7: fcc Pt spinor SOC model "
              f"(18 Wannier functions) - spin texture, band {ib}",
        lattice=lattice, kplane=kplane)
    return out
