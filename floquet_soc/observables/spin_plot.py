# -*- coding: utf-8 -*-
"""2D spin-texture plotting (reference style, BZ-clipped).

Style:
  * S_z background: ``tricontourf`` with ``bwr``, clipped to the 2D BZ,
  * black arrows of uniform length (direction = in-plane spin S_xy),
  only inside the BZ and clipped to the BZ boundary; points with
  |S_xy|/|S_z| below ``dot_ratio`` are drawn as black dots,
  * black BZ boundary, equal aspect, no ticks/spines,
  * colorbar equal in height to the BZ panel.

k points are fractional reciprocal coordinates (``WannierModel.Hk`` convention)
and converted to 1/Angstrom for the axes.  ``kplane`` supports non-orthogonal
lattices (fcc Pt: physical k_z=0 plane = [[1,0,1],[0,1,1]]).
"""
from __future__ import annotations

import numpy as np

from .. import plot_style as ps


# ------------------------------------------------------------------ lattice
def reciprocal_rows(lattice) -> np.ndarray:
    A = np.asarray(lattice, float)
    return 2.0 * np.pi * np.linalg.inv(A).T


def kfrac_to_cart(kfrac, lattice) -> np.ndarray:
    return np.asarray(kfrac, float) @ reciprocal_rows(lattice)


def plane_matrices(lattice, kplane=None):
    """(E, M): E=(2,3) plane basis rows (fractional), M = E @ B."""
    B = reciprocal_rows(lattice)
    if kplane is None:
        E = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    else:
        E = np.atleast_2d(np.asarray(kplane, float))
        if E.shape != (2, 3):
            raise ValueError("kplane must be 2x3 fractional vectors")
    return E, E @ B


def plane_xy(kpts_frac, kplane=None):
    E = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]) if kplane is None \
        else np.atleast_2d(np.asarray(kplane, float))
    return np.asarray(kpts_frac, float) @ E.T @ np.linalg.inv(E @ E.T)


def plane_kpts(kplane, n, kmax):
    """k = x*e1 + y*e2, with n/kmax scalars or (nx,ny)/(kxmax,kymax)."""
    E = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]) if kplane is None \
        else np.atleast_2d(np.asarray(kplane, float))
    if np.ndim(n) == 0:
        nx = ny = int(n)
    else:
        nx, ny = int(n[0]), int(n[1])
    if np.ndim(kmax) == 0:
        kxmax = kymax = float(kmax)
    else:
        kxmax, kymax = float(kmax[0]), float(kmax[1])
    xs = np.linspace(-kxmax, kxmax, nx)
    ys = np.linspace(-kymax, kymax, ny)
    return np.array([[x * E[0] + y * E[1]]
                     for x in xs for y in ys]).reshape(-1, 3)


# ------------------------------------------------------------------ BZ polygon
def _clip_polygon(poly, a, b, c, tol=1e-12):
    out = []
    n = len(poly)
    for i in range(n):
        p = np.asarray(poly[i], float)
        q = np.asarray(poly[(i + 1) % n], float)
        fp = a * p[0] + b * p[1] - c
        fq = a * q[0] + b * q[1] - c
        pin, qin = fp <= tol, fq <= tol
        if pin:
            out.append(tuple(p))
        if pin != qin:
            t = fp / (fp - fq)
            out.append(tuple(p + t * (q - p)))
    return out


def bz_polygon_plane(lattice, kplane=None, nG=3):
    """BZ cross-section in the plane (x,y) coordinates (fractional)."""
    A = np.asarray(lattice, float)
    if A.shape != (3, 3) or abs(np.linalg.det(A)) < 1e-12:
        return np.empty((0, 2))
    B = reciprocal_rows(A)
    E, M = plane_matrices(lattice, kplane)
    poly = [(-2.5, -2.5), (2.5, -2.5), (2.5, 2.5), (-2.5, 2.5)]
    for g1 in range(-nG, nG + 1):
        for g2 in range(-nG, nG + 1):
            for g3 in range(-nG, nG + 1):
                if g1 == g2 == g3 == 0:
                    continue
                Gc = np.array([g1, g2, g3], float) @ B
                rhs = 0.5 * float(Gc @ Gc)
                a = M @ Gc
                if abs(a[0]) < 1e-14 and abs(a[1]) < 1e-14:
                    continue
                poly = _clip_polygon(poly, a[0], a[1], rhs)
                if not poly:
                    return np.empty((0, 2))
    pts = np.asarray(poly, float)
    out = []
    for p in pts:
        if not out or np.linalg.norm(p - out[-1]) > 1e-8:
            out.append(p)
    return np.asarray(out, float)


def bz_polygon_2d(lattice, nG=3):
    """Backward-compatible k_z=0 polygon in fractional (k1,k2)."""
    return bz_polygon_plane(lattice, None, nG=nG)


# ------------------------------------------------------------------ drawing core
def _align_colorbar(ax, cbar, ylo, yhi):
    """Force the colorbar vertical extent to match [ylo, yhi] in data coords."""
    fig = ax.figure
    try:
        fig.canvas.draw()
        inv = ax.transAxes.inverted()
        y0 = inv.transform(ax.transData.transform((0.0, ylo)))[1]
        y1 = inv.transform(ax.transData.transform((0.0, yhi)))[1]
        ap = ax.get_position()
        cp = cbar.ax.get_position()
        cbar.ax.set_position([cp.x0, ap.y0 + y0 * ap.height, cp.width,
                              (y1 - y0) * ap.height])
    except Exception:                                     # noqa: BLE001
        pass


def _panel(ax, kpts_frac, S, lattice=None, kplane=None, vmin=-1.0, vmax=1.0,
           cmap="bwr", label="", colorbar_label=r"$S_z$", quiver_scale=16.0,
           dot_ratio=0.25,
           show_bz=True):
    S = np.asarray(S, float)
    kf = np.atleast_2d(np.asarray(kpts_frac, float))
    pxy = None
    path = None
    if lattice is not None:
        E, M = plane_matrices(lattice, kplane)
        xy = plane_xy(kf, kplane)
        cart = xy @ M
        Bmat = M[:, :2]
        poly_xy = bz_polygon_plane(lattice, kplane)
        if len(poly_xy) >= 3:
            from matplotlib.path import Path
            pxy = poly_xy @ Bmat
            path = Path(pxy)
    else:
        cart = kf
    x, y = cart[:, 0], cart[:, 1]
    z = S[:, 2]

    if path is not None:
        from matplotlib.tri import Triangulation
        tri = Triangulation(x, y)
        cx = x[tri.triangles].mean(axis=1)
        cy = y[tri.triangles].mean(axis=1)
        tri.set_mask(~path.contains_points(np.column_stack([cx, cy])))
        cf = ax.tricontourf(tri, z, levels=np.linspace(vmin, vmax, 121),
                            cmap=cmap, vmin=vmin, vmax=vmax, zorder=1)
        try:
            from matplotlib.patches import PathPatch
            cf.set_clip_path(PathPatch(path, transform=ax.transData,
                                       facecolor="none"))
        except Exception:                                 # noqa: BLE001
            pass
    else:
        cf = ax.tricontourf(x, y, z, levels=np.linspace(vmin, vmax, 121),
                            cmap=cmap, vmin=vmin, vmax=vmax, zorder=1)

    u, v = S[:, 0], S[:, 1]
    sz_abs = np.abs(S[:, 2])
    sel = np.ones(len(x), bool)
    if path is not None:
        sel = path.contains_points(np.column_stack([x, y]))
    uu, vv, zz = u[sel], v[sel], sz_abs[sel]
    s_xy = np.hypot(uu, vv)
    if s_xy.max() <= 1e-12:
        # fully out-of-plane: draw dots everywhere (the |S_xy|/|S_z| rule
        # has set the in-plane components to zero)
        ax.scatter(x[sel], y[sel], s=6.5, c="black", marker="o",
                   linewidths=0.0, zorder=3)
    if s_xy.max() > 1e-12:
        # reference style (Desktop/plot_spintexture.py): in-plane spin is
        # normalised to a unit vector, all arrows have the same physical size;
        # |S_xy|/|Sz| < dot_ratio is drawn as a black dot instead.
        ratio = s_xy / np.maximum(zz, 1e-12)
        dots = (ratio < float(dot_ratio)) | (s_xy <= 1e-12)
        clip = None
        if path is not None:
            try:
                from matplotlib.patches import PathPatch
                clip = PathPatch(path, transform=ax.transData,
                                 facecolor="none")
            except Exception:                             # noqa: BLE001
                clip = None
        if dots.any():
            scat = ax.scatter(x[sel][dots], y[sel][dots], s=6.5, c="black",
                              marker="o", linewidths=0.0, zorder=3)
            if clip is not None:
                scat.set_clip_path(clip)
        arr = ~dots
        if arr.any():
            ux = uu[arr] / s_xy[arr]
            uy = vv[arr] / s_xy[arr]
            # reference Desktop/plot_spintexture.py: width=0.008,
            # scale=16, pivot="mid", scale_units="inches", default heads
            q = ax.quiver(x[sel][arr], y[sel][arr], ux, uy,
                          color="black", width=0.008,
                          scale=float(quiver_scale), pivot="mid",
                          scale_units="inches", zorder=3)
            if clip is not None:
                q.set_clip_path(clip)

    if pxy is not None:
        from matplotlib.patches import Polygon
        ax.add_patch(Polygon(pxy, closed=True, fill=False, edgecolor="black",
                             linewidth=1.8, zorder=4))
        xlo, xhi = np.min(pxy[:, 0]), np.max(pxy[:, 0])
        ylo, yhi = np.min(pxy[:, 1]), np.max(pxy[:, 1])
        mx = 0.02 * max(xhi - xlo, yhi - ylo)
        ax.set_xlim(xlo - mx, xhi + mx)
        ax.set_ylim(ylo - mx, yhi + mx)
    else:
        ax.set_xlim(np.min(x), np.max(x))
        ax.set_ylim(np.min(y), np.max(y))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = ax.figure.colorbar(cf, ax=ax, fraction=0.050, pad=0.025,
                              shrink=1.0, aspect=22)
    cbar.set_ticks(np.linspace(vmin, vmax, 5))
    cbar.set_label(colorbar_label, fontsize=ps.FONTS["cbar_spin"],
                   labelpad=4)
    cbar.outline.set_linewidth(1.0)
    cbar.ax.tick_params(labelsize=ps.FONTS["cbar_spin_tick"], width=1.0)
    if pxy is not None:
        ax._fsoc_cbar = cbar
        ax._fsoc_y = (float(np.min(pxy[:, 1])), float(np.max(pxy[:, 1])))
    return cf


# ------------------------------------------------------------------ public API
def plot_spin_texture_2d(kpts_frac, spins, path, title="", lattice=None,
                         kplane=None, vmin=-1.0, vmax=1.0, level_count=121,
                         dpi=300, figsize=(3.2, 3.0), show_bz=True,
                         colorbar_label=r"$S_z$", quiver_scale=16.0,
                         dot_ratio=0.25):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps.enforce_style()
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    _panel(ax, kpts_frac, spins, lattice=lattice, kplane=kplane,
           vmin=vmin, vmax=vmax, label="", colorbar_label=colorbar_label,
           quiver_scale=quiver_scale, dot_ratio=dot_ratio,
           show_bz=show_bz)
    fig.tight_layout()
    if hasattr(ax, "_fsoc_cbar"):
        _align_colorbar(ax, ax._fsoc_cbar, *ax._fsoc_y)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_spin_texture_panels(panels, path, suptitle="", vmin=-1.0, vmax=1.0,
                             dpi=300, figsize=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps.enforce_style()
    npan = len(panels)
    if figsize is None:
        figsize = (3.2 * npan + 0.35, 3.05)
    fig, axes = plt.subplots(1, npan, figsize=figsize, dpi=dpi,
                             constrained_layout=True)
    if npan == 1:
        axes = [axes]
    for ax, p in zip(axes, panels):
        _panel(ax, p["kpts"], p["spins"], lattice=p.get("lattice"),
               kplane=p.get("kplane"), vmin=p.get("vmin", vmin),
               vmax=p.get("vmax", vmax), cmap=p.get("cmap", "bwr"),
               label=p.get("label", ""),
               colorbar_label=p.get("colorbar_label", r"$S_z$"),
               quiver_scale=p.get("quiver_scale", 16.0),
               dot_ratio=p.get("dot_ratio", 0.25),
               show_bz=p.get("show_bz", True))
    fig.canvas.draw()
    try:
        fig.set_layout_engine("none")
    except Exception:                                     # noqa: BLE001
        pass
    for ax in axes:
        if hasattr(ax, "_fsoc_cbar"):
            _align_colorbar(ax, ax._fsoc_cbar, *ax._fsoc_y)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path
