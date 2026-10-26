# -*- coding: utf-8 -*-
"""Light-dressed Wannier export: H_p(R) files, round-trip and Peierls checks.

``export_dressed`` writes one Wannier90 ``_hr.dat`` per photon order p, whose
Fourier sum is H_p(k).  Downstream codes can read them directly; this test
locks the file format and the physical convention:

    H_p(k) = sum_R H_p(R)/ndegen exp(2*pi*i k.R)
    H(k,t) = sum_p H_p(k) exp(-i p omega t)
"""
from __future__ import annotations

import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import WannierModel, read_hr, write_hr  # noqa: E402
from floquet_soc.drives.peierls import (Drive, hkt,                   # noqa: E402
                                        peierls_hr_components,
                                        peierls_fourier_components)
from floquet_soc.core.model_factory import load_config, build_model   # noqa: E402


def _toy():
    win = dict(num_wann=2, num_bands=2, spinors=True,
               lattice=np.eye(3) * 4.0, mp_grid=np.array([2, 2, 1]),
               kpoints=np.zeros((1, 3)))
    h0 = np.array([[0.0, 0.25], [0.25, 0.3]], complex)
    h1 = np.array([[0.1, 0.2j], [-0.2j, -0.1]], complex)
    hr = dict(num_wann=2, nrpts=3,
              R=np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0)]),
              H=np.asarray([h0, h1, h1.conj().T], complex),
              ndegen=np.ones(3, int))
    return WannierModel(win, hr, tau=np.zeros((2, 3)))


class TestWriteHrRoundTrip(unittest.TestCase):

    def test_write_then_read(self):
        m = _toy()
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "hr.dat")
            write_hr(p, m.nw, m.Rlist, m.ndegen, m.HR)
            back = read_hr(p)
        self.assertEqual(back["num_wann"], 2)
        self.assertEqual(back["nrpts"], 3)
        np.testing.assert_allclose(back["R"], m.Rlist)
        np.testing.assert_allclose(back["H"], m.HR, atol=1e-12)


class TestDressedComponents(unittest.TestCase):

    def test_hr_components_fourier_matches_kspace(self):
        m = _toy()
        drive = Drive(hw=1.0, gamma=0.2, a_ref=4.0)
        hr = peierls_hr_components(m, drive, pmax=3)
        k = np.array([0.13, 0.07, 0.0])
        Hp_k = peierls_fourier_components(m, drive, k, pmax=3)
        phase = np.exp(2j * np.pi * (m.Rlist @ k)) / m.ndegen
        for p, Hp_R in hr["Hp_R"].items():
            np.testing.assert_allclose(
                np.einsum("r,rij->ij", phase, Hp_R), Hp_k[p], atol=1e-12)

    def test_time_zero_sum_reproduces_peierls_hkt(self):
        m = _toy()
        drive = Drive(hw=1.0, gamma=0.2, a_ref=4.0)
        k = np.array([0.13, 0.07, 0.0])
        Hp = peierls_fourier_components(m, drive, k, pmax=6)
        Hsum = sum(Hp.values())
        H0 = hkt(m, drive, k, 0.0)
        np.testing.assert_allclose(Hsum, H0, atol=1e-6)

    def test_zero_drive_leaves_only_p0(self):
        m = _toy()
        drive = Drive(hw=1.0, gamma=0.0, a_ref=4.0)
        hr = peierls_hr_components(m, drive, pmax=3)
        np.testing.assert_allclose(hr["Hp_R"][0], m.HR, atol=1e-12)
        for p in (-1, 1, -2, 2, -3, 3):
            self.assertLess(float(np.max(np.abs(hr["Hp_R"][p]))), 1e-12)


class TestExportOnShippedConfig(unittest.TestCase):

    def test_export_pt_writes_files(self):
        cfg = load_config(os.path.join(HERE, "configs", "pt.json"))
        model = build_model(cfg)
        from floquet_soc.drives.peierls import peierls_hr_components
        drive = Drive(hw=float(cfg["drive"]["hw"]),
                      gamma=float(cfg["drive"]["gamma"]),
                      a_ref=float(model.a_ref))
        hr = peierls_hr_components(model, drive, pmax=2)
        with tempfile.TemporaryDirectory() as d:
            for p, H in hr["Hp_R"].items():
                path = os.path.join(d, f"hr_p{p:+d}.dat")
                write_hr(path, model.nw, hr["Rlist"],
                         hr["ndegen"] if hr["ndegen"] is not None
                         else np.ones(len(hr["Rlist"]), int), H)
                back = read_hr(path)
                self.assertEqual(back["num_wann"], model.nw)
                self.assertEqual(back["nrpts"], len(hr["Rlist"]))


if __name__ == "__main__":
    unittest.main()
