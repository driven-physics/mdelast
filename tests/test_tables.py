"""Tabulated potentials: interpolation option and properties of the quintic tables."""

from __future__ import annotations

import contextlib
import io
import math
import os
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = Path(
    os.environ.get("MDELAST_BUILD_DIR", ROOT / "build" / "Release")
).resolve()
HAS_EXTENSION = BUILD_DIR.is_dir() and any(BUILD_DIR.glob("md_driver_module*.so"))


def load(interpolation=None):
    """ZnCu tables; the C++ constructors print to stdout, which is silenced."""
    import get_potential

    with contextlib.redirect_stdout(io.StringIO()):
        if interpolation is None:
            return get_potential.get_EAM_alloy_file("ZnCu.eam.alloy")
        return get_potential.get_EAM_alloy_file("ZnCu.eam.alloy", interpolation=interpolation)


@unittest.skipUnless(
    HAS_EXTENSION,
    "compiled md_driver_module not found; build first or set MDELAST_BUILD_DIR",
)
class TableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        sys.path.insert(0, str(BUILD_DIR))
        sys.path.insert(0, str(ROOT / "src" / "python"))
        import md_driver_module

        cls.A = md_driver_module.Units.angstrom
        cls.quintic = load()
        cls.bspline = load("bspline")
        # physical range of the pair functions and densities: 1.8 A to just below the cutoff
        cls.r = np.linspace(1.8, 6.7, 997) * cls.A

    def test_default_and_option(self) -> None:
        import get_potential

        _, pots, eams = self.quintic
        self.assertEqual(pots[0].interpolation, "quintic")
        self.assertEqual(eams[0].interpolation, "quintic")
        _, pots, eams = self.bspline
        self.assertEqual(pots[0].interpolation, "bspline")
        self.assertEqual(eams[0].interpolation, "bspline")
        with self.assertRaises(ValueError):
            load("cubic")
        saved = get_potential.DEFAULT_INTERPOLATION
        try:
            get_potential.DEFAULT_INTERPOLATION = "bspline"
            self.assertEqual(load()[1][0].interpolation, "bspline")
        finally:
            get_potential.DEFAULT_INTERPOLATION = saved

    def _pair(self, tables, k=0):
        v = [tables[1][k].calculate(float(x)) for x in self.r]
        return np.array([[p.V, p.dVdr, p.d2Vdr2] for p in v])

    def test_quintic_derivatives_are_consistent(self) -> None:
        # f' and f'' are the exact derivatives of the interpolant: central differences agree
        # to the truncation error of the difference quotient
        pot = self.quintic[1][2]
        d = 1e-5 * self.A
        for x in self.r[::10]:
            p, m, c = (pot.calculate(float(x + s)) for s in (d, -d, 0.0))
            self.assertAlmostEqual((p.V - m.V) / (2 * d) / c.dVdr, 1.0, delta=1e-7)
            self.assertAlmostEqual((p.dVdr - m.dVdr) / (2 * d), c.d2Vdr2,
                                   delta=1e-6 * abs(c.d2Vdr2) + 1e-6 * self._scale2)

    @property
    def _scale2(self) -> float:
        return float(np.max(np.abs(self._pair(self.quintic, 2)[:, 2])))

    def test_continuity_across_nodes(self) -> None:
        pot = self.quintic[1][1]
        h = 3.4017008504252126e-03 * self.A          # grid spacing of ZnCu.eam.alloy
        scale = np.max(np.abs(self._pair(self.quintic, 1)), axis=0)
        for i in range(600, 1950, 37):
            x = i * h
            lo, hi = pot.calculate(float(np.nextafter(x, 0))), pot.calculate(float(np.nextafter(x, 1)))
            jumps = np.abs([hi.V - lo.V, hi.dVdr - lo.dVdr, hi.d2Vdr2 - lo.d2Vdr2]) / scale
            self.assertTrue(np.all(jumps < 1e-12), (i, jumps))

    def test_quintic_and_bspline_agree(self) -> None:
        for k in range(3):
            q, b = self._pair(self.quintic, k), self._pair(self.bspline, k)
            rel = np.max(np.abs(q - b), axis=0) / np.max(np.abs(q), axis=0)
            self.assertTrue(np.all(rel < [1e-8, 1e-6, 1e-4]), (k, rel))
        for e in range(2):
            eq, eb = self.quintic[2][e], self.bspline[2][e]
            for f in ("rho", "drhodr", "d2rhodr2"):
                a = np.array([getattr(eq, f)(float(x)) for x in self.r])
                b = np.array([getattr(eb, f)(float(x)) for x in self.r])
                self.assertLess(np.max(np.abs(a - b)) / np.max(np.abs(a)), 1e-4, f)

    def test_nan_argument_returns_nan(self) -> None:
        # an invalid argument must not read outside the table
        self.assertTrue(math.isnan(self.quintic[1][0].calculate(float("nan")).V))
        self.assertTrue(math.isnan(self.quintic[2][0].F(float("nan"))))


if __name__ == "__main__":
    unittest.main()
