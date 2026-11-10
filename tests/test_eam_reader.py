"""SETFL reader: pair blocks of three or more elements reach the right driver slots.

A SETFL file lists the pair functions as 11, 21, 22, 31, 32, 33, ...; MDDriver
looks up the pair (i, j), i <= j, in slot j(j+1)/2 + i.  For two elements every
ordering of the blocks gives the same slots, so the Cu-Zn tests cannot detect a
wrong order.  These tests use a synthetic file whose pair functions are
distinct constants, so a block in the wrong slot changes the result.
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = Path(
    os.environ.get("MDELAST_BUILD_DIR", ROOT / "build" / "Release")
).resolve()
HAS_EXTENSION = BUILD_DIR.is_dir() and any(BUILD_DIR.glob("md_driver_module*.so"))

NRHO, DRHO = 32, 0.1      # embedding grid
NR, DR = 64, 0.1          # radial grid [A]; the last node, 6.3 A, is the cutoff
RHO_CONST = 1.0e-3        # density function, far below the embedding-grid limit


def pair_value(i, j):
    """Pair energy [eV] of elements i <= j: 11, 21, 22, 31, 32, 33, ..."""
    return 10.0 * (j + 1) + (i + 1)


def write_setfl(path, n_elements):
    """Synthetic SETFL file: zero embedding, constant density, constant pair functions."""
    names = ["Aa", "Bb", "Cc", "Dd"][:n_elements]
    r = DR * np.arange(NR)
    out = ["synthetic potential for the reader test", "pair function ij = 10 j + i [eV]", "",
           f"{n_elements} " + " ".join(names),
           f"{NRHO} {DRHO:.6f} {NR} {DR:.6f} {DR * (NR - 1):.6f}"]
    for k in range(n_elements):
        out.append(f"{k + 1} {10.0 * (k + 1):.4f} 3.6000 fcc")
        out += ["0.0"] * NRHO                    # F(rho) = 0
        out += [f"{RHO_CONST:.6e}"] * NR         # rho(r)
    for j in range(n_elements):                 # blocks 11, 21, 22, 31, 32, 33, ...
        for i in range(j + 1):
            out += [f"{pair_value(i, j) * x:.10e}" for x in r]   # r * phi(r)
    Path(path).write_text("\n".join(out) + "\n")


@contextlib.contextmanager
def quiet():
    """Silence the C++ constructor messages (file descriptor 1) and the 0/0 at r = 0
    that the reader overwrites."""
    sys.stdout.flush()
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with np.errstate(divide="ignore", invalid="ignore"):
                yield
    finally:
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


@unittest.skipUnless(
    HAS_EXTENSION,
    "compiled md_driver_module not found; build first or set MDELAST_BUILD_DIR",
)
class SetflReaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        sys.path.insert(0, str(BUILD_DIR))
        sys.path.insert(0, str(ROOT / "src" / "python"))
        import get_potential
        import md_driver_module

        cls.gp = get_potential
        cls.m = md_driver_module
        cls.U = md_driver_module.Units
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def read(self, n_elements, interpolation):
        path = Path(self.tmp.name) / f"synthetic{n_elements}.eam.alloy"
        if not path.exists():
            write_setfl(path, n_elements)
        with quiet():
            # an absolute path bypasses the potentials/ directory of the package
            return self.gp.get_EAM_alloy_file(str(path), interpolation=interpolation)

    def test_pair_blocks_reach_the_driver_slots(self) -> None:
        r = 3.0 * self.U.angstrom
        for n in (3, 4):
            for interpolation in ("quintic", "bspline"):
                elements, potentials, eams = self.read(n, interpolation)
                self.assertEqual([e[0] for e in elements], ["Aa", "Bb", "Cc", "Dd"][:n])
                self.assertEqual(len(potentials), n * (n + 1) // 2)
                self.assertEqual(len(eams), n)
                for j in range(n):
                    for i in range(j + 1):
                        with self.subTest(n=n, interpolation=interpolation, pair=(i, j)):
                            slot = j * (j + 1) // 2 + i
                            v = potentials[slot].calculate(r).V / self.U.ev
                            self.assertAlmostEqual(v, pair_value(i, j), places=8)

    def test_three_element_energy(self) -> None:
        # five atoms of types 0, 1, 1, 2, 2, all within the cutoff of one another
        elements, potentials, eams = self.read(3, "quintic")
        types = [0, 1, 1, 2, 2]
        pos = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, 2.6, 0.0],
                        [0.0, 0.0, 2.7], [2.0, 2.1, 2.2]])            # [A]
        box = 30.0                                                    # [A], no periodic images in range
        expected = sum(pair_value(*sorted((types[a], types[b])))
                       for a in range(5) for b in range(a + 1, 5))    # [eV]
        m, U = self.m, self.U
        P = [m.Vector3d(*(p * U.angstrom)) for p in pos]
        V = [m.Vector3d(0.0, 0.0, 0.0) for _ in pos]
        th = m.Thermostat(300.0, [1e-30], [0.0], 3.0 * len(pos))
        baro = m.Barostat(0.0, 1e-30, 0.0, m.Thermostat(300.0, [1e-30], [0.0], 1.0))
        rg = [(k + 0.5) * 0.01 * U.angstrom for k in range(10)]
        mass = [e[2] / U.Navogadro for e in elements]
        with quiet():
            d = m.MDDriver(P, V, m.PotentialVector(potentials), m.EAMVector(eams),
                           m.MDDriver.NVE, m.MDDriver.THERMOSTAT_KINETIC, th, baro,
                           m.Box(np.diag(np.full(3, box * U.angstrom))), rg,
                           3.0 * len(pos), mass, types, 1e-15, False, 1.0, 1e3)
        energy = d.epot_atom * len(pos) / U.ev
        self.assertAlmostEqual(energy, expected, places=6)
        # the block order 11, 12, 13, 22, 23, 33 would give a different energy
        wrong = {(0, 0): 11, (0, 1): 21, (0, 2): 22, (1, 1): 31, (1, 2): 32, (2, 2): 33}
        mistaken = sum(wrong[tuple(sorted((types[a], types[b])))]
                       for a in range(5) for b in range(a + 1, 5))
        self.assertNotAlmostEqual(expected, mistaken, places=3)


if __name__ == "__main__":
    unittest.main()
