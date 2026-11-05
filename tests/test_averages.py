"""Unit tests for statistical helpers used by the MD output layer."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

from averages import gaussian, get_mean_cum, get_mean_std, get_mean_std_cum

# NumPy 2 renamed np.trapz to np.trapezoid; NumPy 1.x has only np.trapz.
trapezoid = getattr(np, "trapezoid", None) or np.trapz


class AverageTests(unittest.TestCase):
    def test_mean_and_population_standard_deviation(self) -> None:
        values = np.array([1.0, 2.0, 3.0, 4.0])
        mean, standard_deviation = get_mean_std(values)
        self.assertAlmostEqual(float(mean), 2.5)
        self.assertAlmostEqual(float(standard_deviation), float(np.std(values)))

    def test_cumulative_mean(self) -> None:
        values = np.array([2.0, 4.0, 8.0])
        np.testing.assert_allclose(get_mean_cum(values), [2.0, 3.0, 14.0 / 3.0])

    def test_cumulative_statistics_end_at_full_statistics(self) -> None:
        values = np.array([1.5, -2.0, 4.0, 6.5])
        means, standard_deviations = get_mean_std_cum(values)
        mean, standard_deviation = get_mean_std(values)
        self.assertAlmostEqual(float(means[-1]), float(mean))
        self.assertAlmostEqual(float(standard_deviations[-1]), float(standard_deviation))

    def test_gaussian_is_normalized_numerically(self) -> None:
        x = np.linspace(-8.0, 8.0, 100_001)
        density = gaussian(x, 0.0, 1.0)
        self.assertAlmostEqual(float(trapezoid(density, x)), 1.0, places=6)


if __name__ == "__main__":
    unittest.main()
