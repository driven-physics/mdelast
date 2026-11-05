"""Integrity tests for the shipped manuscript reference artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "sample_case" / "expected_output"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ReferenceDataTests(unittest.TestCase):
    def test_reference_metadata_is_complete(self) -> None:
        result = json.loads((EXPECTED / "result.json").read_text(encoding="utf-8"))
        self.assertEqual(result["placement_seed"], 42)
        self.assertEqual(result["npt_velocity_seed"], 0)
        self.assertEqual(result["nvt_velocity_seed"], 0)
        self.assertEqual(result["cycles_thermalization"], 50_000)
        self.assertEqual(result["production_cycles"], 100_000)

    def test_decimated_series_are_monotonic_and_complete(self) -> None:
        for name in ("npt_series_decimated.csv", "nvt_series_decimated.csv"):
            with (EXPECTED / name).open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            steps = [int(row["step"]) for row in rows]
            times = [float(row["time_ps"]) for row in rows]
            self.assertEqual(len(rows), 301, name)
            self.assertEqual(steps[0], 0, name)
            self.assertEqual(steps[-1], 149_999, name)
            self.assertTrue(all(b > a for a, b in zip(steps, steps[1:])), name)
            self.assertTrue(all(b > a for a, b in zip(times, times[1:])), name)

    def test_nvt_csv_final_elastic_constants_match_json(self) -> None:
        result = json.loads((EXPECTED / "result.json").read_text(encoding="utf-8"))
        with (EXPECTED / "nvt_series_decimated.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            final_row = list(csv.DictReader(handle))[-1]
        for component in ("C11", "C12", "C44"):
            csv_value = float(final_row[f"{component}_running_mean_GPa"])
            json_value = float(result["nvt"][f"{component}_GPa"])
            self.assertTrue(math.isclose(csv_value, json_value, abs_tol=1e-12))

    def test_potential_checksums_match_notices(self) -> None:
        self.assertEqual(
            sha256(ROOT / "potentials" / "ZnCu.eam.alloy"),
            "bdb0497544fa28eff8ed81fc799da9bf8d6cf8b4e907658f40b18c13652ee2ce",
        )
        self.assertEqual(
            sha256(ROOT / "potentials" / "Cu01.eam.alloy.txt"),
            "5e374260448eea5c26de59432d75e702c1778cc6999d2527a789e5387ce9dc9e",
        )


if __name__ == "__main__":
    unittest.main()
