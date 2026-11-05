"""Short end-to-end smoke tests for every shipped example."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py


ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = Path(
    os.environ.get("MDELAST_BUILD_DIR", ROOT / "build" / "Release")
).resolve()
HAS_EXTENSION = BUILD_DIR.is_dir() and any(BUILD_DIR.glob("md_driver_module*.so"))


@unittest.skipUnless(
    HAS_EXTENSION,
    "compiled md_driver_module not found; build first or set MDELAST_BUILD_DIR",
)
class ExampleSmokeTests(unittest.TestCase):
    CASES = (
        ("ar_lj_gas", "examples/ar_lj_gas/run.py"),
        ("cu_300k", "examples/cu_300k/run.py"),
        ("zn06_300k", "examples/zn06_300k/run.py"),
    )

    def test_all_examples_complete_short_run(self) -> None:
        for case, script in self.CASES:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temp:
                output_dir = Path(temp) / "output"
                command = [
                    sys.executable,
                    str(ROOT / script),
                    "--build-dir",
                    str(BUILD_DIR),
                    "--output-dir",
                    str(output_dir),
                    "--cycles-thermalization",
                    "10",
                    "--production-cycles",
                    "20",
                ]
                completed = subprocess.run(
                    command,
                    cwd=ROOT,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    timeout=180,
                    check=False,
                )
                self.assertEqual(
                    completed.returncode,
                    0,
                    f"{case} failed:\n{completed.stdout[-4000:]}",
                )

                result_path = output_dir / "result.json"
                h5_path = output_dir / "series.h5"
                self.assertTrue(result_path.is_file())
                self.assertTrue(h5_path.is_file())
                self.assertTrue((output_dir / "dump" / "0.xyz").is_file())

                result = json.loads(result_path.read_text(encoding="utf-8"))
                self.assertEqual(result["case"], case)
                self.assertEqual(result["cycles_thermalization"], 10)
                self.assertEqual(result["production_cycles"], 20)
                for field in ("temperature_kinetic_K", "pressure_GPa"):
                    self.assertTrue(math.isfinite(float(result[field])), field)

                with h5py.File(h5_path, "r") as h5:
                    self.assertIn("series/time", h5)
                    self.assertIn("thermal_results", h5)
                    self.assertGreater(len(h5["series/time"]), 0)

    def test_sample_case_completes_short_npt_to_nvt_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            output_dir = Path(temp) / "sample"
            command = [
                sys.executable,
                str(ROOT / "sample_case" / "run.py"),
                "--build-dir",
                str(BUILD_DIR),
                "--output-dir",
                str(output_dir),
                "--cycles-thermalization",
                "10",
                "--production-cycles",
                "20",
            ]
            completed = subprocess.run(
                command,
                cwd=ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=180,
                check=False,
            )
            self.assertEqual(
                completed.returncode,
                0,
                f"sample case failed:\n{completed.stdout[-4000:]}",
            )
            result = json.loads(
                (output_dir / "result.json").read_text(encoding="utf-8")
            )
            self.assertIn("npt", result)
            self.assertIn("nvt", result)
            for stage in ("npt", "nvt"):
                self.assertTrue((output_dir / stage / "series.h5").is_file())
            for field in ("C11_GPa", "C12_GPa", "C44_GPa"):
                self.assertTrue(math.isfinite(float(result["nvt"][field])), field)


if __name__ == "__main__":
    unittest.main()
