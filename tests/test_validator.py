"""Unit tests for the manuscript result validator."""

from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

from sample_case.validate import validate_results


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_PATH = ROOT / "sample_case" / "expected_output" / "result.json"


class ValidatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))

    def test_reference_validates_itself(self) -> None:
        checks = validate_results(self.reference, self.reference)
        self.assertTrue(checks)
        self.assertTrue(all(check.passed for check in checks))

    def test_value_outside_tolerance_fails(self) -> None:
        actual = copy.deepcopy(self.reference)
        actual["nvt"]["C11_GPa"] += (
            self.reference["acceptance_tolerance"]["C11_GPa"] + 0.01
        )
        checks = validate_results(actual, self.reference)
        c11_check = next(check for check in checks if check.name == "nvt.C11_GPa")
        self.assertFalse(c11_check.passed)

    def test_value_at_tolerance_passes(self) -> None:
        actual = copy.deepcopy(self.reference)
        actual["nvt"]["C44_GPa"] += self.reference["acceptance_tolerance"]["C44_GPa"]
        checks = validate_results(actual, self.reference)
        c44_check = next(check for check in checks if check.name == "nvt.C44_GPa")
        self.assertTrue(c44_check.passed)

    def test_metadata_mismatch_fails(self) -> None:
        actual = copy.deepcopy(self.reference)
        actual["placement_seed"] = 7
        checks = validate_results(actual, self.reference)
        seed_check = next(check for check in checks if check.name == "placement_seed")
        self.assertFalse(seed_check.passed)

    def test_missing_required_field_fails(self) -> None:
        actual = copy.deepcopy(self.reference)
        del actual["nvt"]["temperature_kinetic_K"]
        checks = validate_results(actual, self.reference)
        temperature_check = next(
            check for check in checks if check.name == "nvt.temperature_kinetic_K"
        )
        self.assertFalse(temperature_check.passed)


if __name__ == "__main__":
    unittest.main()
