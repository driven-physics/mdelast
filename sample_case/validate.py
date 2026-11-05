#!/usr/bin/env python3
"""Validate a completed manuscript sample run against the shipped tolerances."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACTUAL = PACKAGE_ROOT / "sample_case_output" / "result.json"
DEFAULT_REFERENCE = Path(__file__).resolve().parent / "expected_output" / "result.json"


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str


EXACT_FIELDS = (
    "case",
    "cycles_thermalization",
    "production_cycles",
    "placement_seed",
    "npt_velocity_seed",
    "nvt_velocity_seed",
)

NUMERIC_RULES = (
    ("npt.box_length_A", "box_length_A"),
    ("nvt.box_length_A", "box_length_A"),
    ("nvt.C11_GPa", "C11_GPa"),
    ("nvt.C12_GPa", "C12_GPa"),
    ("nvt.C44_GPa", "C44_GPa"),
    ("npt.temperature_kinetic_K", "temperature_K"),
    ("npt.temperature_configurational_reciprocal_mean_K", "temperature_K"),
    ("npt.temperature_configurational_sample_mean_K", "temperature_K"),
    ("nvt.temperature_kinetic_K", "temperature_K"),
    ("nvt.temperature_configurational_reciprocal_mean_K", "temperature_K"),
    ("nvt.temperature_configurational_sample_mean_K", "temperature_K"),
)


def nested_value(data: dict[str, Any], dotted_path: str) -> Any:
    value: Any = data
    for part in dotted_path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(dotted_path)
        value = value[part]
    return value


def validate_results(
    actual: dict[str, Any], reference: dict[str, Any]
) -> list[Check]:
    checks: list[Check] = []

    for field in EXACT_FIELDS:
        try:
            actual_value = nested_value(actual, field)
            reference_value = nested_value(reference, field)
        except KeyError:
            checks.append(Check(field, False, "required field is missing"))
            continue
        passed = actual_value == reference_value
        checks.append(
            Check(
                field,
                passed,
                f"actual={actual_value!r}, expected={reference_value!r}",
            )
        )

    try:
        tolerances = nested_value(reference, "acceptance_tolerance")
    except KeyError:
        tolerances = {}

    for field, tolerance_key in NUMERIC_RULES:
        try:
            actual_value = float(nested_value(actual, field))
            reference_value = float(nested_value(reference, field))
            tolerance = float(tolerances[tolerance_key])
        except (KeyError, TypeError, ValueError):
            checks.append(
                Check(field, False, f"field or tolerance {tolerance_key!r} is missing/invalid")
            )
            continue

        difference = abs(actual_value - reference_value)
        passed = (
            math.isfinite(actual_value)
            and math.isfinite(reference_value)
            and math.isfinite(tolerance)
            and tolerance >= 0.0
            and difference <= tolerance
        )
        checks.append(
            Check(
                field,
                passed,
                f"actual={actual_value:.10g}, expected={reference_value:.10g}, "
                f"|delta|={difference:.4g}, tolerance={tolerance:.4g}",
            )
        )

    return checks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "actual",
        nargs="?",
        type=Path,
        default=DEFAULT_ACTUAL,
        help=f"result JSON to validate (default: {DEFAULT_ACTUAL})",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_REFERENCE,
        help=f"reference JSON (default: {DEFAULT_REFERENCE})",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise SystemExit(f"file not found: {path}") from error
    except json.JSONDecodeError as error:
        raise SystemExit(f"invalid JSON in {path}: {error}") from error
    if not isinstance(value, dict):
        raise SystemExit(f"expected a JSON object in {path}")
    return value


def main() -> int:
    args = parse_args()
    actual = load_json(args.actual.resolve())
    reference = load_json(args.reference.resolve())
    checks = validate_results(actual, reference)

    for check in checks:
        stream = sys.stdout if check.passed else sys.stderr
        status = "PASS" if check.passed else "FAIL"
        print(f"{status:4} {check.name}: {check.detail}", file=stream)

    failed = [check for check in checks if not check.passed]
    if failed:
        print(
            f"Validation failed: {len(failed)} of {len(checks)} checks failed.",
            file=sys.stderr,
        )
        return 1

    print(f"Validation passed: all {len(checks)} checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
