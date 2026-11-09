#!/usr/bin/env python3
"""Replace the expected values of an MD benchmark by those of a new reference run.

Usage (from the repository root):
    python3 benchmarks/update_expected.py zn06_300K    <reference_output_dir> --label "..."
    python3 benchmarks/update_expected.py l12_cuzn3_5K <reference_output_dir> --label "..."

<reference_output_dir> is the output of the benchmark's run script (npt/series.h5 and
nvt/series.h5). Every "expected" value in benchmarks/expected/<benchmark>.json is replaced by
the value of the reference run, extracted exactly as verify.py extracts it. Tolerances, bounds,
SE_run and s_vel are left unchanged, the previous value is kept as "expected_previous", and
the label is stored under source.expected_values. Run lengths must match the reference schedule.
verify.py reads the JSON directly; benchmarks/expected/EXPECTED.md is updated separately, from
the table this script prints.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bench_lib  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("benchmark", choices=("zn06_300K", "l12_cuzn3_5K"))
    p.add_argument("reference_dir", type=Path)
    p.add_argument("--label", required=True, help="provenance of the reference run (code version, options)")
    a = p.parse_args()

    path = HERE / "expected" / f"{a.benchmark}.json"
    exp = bench_lib.load_json(path)
    rows = []
    for stage, st in exp["stages"].items():
        s = bench_lib.summarize_run(a.reference_dir / stage / "series.h5")
        if (s["n_steps"], s["n_thermalization"]) != (st["n_steps"], st["n_thermalization"]):
            raise SystemExit(f"{stage}: reference run length differs from the benchmark schedule")
        s["econs_drift_abs_per_ns"] = abs(s["econs_drift_per_ns"])
        for q, spec in st["quantities"].items():
            if "expected" not in spec:
                continue
            new = float(s[q])
            spec.setdefault("expected_previous", spec["expected"])
            rows.append((f"{stage}.{q}", spec["expected_previous"], new, spec.get("tol")))
            spec["expected"] = new
    exp.setdefault("source", {})["expected_values"] = a.label
    path.write_text(json.dumps(exp, indent=1) + "\n", encoding="utf-8")

    print(f"{'quantity':24s} {'previous':>14s} {'new':>14s} {'diff/tol':>9s}")
    for q, old, new, tol in rows:
        r = (new - old) / tol if tol else float("nan")
        print(f"{q:24s} {old:14.8g} {new:14.8g} {r:9.2f}")
    print(f"written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
