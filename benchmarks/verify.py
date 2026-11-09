#!/usr/bin/env python3
"""Compare benchmark outputs with the expected values; print PASS/FAIL per quantity.

Usage (from the repository root):
    python3 benchmarks/verify.py smoke        <output_dir>   # benchmarks/smoke/run_smoke.py
    python3 benchmarks/verify.py zn06_300K    <output_dir>   # sample_case/run.py output
    python3 benchmarks/verify.py l12_cuzn3_5K <output_dir>   # benchmarks/l12_cuzn3_5K/run.py
    python3 benchmarks/verify.py all          <root>         # every <root>/<benchmark> that exists

For zn06_300K and l12_cuzn3_5K, <output_dir> must contain npt/series.h5 and
nvt/series.h5 (the layout written by the run scripts). Exit status is 0 only if
no check fails. Only NumPy and h5py are needed (no compiled extension).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bench_lib  # noqa: E402

EXPECTED = HERE / "expected"
BENCHMARKS = ("smoke", "zn06_300K", "l12_cuzn3_5K")


class Report:
    def __init__(self, title):
        self.title = title
        self.rows = []

    def add(self, status, msg):
        self.rows.append((status, msg))
        print(f"  {status:4}  {msg}")

    def check(self, name, actual, spec):
        self.add(*bench_lib.check_value(name, actual, spec))

    @property
    def n_fail(self):
        return sum(1 for s, _ in self.rows if s == "FAIL")

    @property
    def n_pass(self):
        return sum(1 for s, _ in self.rows if s == "PASS")


def dotted(d, path):
    for part in path.split("."):
        if isinstance(d, list) and part.isdigit() and int(part) < len(d):
            d = d[int(part)]
        elif isinstance(d, dict) and part in d:
            d = d[part]
        else:
            return None
    return d


def verify_json_benchmark(name, out_dir, expected_path):
    exp = bench_lib.load_json(expected_path)
    rep = Report(name)
    res_path = Path(out_dir) / "result.json"
    if not res_path.is_file():
        rep.add("FAIL", f"{res_path} not found")
        return rep
    res = bench_lib.load_json(res_path)
    # Fixed-configuration values depend on the table interpolation of the potential; outputs
    # written before version 1.2 carry no "interpolation" field and used the B-spline tables.
    interp = res.get("interpolation", "bspline")
    print(f" table interpolation: {interp}")
    for key, spec in exp["checks"].items():
        by = spec.get("expected_by_interpolation")
        if by is not None:
            if interp not in by:
                rep.add("FAIL", f"{key}: no expected value for interpolation '{interp}'")
                continue
            spec = dict(spec, expected=by[interp], tol=spec["rel_tol"] * abs(by[interp]))
        val = dotted(res, key)
        rep.check(key, None if val is None else float(val), spec)
    return rep


def verify_md_benchmark(name, out_dir, expected_path):
    exp = bench_lib.load_json(expected_path)
    rep = Report(name)
    summaries = {}
    for stage, st in exp["stages"].items():
        h5 = Path(out_dir) / stage / "series.h5"
        print(f" [{stage.upper()}] {h5}")
        if not h5.is_file():
            rep.add("FAIL", f"{h5} not found")
            continue
        s = bench_lib.summarize_run(h5)
        summaries[stage] = s
        if (s["n_steps"], s["n_thermalization"]) != (st["n_steps"], st["n_thermalization"]):
            rep.add("FAIL", f"{stage}: run length {s['n_thermalization']}+"
                    f"{s['n_steps'] - s['n_thermalization']} steps differs from the reference "
                    f"{st['n_thermalization']}+{st['n_steps'] - st['n_thermalization']}; "
                    "the tolerances apply to the full-length run only")
        s = dict(s, econs_drift_abs_per_ns=abs(s["econs_drift_per_ns"]))
        for q, spec in st["quantities"].items():
            rep.check(f"{stage}.{q}", s.get(q), spec)
    zk = exp.get("zero_kelvin_reference")
    if zk and "nvt" in summaries:
        print(" [0 K reference]")
        for q, ref in zk["values_GPa"].items():
            dev = (summaries["nvt"][q] - ref) / ref
            rep.check(f"nvt.{q} vs 0 K reference {ref} GPa (relative deviation)", dev,
                      {"min": -zk["max_rel_deviation"], "max": zk["max_rel_deviation"]})
    return rep


def run_one(name, out_dir):
    print(f"=== benchmark {name}: {out_dir}")
    if name == "smoke":
        return verify_json_benchmark(name, out_dir, EXPECTED / f"{name}.json")
    return verify_md_benchmark(name, out_dir, EXPECTED / f"{name}.json")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("benchmark", choices=BENCHMARKS + ("all",))
    p.add_argument("output_dir", type=Path)
    a = p.parse_args()

    if a.benchmark == "all":
        reports = [run_one(b, a.output_dir / b) for b in BENCHMARKS
                   if (a.output_dir / b).is_dir()]
        if not reports:
            raise SystemExit(f"no benchmark output directories under {a.output_dir}")
    else:
        reports = [run_one(a.benchmark, a.output_dir)]

    print("=== summary")
    n_fail = 0
    for r in reports:
        state = "PASS" if r.n_fail == 0 else "FAIL"
        print(f"  {state:4}  {r.title}: {r.n_pass} passed, {r.n_fail} failed")
        n_fail += r.n_fail
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
