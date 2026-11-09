#!/usr/bin/env python3
"""Render benchmarks/expected/*.json as the human-readable benchmarks/expected/EXPECTED.md.

Usage (from the repository root): python3 benchmarks/make_expected_md.py
"""

from __future__ import annotations

import json
from pathlib import Path

EXP = Path(__file__).resolve().parent / "expected"

LABEL = {
    "L_L0": "L/L0 (mean box length / initial)", "C11": "C11", "C12": "C12", "C44": "C44",
    "CB11": "C^B_11 (Born)", "CB12": "C^B_12 (Born)", "CB44": "C^B_44 (Born)",
    "fluc11": "stress-fluctuation term, 11", "fluc12": "stress-fluctuation term, 12",
    "fluc44": "stress-fluctuation term, 44", "kin11": "kinetic (ideal-gas) term, 11",
    "T_kin": "T_K", "T_conf_recip": "T_conf, 1/<1/T_conf>", "T_conf_sample": "T_conf, <T_conf>",
    "pressure": "mean pressure", "econs_err_max": "max abs(E_cons(t)/E_cons(0) - 1)",
    "econs_drift_abs_per_ns": "abs(linear drift of E_cons/E_cons(0)) over production",
    "K_npt": "NPT K term",
}


def fmt(x, nd=6):
    if x is None:
        return "-"
    if x == 0:
        return "0"
    a = abs(x)
    if 1e-3 <= a < 1e5:
        return f"{x:.{nd}g}"
    return f"{x:.3g}"


def fmt_values(values):
    """Format a {name: value} dict of elastic constants in GPa as plain text."""
    return ", ".join(f"{k} = {v}" for k, v in values.items()) + " GPa"


def md_table(spec, what):
    out = []
    for stage, st in spec["stages"].items():
        out.append(f"\n**{stage.upper()} stage** ({st['n_thermalization']:,} + "
                   f"{st['n_steps'] - st['n_thermalization']:,} steps; reference: {what}; "
                   f"{st['n_velocity_replicas']} velocity replicas)\n")
        out.append("| quantity | unit | expected | tolerance | per-run SE | velocity-replica std | version 1.0 reference |")
        out.append("|---|---|---|---|---|---|---|")
        for q, s in st["quantities"].items():
            name = LABEL.get(q, q)
            u = s.get("unit", "")
            prev = fmt(s.get("expected_previous"))
            if "info" in s:
                out.append(f"| {name} | {u} | {fmt(s['expected'])} | not checked (INFO) | "
                           f"{fmt(s.get('se_run'), 2)} | - | {prev} |")
            elif "max" in s:
                out.append(f"| {name} | {u} | <= {fmt(s['max'], 2)} | one-sided | "
                           f"reference {fmt(s['reference'], 2)} | max of set {fmt(s['replica_max'], 2)} | - |")
            else:
                out.append(f"| {name} | {u} | {fmt(s['expected'])} | +/- {fmt(s['tol'], 2)} | "
                           f"{fmt(s['se_run'], 2)} | {fmt(s['s_vel'], 2)} | {prev} |")
    return "\n".join(out)


def json_table(spec):
    out = ["| check | unit | expected / bound, quintic tables | expected / bound, bspline tables | "
           "observed, quintic | observed, bspline |", "|---|---|---|---|---|---|"]
    for k, s in spec["checks"].items():
        u = s.get("unit") or ""
        if "expected_by_interpolation" in s:
            by, r = s["expected_by_interpolation"], s["rel_tol"]
            q = f"{by['quintic']:.10g} +/- {r * abs(by['quintic']):.2g}"
            b = f"{by['bspline']:.10g} +/- {r * abs(by['bspline']):.2g}"
            out.append(f"| `{k}` | {u} | {q} | {b} | = expected | = expected |")
            continue
        parts = []
        if "min" in s:
            parts.append(f">= {fmt(s['min'], 3)}")
        if "max" in s:
            parts.append(f"<= {fmt(s['max'], 3)}")
        bound = " and ".join(parts)
        out.append(f"| `{k}` | {u} | {bound} | {bound} | {fmt(s.get('observed_quintic'), 4)} | "
                   f"{fmt(s.get('observed_bspline', s.get('observed')), 4)} |")
    return "\n".join(out)


def main():
    z = json.loads((EXP / "zn06_300K.json").read_text())
    l12 = json.loads((EXP / "l12_cuzn3_5K.json").read_text())
    sm = json.loads((EXP / "smoke.json").read_text())
    zk = l12["zero_kelvin_reference"]
    txt = [
        "# Expected benchmark values",
        "",
        "Generated from `expected/*.json` (the files `verify.py` reads) by `make_expected_md.py`. "
        "All stresses and elastic constants in GPa, temperatures in K. "
        "Cubic averages: C11 = (C11+C22+C33)/3, C12 = (C12+C13+C23)/3, C44 = (C44+C55+C66)/3.",
        "",
        "Tolerance rule for Z1 and L1: tol = 4 x sqrt(2) x max(per-run SE, velocity-replica std); "
        "conserved-energy checks are upper bounds of 3 x the largest value in the reference set. "
        "See README.md, Section 'How the tolerances were set'.",
        "",
        "## Benchmark S: smoke test (`smoke/run_smoke.py`)",
        "",
        "The fixed-cell values depend on the table interpolation of the potential; `verify.py` uses "
        "the column that matches the `interpolation` field of `result.json` (outputs of versions "
        "<= 1.1 have no such field and used the B-spline tables). The other bounds apply to both.",
        "",
        json_table(sm),
        "",
        "## Benchmark Z1: 6 at.% Zn, 300 K, NPT -> NVT (`sample_case/run.py`)",
        md_table(z, z.get("source", {}).get("expected_values", "the run of manuscript Sec. 5.1 (Figs. 1-5)")),
        "",
        "## Benchmark L1: L1_2 CuZn3, 5 K, NPT -> NVT (`l12_cuzn3_5K/run.py`)",
        md_table(l12, l12.get("source", {}).get("expected_values", "replica 0 of the set behind Table 1")),
        "",
        f"Additional check: |C(5 K) - C(0 K)| / C(0 K) <= {zk['max_rel_deviation']} against the 0 K "
        f"reference values {fmt_values(zk['values_GPa'])} ({zk['source']}). Table 1 lists the "
        f"nine-replica means {fmt_values(l12['table_replica_mean_GPa'])}.",
        "",
    ]
    (EXP / "EXPECTED.md").write_text("\n".join(txt))
    print("wrote", EXP / "EXPECTED.md")


if __name__ == "__main__":
    main()
