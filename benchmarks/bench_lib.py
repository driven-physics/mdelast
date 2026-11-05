"""Shared helpers for the mdelast benchmarks (reading ``series.h5``, comparing).

Only NumPy and h5py are needed; the compiled extension is not imported, so the
verification step can run on any machine that has the output files.

Units in ``series.h5`` are CGS: stresses and elastic constants in dyn/cm^2
(1 GPa = 1e10 dyn/cm^2), temperatures in K. Every ``series/<q>/data`` array
holds all MD steps (thermalization + production); production starts at
``series/<q>/cycles_thermalization``. All values returned here are in GPa, K,
or dimensionless.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

GPA = 1e10  # dyn/cm^2 per GPa

# Voigt index groups for the cubic averages used throughout the paper
CUBIC = {
    "11": ((1, 1), (2, 2), (3, 3)),
    "12": ((1, 2), (1, 3), (2, 3)),
    "44": ((4, 4), (5, 5), (6, 6)),
}


def _mean(h5, name):
    return float(h5[f"series/{name}/mean"][()])


def _data(h5, name):
    return np.asarray(h5[f"series/{name}/data"][()], dtype=float)


def _cubic_mean(h5, prefix):
    out = {}
    for key, comps in CUBIC.items():
        out[key] = float(np.mean([_mean(h5, f"{prefix}_{i}_{j}") for i, j in comps]))
    return out


def summarize_run(h5_path) -> dict:
    """Return the benchmark quantities of one mdelast run (NVT or NPT).

    Keys (GPa / K / dimensionless):
      C11, C12, C44           total isothermal constants (cubic averages)
      CB11, CB12, CB44        Born term
      fluc11, fluc12, fluc44  stress-fluctuation term (the subtracted term)
      kin11, kin44            ideal-gas (kinetic) term
      K_npt                   NPT bulk-modulus term (NPT only)
      T_kin                   kinetic temperature, production mean
      T_conf_recip            configurational temperature, 1/<1/T_conf>
      T_conf_sample           configurational temperature, <T_conf> per step
      pressure                mean pressure
      L_L0, V_V0              mean box-length / volume ratio to the initial box
      econs_err_max           max |E_cons(t)/E_cons(0) - 1| over the whole run
      econs_drift_per_ns      slope of E_cons(t)/E_cons(0) - 1 over production [1/ns]
      cubic_spread_C11        max-min of (C_11, C_22, C_33), a symmetry diagnostic
    """
    import h5py

    with h5py.File(Path(h5_path), "r") as h5:
        tr = h5["thermal_results"]
        n0 = int(h5["series/temperature_kin/cycles_thermalization"][()])
        n_steps = int(len(h5["series/time"][()]))
        out = {
            "n_steps": n_steps,
            "n_thermalization": n0,
            "C11": float(tr["elasticity_tensor_11"][()]) / GPA,
            "C12": float(tr["elasticity_tensor_12"][()]) / GPA,
            "C44": float(tr["elasticity_tensor_44"][()]) / GPA,
        }
        full = np.asarray(tr["elasticity_tensor"][()]) / GPA
        out["cubic_spread_C11"] = float(np.ptp(np.diag(full)[:3]))
        for key in CUBIC:
            out[f"CB{key}"] = _mean(h5, f"born_elasticity_tensor_{key}") / GPA
        fluc = _cubic_mean(h5, "stress_fluc")
        kin = _cubic_mean(h5, "elasticity_ideal_gas")
        for key in CUBIC:
            out[f"fluc{key}"] = fluc[key] / GPA
        out["kin11"] = kin["11"] / GPA
        out["kin44"] = kin["44"] / GPA
        kb = float(tr["bulk_modulus"][()]) / GPA if "bulk_modulus" in tr else 0.0
        out["K_npt"] = kb

        out["T_kin"] = _mean(h5, "temperature_kin")
        inv = _data(h5, "temperature_conf_inverse")
        out["T_conf_recip"] = 1.0 / _mean(h5, "temperature_conf_inverse")
        out["T_conf_sample"] = float(np.mean(1.0 / inv[n0:]))
        out["pressure"] = _mean(h5, "pressure") / GPA
        out["L_L0"] = _mean(h5, "L_L0")
        out["V_V0"] = _mean(h5, "V_V0")

        econs = _data(h5, "etot_cons")
        rel = econs / econs[0] - 1.0
        out["econs_err_max"] = float(np.max(np.abs(rel)))
        t_ns = np.asarray(h5["series/time"][()], dtype=float)[n0:] * 1e9
        out["econs_drift_per_ns"] = float(np.polyfit(t_ns, rel[n0:], 1)[0])
    return out


def production_series(h5_path) -> dict:
    """Per-step production series of the scalar quantities (for standard errors)."""
    import h5py

    with h5py.File(Path(h5_path), "r") as h5:
        n0 = int(h5["series/temperature_kin/cycles_thermalization"][()])
        s = {}
        for key in CUBIC:
            s[f"CB{key}"] = _data(h5, f"born_elasticity_tensor_{key}")[n0:] / GPA
        s["T_kin"] = _data(h5, "temperature_kin")[n0:]
        inv = _data(h5, "temperature_conf_inverse")[n0:]
        s["T_conf_sample"] = 1.0 / inv
        s["inv_T_conf"] = inv
        s["pressure"] = _data(h5, "pressure")[n0:] / GPA
        s["L_L0"] = _data(h5, "L_L0")[n0:]
    return s


# ----------------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------------

def check_value(name, actual, spec):
    """Compare one value with its spec; return (status, message).

    spec keys: ``expected`` and ``tol`` (|actual - expected| <= tol), or
    ``max`` / ``min`` (one-sided bounds), or ``info`` (report only).
    """
    if actual is None or (isinstance(actual, float) and not math.isfinite(actual)):
        return "FAIL", f"{name}: missing or non-finite value"
    unit = spec.get("unit", "")
    if spec.get("info"):
        ref = spec.get("expected")
        ref_s = f", reference {ref:.6g}" if ref is not None else ""
        return "INFO", f"{name} = {actual:.6g} {unit}{ref_s} (not checked: {spec['info']})"
    if "expected" in spec:
        exp, tol = float(spec["expected"]), float(spec["tol"])
        d = actual - exp
        status = "PASS" if abs(d) <= tol else "FAIL"
        return status, (f"{name} = {actual:.6g} {unit}, expected {exp:.6g} +/- {tol:.3g}"
                        f" (diff {d:+.3g}, {abs(d) / tol if tol else float('inf'):.2f} of tol)")
    lo, hi = spec.get("min", -math.inf), spec.get("max", math.inf)
    status = "PASS" if lo <= actual <= hi else "FAIL"
    rng = []
    if math.isfinite(lo):
        rng.append(f">= {lo:.3g}")
    if math.isfinite(hi):
        rng.append(f"<= {hi:.3g}")
    return status, f"{name} = {actual:.6g} {unit}, required {' and '.join(rng)}"


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
