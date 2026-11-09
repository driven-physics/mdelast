#!/usr/bin/env python3
"""
Shape, continuity and accuracy tests of the piecewise-quintic potential tables on real
potential files, with plots. Compares with the cubic B-spline tables of versions <= 1.1, which
the same build provides as interpolation="bspline", where relevant.

For every table of a SETFL file (phi_ij(r), rho_i(r), F_i(rho)) the node data are built
exactly as get_potential.py builds them (raw file values, f' and f'' at the nodes from the
scipy quintic spline k=5, s=0), and the per-interval coefficients are reconstructed in
extended precision with the formulas of src/poly_table.hpp.

Tests
  A  reconstruction     the reconstructed polynomials equal the C++ evaluation (so the exact
                        root analysis below applies to the compiled tables)
  B  roots / extrema / inflections
                        real roots of f, f', f'' of every interval polynomial in [0, 1]
                        (eigenvalues of the companion matrix, cross-checked by sampling);
                        each is compared with the features implied by the raw data alone:
                        sign changes of f_i, of the first differences (extrema) and of the
                        second differences (inflections). An interpolant feature with no data
                        feature within a few grid spacings is SPURIOUS.
                        Also: minimum of rho(r) on [0, r_c] (must stay >= 0).
  C  continuity         jumps of derivatives 0..5 across every interior node (coefficients,
                        long double) and of f, f', f'' through the C++ API at x_i -+ 1 ulp
  D  ends               node values = input arrays; first/last intervals vs the scipy spline;
                        values left at the cut-off (truncation jump of the potential, the same
                        for any interpolant); end-node derivatives vs one-sided differences
  E  accuracy vs data   leave-every-other-node-out: tables built from the even nodes (2h),
                        evaluated at the odd nodes, against the raw values (f) and against the
                        full-grid spline (f', f''), new vs old implementation
  F  data noise         rounding noise of the file values from 7th differences, and the f''
                        noise it implies for ANY interpolant (5.77 sigma / h^2)
  H  data kinks         slope/curvature breaks in the file values (third-difference spikes);
                        no C2 interpolant can follow them, so features next to a kink are
                        reported separately as ringing at a data kink
  G  r = 0 handling     the loader sets phi(0) = phi(h); effect of that node on f'' in the
                        physical range (refit without node 0)

Usage: test_table_shape.py --file ZnCu.eam.alloy --out-dir DIR [--repo .] [--build-dir build/Release]
"""
import argparse
import json
import math
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import UnivariateSpline  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--file", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--repo", default=os.path.dirname(HERE))
ap.add_argument("--build-dir", default=None, help="compiled module (default <repo>/build/Release)")
ap.add_argument("--r-core", type=float, default=1.5, help="physical range of phi, rho starts here [A]")
args = ap.parse_args()
args.build_dir = os.path.abspath(args.build_dir or os.path.join(args.repo, "build", "Release"))
os.makedirs(args.out_dir, exist_ok=True)
TAG = os.path.splitext(os.path.basename(args.file))[0].replace(".eam", "")

sys.path.insert(0, args.build_dir)
sys.path.insert(0, os.path.join(args.repo, "src", "python"))
from md_driver_module import Units  # noqa: E402
import get_potential  # noqa: E402

A, EV = Units.angstrom, Units.ev
C_NEW, C_OLD, C_DATA, C_BAD, C_INK = "#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C_INK, "axes.linewidth": 0.6,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.5,
                     "lines.linewidth": 1.2, "legend.frameon": False})


# ---------------------------------------------------------------- read the SETFL file -----
def read_setfl(fname):
    path = os.path.join(args.repo, "potentials", fname)
    lines = open(path).read().splitlines()
    ne = int(lines[3].split()[0])
    sym = lines[3].split()[1:1 + ne]
    p = lines[4].split()
    Nrho, drho, Nr, dr = int(p[0]), float(p[1]) / A**3, int(p[2]), float(p[3]) * A
    r = np.linspace(0.0, dr * (Nr - 1.0), Nr)
    rho = np.linspace(0.0, drho * (Nrho - 1.0), Nrho)
    tabs, index = {}, 5
    for e in range(ne):
        index += 1
        d = np.loadtxt(path, skiprows=index, max_rows=Nrho + Nr)
        tabs[f"F_{sym[e]}"] = dict(x=rho, h=drho, f=d[:Nrho] * EV, kind="F")
        tabs[f"rho_{sym[e]}"] = dict(x=r, h=dr, f=d[Nrho:Nrho + Nr] / A**3, kind="rho")
        index += Nrho + Nr
    i = j = 0
    for _ in range(ne * (ne + 1) // 2):
        phir = np.loadtxt(path, skiprows=index, max_rows=Nr)
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = phir / r * EV * A
        phi[0] = phi[1]                     # as get_potential.py
        tabs[f"phi_{sym[i]}{sym[j]}"] = dict(x=r, h=dr, f=phi, kind="phi")
        index += Nr
        if i == j:
            j, i = j + 1, 0
        else:
            i += 1
    for t in tabs.values():
        t["f1"], t["f2"] = get_potential.get_derivatives(x=t["x"], y=t["f"])
        # display units: x in A (r) or A^-3 (rho); f in eV (phi, F) or A^-3 (rho)
        t["xs"] = 1 / A if t["kind"] != "F" else A**3
        t["fs"] = A**3 if t["kind"] == "rho" else 1 / EV
        t["xu"] = "r [Å]" if t["kind"] != "F" else "ρ [Å⁻³]"
        t["fu"] = {"phi": "eV", "rho": "Å⁻³", "F": "eV"}[t["kind"]]
    return sym, tabs, r[-1]


sym, T, RC = read_setfl(args.file)
print(f"=== {args.file}: {len(T)} tables, r_c = {RC / A:.6f} A")


# ---------------------------------------------------------------- coefficients ------------
def coefficients(t):
    ld = np.longdouble
    f, f1, f2, h = (np.asarray(t[k], ld) for k in ("f", "f1", "f2", "h"))
    dy = f[1:] - f[:-1]
    D0, D1 = h * f1[:-1], h * f1[1:]
    S0, S1 = h * h * f2[:-1], h * h * f2[1:]
    c = np.stack([f[:-1], D0, 0.5 * S0,
                  10 * dy - 6 * D0 - 4 * D1 - 0.5 * (3 * S0 - S1),
                  -15 * dy + 8 * D0 + 7 * D1 + 0.5 * (3 * S0 - 2 * S1),
                  6 * dy - 3 * (D0 + D1) - 0.5 * (S0 - S1)], axis=1)
    return c.astype(np.float64), c


def deriv_coef(c, k):
    """coefficients (ascending in t) of d^k/dt^k of every interval polynomial"""
    m = np.arange(k, 6)
    fac = np.array([math.factorial(int(q)) // math.factorial(int(q) - k) for q in m], float)
    return c[:, k:] * fac


def eval_poly(c, k, i, tt, h):
    d = deriv_coef(c, k)
    s = np.zeros_like(tt, dtype=c.dtype)
    for q in range(d.shape[1] - 1, -1, -1):
        s = s * tt + d[i, q]
    return s / np.asarray(h, c.dtype) ** k


def roots_unit(d):
    """real roots in [0, 1] of each row's polynomial (ascending coefficients)."""
    n, m = d.shape
    sc = np.max(np.abs(d), axis=1)
    zero = sc == 0
    out = []
    deg = m - 1
    lead = d[:, -1]
    good = (~zero) & (np.abs(lead) > 1e-12 * sc)
    if np.any(good):
        a = d[good] / lead[good, None]
        comp = np.zeros((a.shape[0], deg, deg))
        comp[:, 1:, :-1] = np.eye(deg - 1)
        comp[:, :, -1] = -a[:, :-1]
        ev = np.linalg.eigvals(comp)
        idx = np.nonzero(good)[0]
        for row, e in zip(idx, ev):
            r = e.real[(np.abs(e.imag) < 1e-7) & (e.real >= -1e-12) & (e.real < 1 - 1e-12)]
            out += [(row, float(v)) for v in np.sort(r)]
    for row in np.nonzero((~zero) & (~good))[0]:    # degenerate leading coefficient
        r = np.roots(d[row][::-1])
        r = r.real[(np.abs(r.imag) < 1e-7) & (r.real >= -1e-12) & (r.real < 1 - 1e-12)]
        out += [(row, float(v)) for v in np.sort(r)]
    return out, int(np.sum(zero))


def sign_changes(v):
    """indices j such that v[j] and the next nonzero value after j have opposite signs."""
    s = np.sign(v)
    nz = np.nonzero(s)[0]
    return nz[:-1][s[nz[:-1]] != s[nz[1:]]]


def call_eval(interp, mode, data):
    inp = os.path.join(args.out_dir, f"_in_{os.getpid()}.npz")
    out = os.path.join(args.out_dir, f"_out_{os.getpid()}.npz")
    np.savez(inp, **data)
    cmd = [sys.executable, os.path.join(HERE, "eval_tables.py"), "--repo", args.repo,
           "--build-dir", args.build_dir, "--interpolation", interp, "--mode", mode, "--inp", inp, "--out", out] + (["--file", args.file] if mode == "native" else [])
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        sys.exit(res.stderr[-3000:])
    r = dict(np.load(out))
    os.remove(inp), os.remove(out)
    return {k[2:]: v for k, v in r.items()}


def phys_mask(t, x):
    if t["kind"] == "F":
        return np.ones_like(x, bool)
    return x >= args.r_core * A


# ---------------------------------------------------------------- evaluation points --------
rng = np.random.default_rng(3)
pts = {}
for n, t in T.items():
    x, h = t["x"], t["h"]
    N = len(x) - 1
    rand = np.sort(rng.uniform(x[0], x[-1], 20000))
    nodes = x[1:-1]
    ulp = np.concatenate([np.nextafter(nodes, -np.inf), np.nextafter(nodes, np.inf)])
    dense = (np.arange(N)[:, None] + (np.arange(8) / 8)[None, :]).ravel() * h  # 8 / interval
    ends = np.concatenate([x[0] + h * np.linspace(0, 3, 61), x[-1] - h * np.linspace(0, 3, 61)[::-1],
                           [np.nextafter(x[-1], -np.inf)]])
    pts[n] = dict(rand=rand, ulp=ulp, dense=dense, ends=ends, nodes=x.copy())
seg = {n: {k: len(v) for k, v in p.items()} for n, p in pts.items()}
allx = {f"x_{n}": np.concatenate([p[k] for k in ("rand", "ulp", "dense", "ends", "nodes")]) for n, p in pts.items()}


def split(v, n):
    o, out = 0, {}
    for k in ("rand", "ulp", "dense", "ends", "nodes"):
        out[k] = v[o:o + seg[n][k]]
        o += seg[n][k]
    return out


print("evaluating new build (native objects) ...", flush=True)
NEWV = {n: split(v, n) for n, v in call_eval("quintic", "native", allx).items()}
print("evaluating old build (native objects) ...", flush=True)
OLDV = {n: split(v, n) for n, v in call_eval("bspline", "native", allx).items()}

# leave-every-other-node-out tables (E)
half = {}
for n, t in T.items():
    xe, fe = t["x"][::2], t["f"][::2]
    f1e, f2e = get_potential.get_derivatives(x=xe, y=fe)
    half.update({f"h_{n}": 2 * t["h"], f"f_{n}": fe, f"f1_{n}": f1e, f"f2_{n}": f2e,
                 f"x_{n}": t["x"][1:-1:2]})
print("evaluating leave-out tables (new, old) ...", flush=True)
HALF_NEW = call_eval("quintic", "custom", half)
HALF_OLD = call_eval("bspline", "custom", half)

# ---------------------------------------------------------------- per-table analysis -------
R = {}
for n, t in T.items():
    x, h, f = t["x"], t["h"], t["f"]
    xs, fs = t["xs"], t["fs"]
    c, cl = coefficients(t)
    N = len(x) - 1
    res = {}
    nv, ov = NEWV[n], OLDV[n]
    spl = UnivariateSpline(x, f, k=5, s=0)

    # A: reconstruction vs C++
    i = np.clip((pts[n]["rand"] / h).astype(int), 0, N - 1)
    tt = pts[n]["rand"] / h - i
    rec = np.stack([eval_poly(c, k, i, tt, h) for k in range(3)], 1)
    scale = np.max(np.abs(nv["nodes"]), axis=0)
    res["A_recon_vs_cpp"] = (np.max(np.abs(rec - nv["rand"]), axis=0) / scale).tolist()

    # B: roots of f, f', f'' and data-implied features
    dfd = [np.concatenate([sign_changes(f) + 0.5, np.nonzero(f == 0)[0].astype(float)]),
           np.concatenate([sign_changes(np.diff(f)) + 1.0, np.nonzero(np.diff(f) == 0)[0] + 0.5]),
           sign_changes(np.diff(f, 2)) + 1.0]
    # H: kinks in the data (slope or curvature breaks): spikes of |third difference| over its
    # running median (+-50 nodes). No smooth interpolant can follow a kink; features within
    # KW grid spacings of one are reported as ringing at a data kink, not as spurious.
    d3 = np.abs(np.diff(f, 3))
    from scipy.ndimage import median_filter
    rmed = median_filter(d3, size=101, mode="nearest")
    z = d3 / np.maximum(rmed, 1e-300)
    cand = np.nonzero((z > 300) & (d3 > 1e-6 * np.max(np.abs(f)) * (h / h)))[0]
    kinks = []
    if len(cand):
        for g in np.split(cand, np.nonzero(np.diff(cand) > 3)[0] + 1):
            j = int(g[np.argmax(z[g])])
            # slope jump: quadratic fits on each side, evaluated at the kink node j + 1.5
            xl, xr_ = np.arange(j - 14, j - 1), np.arange(j + 4, j + 17)
            if xl[0] < 0 or xr_[-1] >= len(f):
                continue
            pl, pr = np.polyfit(xl - (j + 1.5), f[xl], 2), np.polyfit(xr_ - (j + 1.5), f[xr_], 2)
            kinks.append(dict(x=float((j + 1.5) * h * xs), z=float(z[j]),
                              dslope=float((pr[1] - pl[1]) / h * fs / xs),
                              df2=float(2 * (pr[0] - pl[0]) / h**2 * fs / xs**2), node=j + 1.5))
    res["H_data_kinks"] = kinks
    KW = 8.0
    knodes = np.array([k_["node"] for k_ in kinks])
    # local rounding noise of the file values (running median of |7th difference|) and the
    # noise floor it implies for f, f', f'' of any interpolant (max of the Hermite value-basis
    # derivatives: 1, 1.875/h, 5.774/h^2). A feature where |f^(k)| stays below NOISE_FAC times
    # this floor within +-2 intervals is classified as noise-level (it reflects the last digits
    # of the file, not the shape of the function).
    NOISE_FAC = 30.0
    d7 = np.abs(np.diff(f, 7)) / 0.6745 / math.sqrt(3432.0)
    sig = np.concatenate([np.full(4, np.nan), median_filter(d7, size=101, mode="nearest"), np.full(3, np.nan)])
    sig = np.where(np.isnan(sig), np.interp(np.arange(len(f)), np.nonzero(~np.isnan(sig))[0],
                                            sig[~np.isnan(sig)]), sig)
    sig = np.maximum(sig, np.finfo(float).eps * np.max(np.abs(f)))   # absolute floor: last digit of max|f|
    nodal = [np.abs(f), np.abs(t["f1"]), np.abs(t["f2"])]
    cfac = [1.0, 1.875 / h, 5.774 / h**2]

    def noise_level(k, v):
        a, b = max(0, int(v) - 2), min(len(f), int(v) + 4)
        return np.max(nodal[k][a:b]) < NOISE_FAC * cfac[k] * np.max(sig[a:b])

    def at_end(v):
        return v < 3 or v > N - 3
    tol = [1.0, 1.5, 2.5]                        # matching window in grid spacings
    res["B"] = {}
    feats = {}
    for k in range(3):
        rts, nzero = roots_unit(deriv_coef(c, k))
        xr = np.array([(row + tr) for row, tr in rts])           # in units of h
        last = eval_poly(c, k, np.array([N - 1]), np.array([1.0]), h)[0]
        if last == 0.0:                       # root exactly at the last node (t = 1 of last interval)
            xr = np.append(xr, float(N))
        # exclude x = 0 node for phi (flat by construction) only through the region split below
        dat = dfd[k]
        sup = np.array([np.any(np.abs(dat - v) <= tol[k]) for v in xr], bool) if len(xr) else np.zeros(0, bool)
        matched = np.array([np.any(np.abs(xr - v) <= tol[k]) for v in dat], bool) if len(dat) else np.zeros(0, bool)
        # sampling cross-check (64 points per interval, exact polynomials)
        tg = np.linspace(0, 1, 65)[:-1]
        dk = deriv_coef(c, k)
        vals = np.zeros((N, 64))
        for q in range(dk.shape[1] - 1, -1, -1):
            vals = vals * tg[None, :] + dk[:, q:q + 1]
        vals = vals.ravel()
        n_samp = len(sign_changes(vals))
        reg = {}
        for rname, lo, hi in ((("core", 0, args.r_core * A), ("phys", args.r_core * A, np.inf))
                              if t["kind"] != "F" else (("full", -np.inf, np.inf),)):
            mi = (xr * h >= lo) & (xr * h < hi)
            md = (dat * h >= lo) & (dat * h < hi)
            dnoise = np.array([noise_level(k, v) for v in dat], bool)
            atk = np.array([len(knodes) > 0 and np.any(np.abs(knodes - v) <= KW) for v in xr], bool)
            nsl = np.array([noise_level(k, v) for v in xr], bool) & ~sup & ~atk
            end = np.array([at_end(v) for v in xr], bool) & ~sup & ~atk & ~nsl
            atk = atk | nsl | end      # everything below is "explained"; counted separately
            reg_extra = dict(noise_level=int((mi & nsl).sum()), table_end=int((mi & end).sum()))
            reg[rname] = dict(interp=int(mi.sum()), data=int(md.sum()),
                              spurious=int((mi & ~sup & ~atk).sum()),
                              at_data_kink=int((mi & ~sup & atk).sum()) - reg_extra["noise_level"] - reg_extra["table_end"],
                              **reg_extra,
                              missing=int((md & ~matched & ~dnoise).sum()),
                              missing_noise_level=int((md & ~matched & dnoise).sum()),
                              spurious_x=(xr[mi & ~sup & ~atk] * h * xs).round(5).tolist()[:20],
                              kink_x=(xr[mi & ~sup & atk] * h * xs).round(5).tolist()[:20])
        res["B"][["roots_f", "extrema", "inflections"][k]] = dict(
            regions=reg, total_interp=len(xr), total_sampled=n_samp, zero_intervals=nzero)
        atk_all = np.array([(len(knodes) > 0 and np.any(np.abs(knodes - v) <= KW)) or noise_level(k, v) or at_end(v)
                            for v in xr], bool)
        feats[k] = (xr, sup | atk_all, dat, atk_all & ~sup)
    if t["kind"] == "rho":
        rt1 = feats[1][0]
        cand = np.concatenate([[0, N], rt1])
        ii = np.clip(np.floor(cand).astype(int), 0, N - 1)
        vv = eval_poly(c, 0, ii, cand - ii, h)
        res["B"]["rho_min"] = dict(value=float(np.min(vv) * fs), at=float(cand[np.argmin(vv)] * h * xs),
                                   data_min=float(np.min(f) * fs), data_min_at=float(x[np.argmin(f)] * xs))

    # C: continuity, k = 0..5, relative to max |f^(k)| at the nodes of the physical range
    pm = phys_mask(t, x[1:-1])
    jumps = []
    for k in range(6):
        L = eval_poly(cl, k, np.arange(0, N - 1), np.ones(N - 1, np.longdouble), np.longdouble(h))
        Rr = eval_poly(cl, k, np.arange(1, N), np.zeros(N - 1, np.longdouble), np.longdouble(h))
        J = np.abs((L - Rr).astype(float))
        sc = np.max(np.abs(Rr.astype(float))[pm])
        jumps.append(J / sc)
    res["C_jump_phys"] = [float(np.max(j[pm])) for j in jumps]
    res["C_jump_core"] = [float(np.max(j[~pm])) if np.any(~pm) else None for j in jumps]
    u = nv["ulp"]
    m = len(x) - 2
    res["C_api_ulp"] = (np.max(np.abs(u[:m] - u[m:])[np.r_[pm]], axis=0) / scale).tolist()
    uo = ov["ulp"]
    res["C_api_ulp_old"] = (np.max(np.abs(uo[:m] - uo[m:])[np.r_[pm]], axis=0) / scale).tolist()

    # D: ends
    inp = np.stack([f, t["f1"], t["f2"]], 1)
    res["D_nodes_vs_input"] = (np.max(np.abs(nv["nodes"] - inp), axis=0) / scale).tolist()
    xe = pts[n]["ends"][:-1]
    ref = np.stack([spl(xe), spl.derivative(1)(xe), spl.derivative(2)(xe)], 1)
    res["D_end_intervals_vs_spline"] = (np.max(np.abs(nv["ends"][:-1] - ref), axis=0) / scale).tolist()
    res["D_end_intervals_vs_spline_old"] = (np.max(np.abs(ov["ends"][:-1] - ref), axis=0) / scale).tolist()
    lastv, lastv_old = nv["ends"][-1], ov["ends"][-1]
    if t["kind"] != "F":
        u3 = np.array([fs, fs / xs, fs / xs**2])
        res["D_value_at_cutoff_new"] = (lastv * u3).tolist()
        res["D_value_at_cutoff_old"] = (lastv_old * u3).tolist()
        res["D_value_at_cutoff_rel"] = (np.abs(lastv) / scale).tolist()
    # end-node derivative from the spline vs one-sided 6th-order difference of raw data
    w = np.array([-49 / 20, 6, -15 / 2, 20 / 3, -15 / 4, 6 / 5, -1 / 6])
    d_lo = (w @ f[:7]) / h
    d_hi = -(w @ f[::-1][:7]) / h
    res["D_end_deriv_vs_FD"] = [float(abs(d_lo - t["f1"][0]) / scale[1]), float(abs(d_hi - t["f1"][-1]) / scale[1])]

    # E: leave-out accuracy
    xo = t["x"][1:-1:2]
    pmo = phys_mask(t, xo)
    refo = np.stack([f[1:-1:2], t["f1"][1:-1:2], t["f2"][1:-1:2]], 1)
    en = np.abs(HALF_NEW[n] - refo)
    eo = np.abs(HALF_OLD[n] - refo)
    res["E_leaveout_new"] = (np.max(en[pmo], axis=0) / scale).tolist()
    res["E_leaveout_old"] = (np.max(eo[pmo], axis=0) / scale).tolist()
    io = np.arange(1, len(x) - 1, 2)
    clean = pmo & (io > 10) & (io < len(x) - 11)
    for kn in knodes:
        clean &= np.abs(io - kn) > 2 * KW
    res["E_leaveout_new_clean"] = (np.max(en[clean], axis=0) / scale).tolist()
    res["E_leaveout_old_clean"] = (np.max(eo[clean], axis=0) / scale).tolist()

    # F: data noise
    d7 = np.diff(f, 7)[phys_mask(t, x[:-7])]
    sigma = float(np.median(np.abs(d7)) / 0.6745 / math.sqrt(3432.0))
    res["F_noise_sigma_rel"] = sigma / scale[0]
    res["F_implied_f2_noise_rel"] = 5.774 * sigma / h**2 / scale[2]

    # G: r = 0 handling (phi only)
    if t["kind"] == "phi":
        f1b, f2b = get_potential.get_derivatives(x=x[1:], y=f[1:])
        pmb = phys_mask(t, x[1:])
        res["G_r0_effect_on_f2_phys"] = float(np.max(np.abs(f2b - t["f2"][1:])[pmb]) / scale[2])
        res["G_r0_effect_on_f1_phys"] = float(np.max(np.abs(f1b - t["f1"][1:])[pmb]) / scale[1])
    R[n] = res

    # --------------------------------------------------------- plots ---------------------
    xd = pts[n]["dense"]
    vd, vdo = nv["dense"], ov["dense"]
    labels = ["f", "f′", "f″"]
    units = [t["fu"], f"{t['fu']}/{t['xu'].split()[1][1:-1]}", f"{t['fu']}/{t['xu'].split()[1][1:-1]}²"]
    fig, axs = plt.subplots(3, 2, figsize=(11, 8.5), sharex="col")
    if t["kind"] == "F":
        zooms = [(x[0], x[-1]), (0.15 * x[-1], 0.85 * x[-1])]
    else:
        zooms = [(x[0], x[-1]), (args.r_core * A, x[-1])]
    for col, (lo, hi) in enumerate(zooms):
        msk = (xd >= lo) & (xd <= hi)
        for k in range(3):
            ax = axs[k, col]
            ax.plot(xd[msk] * xs, vd[msk, k] * [fs, fs / xs, fs / xs**2][k],
                    color=C_NEW, lw=1.0, label="new (quintic)")
            if k == 2:
                ax.plot(xd[msk] * xs, vdo[msk, 2] * fs / xs**2, color=C_OLD, lw=0.7, alpha=0.8,
                        label="bspline (cubic B-spline of f″)")
            if k == 0:
                nm = (x >= lo) & (x <= hi)
                step = max(1, nm.sum() // 120)
                ax.plot(x[nm][::step] * xs, f[nm][::step] * fs, "o", ms=2.5, color=C_DATA,
                        label=f"file values (every {step}th)")
            xr, sup, dat, atk = feats[k]
            sel = (xr * h >= lo) & (xr * h <= hi)
            ii = np.clip(np.floor(xr[sel]).astype(int), 0, N - 1)
            yv = eval_poly(c, k, ii, xr[sel] - ii, h) * [fs, fs / xs, fs / xs**2][k]
            ok = sup[sel]
            ak = atk[sel]
            ax.plot(xr[sel][ok & ~ak] * h * xs, yv[ok & ~ak], "|", ms=9, mew=1.2, color=C_INK,
                    label=["zero crossings", "extrema", "inflections"][k] + " (supported by data)")
            if np.any(ak):
                ax.plot(xr[sel][ak] * h * xs, yv[ak], "o", ms=6, mfc="none", mew=1.2, color="#eda100",
                        label="explained: data kink / noise-level / table end")
            for kk in kinks:
                if lo <= kk["node"] * h <= hi:
                    ax.axvline(kk["x"], color="#eda100", lw=0.8, ls="--", zorder=0)
            if np.any(~ok):
                ax.plot(xr[sel][~ok] * h * xs, yv[~ok], "x", ms=8, mew=1.6, color=C_BAD,
                        label="spurious (no data feature nearby)")
            if col == 0 and t["kind"] != "F":
                ax.axvspan(x[0] * xs, args.r_core, color="#f0efec", zorder=0)
            ax.set_ylabel(f"{labels[k]} [{units[k]}]")
            if col == 1 and t["kind"] != "F":
                yr = vd[msk & (xd > args.r_core * A * 1.15), k] * [fs, fs / xs, fs / xs**2][k]
                if len(yr):
                    pad = 0.08 * (yr.max() - yr.min() + 1e-30)
                    ax.set_ylim(yr.min() - pad, yr.max() + pad)
            ax.legend(fontsize=7, loc="best")
        axs[2, col].set_xlabel(t["xu"])
        axs[0, col].set_title(["full table" + (" (grey: r < %.1f Å, core)" % args.r_core if t["kind"] != "F" else ""),
                               "physical range"][col])
    fig.suptitle(f"{TAG}: {n} — new table, its zero crossings / extrema / inflections, and the old f″", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, f"{TAG}_{n}_shape.png"), dpi=130)
    plt.close(fig)

    # continuity + leave-out error figure
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
    ax = axs[0]
    for k, col in zip(range(6), [C_NEW, "#1baf7a", "#eda100", "#e87ba4", "#6250d6", C_INK]):
        ax.semilogy(x[1:-1] * xs, np.maximum(jumps[k], 1e-20), ".", ms=1.5, color=col,
                    label=f"d{k}" + (" (imposed)" if k <= 2 else ""))
    ax.set_xlabel(t["xu"])
    ax.set_ylabel("|jump at node| / max|f^(k)|")
    ax.set_title("continuity across nodes (derivatives 0–5)")
    ax.legend(fontsize=7, ncol=3, markerscale=5)
    ax = axs[1]
    for k, ls in zip(range(3), ["-", "--", ":"]):
        ax.semilogy(xo[pmo] * xs, np.maximum(en[pmo, k] / scale[k], 1e-18), ls, color=C_NEW, lw=0.8,
                    label=f"quintic {labels[k]}")
        ax.semilogy(xo[pmo] * xs, np.maximum(eo[pmo, k] / scale[k], 1e-18), ls, color=C_OLD, lw=0.8,
                    label=f"bspline {labels[k]}")
    ax.set_xlabel(t["xu"])
    ax.set_ylabel("|error| / max|f^(k)|")
    ax.set_title("accuracy, table from even nodes (2h) at odd nodes")
    ax.legend(fontsize=7, ncol=3)
    fig.suptitle(f"{TAG}: {n}", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, f"{TAG}_{n}_continuity_accuracy.png"), dpi=130)
    plt.close(fig)

# ---------------------------------------------------------------- report -------------------
json.dump(R, open(os.path.join(args.out_dir, f"{TAG}_results.json"), "w"), indent=1)
fmt = lambda v: " ".join(f"{a:8.1e}" if a is not None else "     -  " for a in v)  # noqa: E731
for n, r in R.items():
    print(f"--- {n}")
    print(f"  A recon vs C++        f/f'/f''   {fmt(r['A_recon_vs_cpp'])}")
    for key, b in r["B"].items():
        if key == "rho_min":
            print(f"  B min rho on [0,rc]   {b['value']:.3e} A^-3 at r = {b['at']:.4f} A  (file data min {b['data_min']:.3e} at {b['data_min_at']:.4f} A)")
            continue
        regs = "  ".join(f"{rn}: interp {v['interp']} data {v['data']} SPURIOUS {v['spurious']} | at-kink {v['at_data_kink']} noise {v['noise_level']} end {v['table_end']} | missing {v['missing']}"
                         for rn, v in b["regions"].items())
        print(f"  B {key:12s} {regs}  (eig {b['total_interp']}, sampled {b['total_sampled']})")
        for rn, v in b["regions"].items():
            if v["spurious"]:
                print(f"      spurious {key} in {rn} at x = {v['spurious_x']}")
    print(f"  C jump phys  d0..d5   {fmt(r['C_jump_phys'])}")
    if r["C_jump_core"][0] is not None:
        print(f"  C jump core  d0..d5   {fmt(r['C_jump_core'])}")
    print(f"  C API +-1ulp new/old  {fmt(r['C_api_ulp'])} | {fmt(r['C_api_ulp_old'])}")
    print(f"  D nodes vs input      {fmt(r['D_nodes_vs_input'])}")
    print(f"  D end intervals vs spline new/old {fmt(r['D_end_intervals_vs_spline'])} | {fmt(r['D_end_intervals_vs_spline_old'])}")
    if "D_value_at_cutoff_new" in r:
        print(f"  D at cutoff (f, f', f'') new {fmt(r['D_value_at_cutoff_new'])} old {fmt(r['D_value_at_cutoff_old'])}  rel {fmt(r['D_value_at_cutoff_rel'])}")
    print(f"  D end f' vs FD (lo, hi) {fmt(r['D_end_deriv_vs_FD'])}")
    print(f"  E leave-out new / old {fmt(r['E_leaveout_new'])} | {fmt(r['E_leaveout_old'])}")
    print(f"  E   (away from ends and data kinks) {fmt(r['E_leaveout_new_clean'])} | {fmt(r['E_leaveout_old_clean'])}")
    for kk in r["H_data_kinks"]:
        print(f"  H DATA KINK at x = {kk['x']:.4f}: |d3f|/running median {kk['z']:.1e}, slope jump {kk['dslope']:+.3e}, f'' jump {kk['df2']:+.3e}")
    print(f"  F data noise sigma/max|f| {r['F_noise_sigma_rel']:.1e}  -> f'' noise floor {r['F_implied_f2_noise_rel']:.1e}")
    if "G_r0_effect_on_f2_phys" in r:
        print(f"  G phi(0)=phi(h) effect in phys range: f' {r['G_r0_effect_on_f1_phys']:.1e}  f'' {r['G_r0_effect_on_f2_phys']:.1e}")
