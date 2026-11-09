#!/usr/bin/env python3
"""
Helper of test_table_shape.py: evaluate potential tables with one compiled build.

Runs in its own process because the compiled module has the same name in every build.

  native  every table of a SETFL file through the objects MDDriver uses
          (TablePotential.calculate, EAM.rho/drhodr/d2rhodr2, EAM.F/dFdrho/d2Fdrho2)
          at the points x_<name> of the input .npz
  custom  one generic table per name, built from the arrays h_, f_, f1_, f2_ of the input
          .npz with the build's TablePotential (the same interpolant class as EAM in both
          versions), with cut-offs disabled, evaluated at x_<name>

Output .npz: v_<name> with shape (len(x), 3) = (f, f', f''), internal units.
Usage: eval_tables.py --repo R --build-dir B [--interpolation quintic|bspline] --file F
                      --mode native|custom --inp in.npz --out out.npz
"""
import argparse
import contextlib
import os
import sys

import numpy as np

ap = argparse.ArgumentParser()
for a in ("--repo", "--build-dir", "--file", "--mode", "--inp", "--out"):
    ap.add_argument(a, required=a != "--file")
ap.add_argument("--interpolation", choices=("quintic", "bspline"), default="quintic")
args = ap.parse_args()
sys.path.insert(0, os.path.abspath(args.build_dir))
sys.path.insert(0, os.path.join(os.path.abspath(args.repo), "src", "python"))
import md_driver_module as m  # noqa: E402
import get_potential  # noqa: E402

inp = np.load(args.inp)
names = [k[2:] for k in inp.files if k.startswith("x_")]
out = {}


def run(fun, x):
    return np.array([fun(float(v)) for v in x])


# the C++ constructors print one line per object; keep stdout for nothing else
with contextlib.redirect_stdout(open(os.devnull, "w")):
    if args.mode == "native":
        elements, pots, eams = get_potential.get_EAM_alloy_file(args.file, interpolation=args.interpolation)
        sym = [e[0] for e in elements]
        objs = {}
        for j in range(len(sym)):
            for i in range(j + 1):
                p = pots[j * (j + 1) // 2 + i]
                objs[f"phi_{sym[i]}{sym[j]}"] = lambda r, p=p: (lambda v: (v.V, v.dVdr, v.d2Vdr2))(p.calculate(r))
        for e, s in enumerate(sym):
            q = eams[e]
            objs[f"rho_{s}"] = lambda r, q=q: (q.rho(r), q.drhodr(r), q.d2rhodr2(r))
            objs[f"F_{s}"] = lambda x, q=q: (q.F(x), q.dFdrho(x), q.d2Fdrho2(x))
        for n in names:
            out[f"v_{n}"] = run(objs[n], inp[f"x_{n}"])
    else:
        for n in names:
            t = m.TablePotential(1e300, -1e300, float(inp[f"h_{n}"]), inp[f"f_{n}"].tolist(),
                                 inp[f"f1_{n}"].tolist(), inp[f"f2_{n}"].tolist(), args.interpolation)
            out[f"v_{n}"] = run(lambda r, t=t: (lambda v: (v.V, v.dVdr, v.d2Vdr2))(t.calculate(r)),
                                inp[f"x_{n}"])
np.savez(args.out, **out)
