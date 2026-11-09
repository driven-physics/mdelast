#!/usr/bin/env python3
"""Benchmark L1: ordered L1_2 CuZn3 at 5 K (manuscript Table 1).

500 atoms (5 x 5 x 5 cubic L1_2 cells; Cu on cube corners, Zn on face centres),
Clement-Auger Cu-Zn EAM potential, 2 fs timestep, Nose-Hoover chains (M = 10,
tau_T = 50 fs, tau_P = 500 fs), neighbour list with r_skin = 1.2 rebuilt every
20 steps:
  1. NPT at zero pressure: 50,000 thermalization + 100,000 production steps;
  2. NVT at the mean NPT volume, started from the last NPT snapshot:
     50,000 + 100,000 steps, accumulating the elasticity tensor.
This is exactly replica 0 of the nine-replica set behind Table 1
(velocity seeds 0 for NPT and 500 for NVT).

Usage (from the repository root, after building):
    python3 benchmarks/l12_cuzn3_5K/run.py --build-dir build/Release \
        --output-dir benchmark_output/l12_cuzn3_5K
    python3 benchmarks/verify.py l12_cuzn3_5K benchmark_output/l12_cuzn3_5K
Runtime: about 2 x 150,000 steps on one core (see README.md).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

TEMPERATURE = 5.0          # K
A0_ANGSTROM = 3.7816       # initial L1_2 CuZn3 lattice constant [A]
NCELL = 5
DISTORT = 1e-5             # tiny uniform offset of all sites (fraction of a0)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo-root", type=Path, default=HERE.parents[1],
                   help="mdelast repository root (contains src/python and potentials/)")
    p.add_argument("--build-dir", type=Path, default=None,
                   help="directory with md_driver_module*.so (default: <repo>/build/Release)")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--cycles-thermalization", type=int, default=50_000)
    p.add_argument("--production-cycles", type=int, default=100_000)
    p.add_argument("--npt-velocity-seed", type=int, default=0)
    p.add_argument("--nvt-velocity-seed", type=int, default=500)
    p.add_argument("--sample-stride", type=int, default=1,
                   help="compute and store the observables every N MD steps (default 1)")
    p.add_argument("--interpolation", choices=("quintic", "bspline"), default="quintic",
                   help="table interpolation of the potential (default quintic; bspline = version <= 1.1)")
    return p.parse_args()


def main():
    args = parse_args()
    repo = args.repo_root.resolve()
    build = (args.build_dir or repo / "build" / "Release").resolve()
    out = args.output_dir.resolve()
    if out.exists():
        raise SystemExit(f"output directory already exists: {out}")
    sys.path.insert(0, str(build))
    sys.path.insert(0, str(repo / "src" / "python"))
    import h5py  # noqa: F401  (fail early if missing)
    from md_driver_module import MDDriver, Units, Vector3d
    from molecular_dynamics import Material, MolecularDynamics
    import get_potential

    get_potential.DEFAULT_INTERPOLATION = args.interpolation
    elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy")
    name = "-".join(e[0] for e in elements)
    weights = [e[2] for e in elements]
    species_to_type = {e[0]: i for i, e in enumerate(elements)}
    cu, zn = species_to_type["Cu"], species_to_type["Zn"]

    def l12_lattice(a, n):
        L = n * a
        d = a * DISTORT
        basis = [([0.0, 0.0, 0.0], cu), ([0.0, 0.5, 0.5], zn),
                 ([0.5, 0.0, 0.5], zn), ([0.5, 0.5, 0.0], zn)]
        pos, types = [], []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for frac, t in basis:
                        pos.append(Vector3d((i + frac[0]) * a + d - L / 2,
                                            (j + frac[1]) * a + d - L / 2,
                                            (k + frac[2]) * a + d - L / 2))
                        types.append(t)
        return pos, types, L

    def run_stage(material, work_dir, ensemble, seed):
        return MolecularDynamics(
            work_dir=str(work_dir), ensamble=ensemble, material=material,
            dt_global=2e-15, tau_thermostat=5e-14, tau_barostat=5e-13,
            temperature_initial=TEMPERATURE, temperature_thermostat=TEMPERATURE,
            pressure_barostat=1e-10 * Units.GPa, dump_cycles=50_000,
            slurm_plt=True, show_graphs=False, plot_graphs=False,
            fatal_conservation=1e3,
            cycles_thermalization=args.cycles_thermalization,
            plot_freq=args.production_cycles,
            cycles_end=args.cycles_thermalization + args.production_cycles,
            neighbors=True, r_skin=1.2, cycles_neigbors=20, gr_stride=100,
            random_seed=seed, sample_stride=args.sample_stride,
        ).run()

    # ---- stage 1: NPT at zero pressure ------------------------------------
    pos, types, L = l12_lattice(A0_ANGSTROM * Units.angstrom, NCELL)
    mat = Material(potential=potentials, atomic_weight=weights, atom_name=name, eam=eams,
                   lattice=pos, atom_type=types, box_length=(L, L, L))
    npt = run_stage(mat, out / "npt", MDDriver.NPT, args.npt_velocity_seed)

    # ---- hand-off: last NPT snapshot, rescaled to the mean NPT volume ------
    # (identical to the procedure used for the manuscript runs)
    dumps = sorted(glob.glob(str(out / "npt" / "dump" / "*.xyz")),
                   key=lambda p: int(os.path.splitext(os.path.basename(p))[0]))
    lines = Path(dumps[-1]).read_text().splitlines()
    header = lines[1].split()
    a_val = float(header[header.index("a=") + 1])
    ni = header.index("n=") + 1
    box0 = [a_val * float(header[ni + k]) * Units.angstrom for k in range(3)]
    scale = float(npt.series.V_V0.mean) ** (1.0 / 3.0)
    box = [b * scale for b in box0]

    def wrap(x, length):
        while x >= 0.5 * length:
            x -= length
        while x < -0.5 * length:
            x += length
        return x

    pos_nvt, types_nvt = [], []
    for line in lines[2:]:
        f = line.split()
        if len(f) < 7:
            continue
        xyz = [float(v) * Units.angstrom * scale for v in f[1:4]]
        pos_nvt.append(Vector3d(*(wrap(x, b) for x, b in zip(xyz, box))))
        types_nvt.append(species_to_type[f[0]])
    mat_nvt = Material(potential=potentials, atomic_weight=weights, atom_name=name, eam=eams,
                       lattice=pos_nvt, atom_type=types_nvt, box_length=tuple(box))

    # ---- stage 2: NVT elasticity --------------------------------------------
    nvt = run_stage(mat_nvt, out / "nvt", MDDriver.NVT, args.nvt_velocity_seed)

    t = nvt.series.thermal_results
    result = {
        "benchmark": "l12_cuzn3_5K",
        "interpolation": args.interpolation,
        "cycles_thermalization": args.cycles_thermalization,
        "production_cycles": args.production_cycles,
        "npt_velocity_seed": args.npt_velocity_seed,
        "nvt_velocity_seed": args.nvt_velocity_seed,
        "sample_stride": args.sample_stride,
        "npt_L_L0": float(npt.series.L_L0.mean),
        "nvt_C11_GPa": float(t.elasticity_tensor_11 / Units.GPa),
        "nvt_C12_GPa": float(t.elasticity_tensor_12 / Units.GPa),
        "nvt_C44_GPa": float(t.elasticity_tensor_44 / Units.GPa),
    }
    (out / "result.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
