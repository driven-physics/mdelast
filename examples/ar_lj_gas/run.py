#!/usr/bin/env python3
"""Run a short Argon Lennard-Jones gas simulation (500 atoms, NVT at 300 K).

This is a minimal example demonstrating:
  - Setting up a monatomic LJ gas from the built-in Ar potential
  - Running a short NVT equilibration followed by a production phase
  - Collecting basic thermodynamic output (temperature, pressure, energy)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    package_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=package_root / "build" / "Release",
        help="directory containing the compiled md_driver_module",
    )
    parser.add_argument(
        "--python-dir",
        type=Path,
        default=package_root / "src" / "python",
        help="directory containing the Python MD modules",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "validation_outputs" / "ar_lj_gas",
        help="directory for simulation output files",
    )
    parser.add_argument("--cycles-thermalization", type=int, default=5_000)
    parser.add_argument("--production-cycles", type=int, default=10_000)
    parser.add_argument("--velocity-seed", type=int, default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    package_root = Path(__file__).resolve().parents[2]
    build_dir = args.build_dir.resolve()
    python_dir = args.python_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not build_dir.is_dir():
        raise SystemExit(f"build directory not found: {build_dir}")
    if not python_dir.is_dir():
        raise SystemExit(f"Python source directory not found: {python_dir}")
    if output_dir.exists():
        raise SystemExit(f"output directory already exists: {output_dir}")

    sys.path.insert(0, str(build_dir))
    sys.path.insert(0, str(python_dir))

    from md_driver_module import MDDriver, Units
    from molecular_dynamics import MolecularDynamics, get_Ar_LJ

    # ---- Material ----
    # Standard Lennard-Jones argon: epsilon = 0.0104 eV, sigma = 3.405 A.
    material = get_Ar_LJ()

    # ---- Lattice geometry ----
    # FCC lattice with side=5 gives 4*5^3 = 500 atoms.
    # A lattice constant of 8 A (vs solid Ar ~5.3 A) yields a gas-like
    # reduced density rho* = N*sigma^3/V ~ 0.3 at supercritical conditions.
    a_lattice = 8.0 * Units.angstrom
    n_side = 5  # 500 atoms

    # ---- Simulation parameters ----
    total_cycles = args.cycles_thermalization + args.production_cycles

    md = MolecularDynamics(
        work_dir=str(output_dir),
        ensamble=MDDriver.NVT,
        material=material,
        N_lattice_side=n_side,
        a_lattice=a_lattice,
        dt_global=5e-15,
        tau_thermostat=1e-13,
        tau_barostat=5e-13,
        temperature_initial=300.0,
        temperature_thermostat=300.0,
        pressure_barostat=0.0,
        dump_cycles=5_000,
        slurm_plt=False,
        show_graphs=False,
        plot_graphs=True,
        fatal_conservation=1e3,
        cycles_thermalization=args.cycles_thermalization,
        plot_freq=args.production_cycles,
        cycles_end=total_cycles,
        neighbors=True,
        r_skin=1.5,
        cycles_neigbors=10,
        gr_stride=50,
        random_seed=args.velocity_seed,
    ).run()

    # ---- Collect results ----
    result = {
        "case": "ar_lj_gas",
        "build_dir": str(build_dir),
        "python_dir": str(python_dir),
        "atoms": 500,
        "temperature_kinetic_K": float(md.series.temperature_kin.mean),
        "temperature_configurational_K": float(
            1.0 / md.series.temperature_conf_inverse.mean
        ),
        "pressure_GPa": float(md.series.pressure.mean / Units.GPa),
        "velocity_seed": args.velocity_seed,
        "cycles_thermalization": args.cycles_thermalization,
        "production_cycles": args.production_cycles,
    }

    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
