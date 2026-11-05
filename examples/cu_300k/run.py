#!/usr/bin/env python3
"""Run the publication example: pure FCC Cu at 300 K."""

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
        help="directory containing md_driver_module",
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
        default=package_root / "validation_outputs" / "cu_300k",
    )
    parser.add_argument("--cycles-thermalization", type=int, default=50_000)
    parser.add_argument("--production-cycles", type=int, default=100_000)
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
    from molecular_dynamics import MolecularDynamics, get_EAM_file

    total_cycles = args.cycles_thermalization + args.production_cycles
    md = MolecularDynamics(
        work_dir=str(output_dir),
        ensamble=MDDriver.NPT,
        material=get_EAM_file(file_name="Cu01.eam.alloy.txt"),
        N_lattice_side=5,
        dt_global=2e-15,
        tau_thermostat=5e-14,
        tau_barostat=5e-13,
        temperature_initial=300.0,
        temperature_thermostat=300.0,
        pressure_barostat=1e-10 * Units.GPa,
        dump_cycles=50_000,
        slurm_plt=False,
        show_graphs=False,
        plot_graphs=False,
        fatal_conservation=1e3,
        cycles_thermalization=args.cycles_thermalization,
        plot_freq=args.production_cycles,
        cycles_end=total_cycles,
        random_seed=args.velocity_seed,
    ).run()

    result = {
        "case": "cu_300k",
        "build_dir": str(build_dir),
        "python_dir": str(python_dir),
        "atoms": 500,
        "temperature_kinetic_K": float(md.series.temperature_kin.mean),
        "temperature_configurational_K": float(1.0 / md.series.temperature_conf_inverse.mean),
        "pressure_GPa": float(md.series.pressure.mean / Units.GPa),
        "bulk_modulus_GPa": float(md.series.thermal_results.bulk_modulus / Units.GPa),
        "C11_GPa": float(md.series.thermal_results.elasticity_tensor_11 / Units.GPa),
        "C12_GPa": float(md.series.thermal_results.elasticity_tensor_12 / Units.GPa),
        "C44_GPa": float(md.series.thermal_results.elasticity_tensor_44 / Units.GPa),
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
