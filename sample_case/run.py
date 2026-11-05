#!/usr/bin/env python3
"""Run the manuscript's 6 at.% Zn, 300 K NPT-to-NVT workflow."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    package_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=package_root / "build" / "Release",
        help="directory containing the compiled md_driver_module",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "sample_case_output",
        help="new directory for the NPT and NVT outputs",
    )
    parser.add_argument("--cycles-thermalization", type=int, default=50_000)
    parser.add_argument("--production-cycles", type=int, default=100_000)
    parser.add_argument("--placement-seed", type=int, default=42)
    parser.add_argument("--npt-velocity-seed", type=int, default=0)
    parser.add_argument("--nvt-velocity-seed", type=int, default=0)
    return parser.parse_args()


def build_fcc_lattice(Vector3d, Units, a_lattice: float, side: int):
    box_length = side * a_lattice
    basis = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ]
    )
    offset = np.full(3, 0.25)
    positions = []
    for i in range(side):
        for j in range(side):
            for k in range(side):
                for site in basis:
                    xyz = -box_length / 2.0 + (
                        np.array([i, j, k]) + offset + site
                    ) * a_lattice
                    positions.append(Vector3d(*xyz))
    return positions, box_length


def make_material(Material, elements, potentials, eams, positions, atom_types, length):
    return Material(
        potential=potentials,
        atomic_weight=[element[2] for element in elements],
        atom_name="-".join(element[0] for element in elements),
        eam=eams,
        lattice=positions,
        atom_type=atom_types,
        box_length=(length, length, length),
    )


def run_stage(
    MolecularDynamics,
    MDDriver,
    Units,
    material,
    work_dir: Path,
    ensemble,
    thermalization: int,
    production: int,
    velocity_seed: int,
):
    return MolecularDynamics(
        work_dir=str(work_dir),
        ensamble=ensemble,
        material=material,
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
        cycles_thermalization=thermalization,
        plot_freq=production,
        cycles_end=thermalization + production,
        neighbors=True,
        r_skin=1.2,
        cycles_neigbors=20,
        gr_stride=100,
        random_seed=velocity_seed,
    ).run()


def stage_summary(md, Units, include_elasticity: bool) -> dict[str, float]:
    configurational_samples = 1.0 / np.asarray(
        md.series.temperature_conf_inverse.data[md.cycles_thermalization :]
    )
    result = {
        "box_length_A": float(md.Lx * md.series.L_L0.mean / Units.angstrom),
        "pressure_GPa": float(md.series.pressure.mean / Units.GPa),
        "temperature_kinetic_K": float(md.series.temperature_kin.mean),
        "temperature_configurational_reciprocal_mean_K": float(
            1.0 / md.series.temperature_conf_inverse.mean
        ),
        "temperature_configurational_sample_mean_K": float(
            configurational_samples.mean()
        ),
    }
    if include_elasticity:
        result.update(
            {
                "C11_GPa": float(
                    md.series.thermal_results.elasticity_tensor_11 / Units.GPa
                ),
                "C12_GPa": float(
                    md.series.thermal_results.elasticity_tensor_12 / Units.GPa
                ),
                "C44_GPa": float(
                    md.series.thermal_results.elasticity_tensor_44 / Units.GPa
                ),
            }
        )
    return result


def write_decimated_series(path: Path, md, Units, include_elasticity: bool) -> None:
    time = np.asarray(md.series.time)
    stride = max(1, len(time) // 300)
    indices = np.arange(0, len(time), stride)
    if indices[-1] != len(time) - 1:
        indices = np.append(indices, len(time) - 1)

    columns = {
        "step": indices,
        "time_ps": time[indices] / Units.ps,
        "temperature_kinetic_K": np.asarray(md.series.temperature_kin.data)[indices],
        "temperature_configurational_K": (
            1.0 / np.asarray(md.series.temperature_conf_inverse.data)[indices]
        ),
        "pressure_GPa": np.asarray(md.series.pressure.data)[indices] / Units.GPa,
        "L_over_L0": np.asarray(md.series.L_L0.data)[indices],
    }
    if include_elasticity:
        running_means = {}
        for name, series in (
            ("C11", md.series.elasticity_tensor_11),
            ("C12", md.series.elasticity_tensor_12),
            ("C44", md.series.elasticity_tensor_44),
        ):
            values = np.full(len(time), np.nan)
            cumulative = np.asarray(series.mean_cum) / Units.GPa
            values[-len(cumulative) :] = cumulative
            running_means[name] = values[indices]
        columns.update(
            {
                "C11_running_mean_GPa": running_means["C11"],
                "C12_running_mean_GPa": running_means["C12"],
                "C44_running_mean_GPa": running_means["C44"],
            }
        )

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def main() -> int:
    args = parse_args()
    package_root = Path(__file__).resolve().parents[1]
    build_dir = args.build_dir.resolve()
    output_dir = args.output_dir.resolve()
    python_dir = package_root / "src" / "python"

    if not build_dir.is_dir():
        raise SystemExit(f"build directory not found: {build_dir}")
    if output_dir.exists():
        raise SystemExit(f"output directory already exists: {output_dir}")
    if args.cycles_thermalization < 1 or args.production_cycles < 1:
        raise SystemExit("both cycle counts must be positive")

    sys.path.insert(0, str(build_dir))
    sys.path.insert(0, str(python_dir))

    from md_driver_module import MDDriver, Units, Vector3d
    from molecular_dynamics import Material, MolecularDynamics
    import get_potential

    elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy")
    positions, initial_length = build_fcc_lattice(
        Vector3d, Units, 3.6345 * Units.angstrom, side=5
    )
    rng = np.random.RandomState(args.placement_seed)
    atom_types = [1] * len(positions)
    for index in rng.choice(len(positions), size=30, replace=False):
        atom_types[int(index)] = 0

    npt_material = make_material(
        Material,
        elements,
        potentials,
        eams,
        positions,
        atom_types,
        initial_length,
    )
    npt = run_stage(
        MolecularDynamics,
        MDDriver,
        Units,
        npt_material,
        output_dir / "npt",
        MDDriver.NPT,
        args.cycles_thermalization,
        args.production_cycles,
        args.npt_velocity_seed,
    )

    # Reproduce the manuscript handoff exactly: read the final NPT trajectory
    # snapshot, scale it by the production mean (V/V0)^(1/3), and use that
    # production-mean box for NVT.
    dump_paths = sorted(
        (output_dir / "npt" / "dump").glob("*.xyz"),
        key=lambda path: int(path.stem),
    )
    dump_lines = dump_paths[-1].read_text().splitlines()
    species_to_type = {element[0]: index for index, element in enumerate(elements)}
    nvt_positions = []
    nvt_atom_types = []
    for line in dump_lines[2:]:
        fields = line.split()
        if len(fields) < 7:
            continue
        nvt_positions.append(
            Vector3d(*(float(value) * Units.angstrom for value in fields[1:4]))
        )
        nvt_atom_types.append(species_to_type[fields[0]])

    scale = npt.series.V_V0.mean ** (1.0 / 3.0)
    mean_length = initial_length * scale
    nvt_positions = [
        Vector3d(position.x * scale, position.y * scale, position.z * scale)
        for position in nvt_positions
    ]
    nvt_material = make_material(
        Material,
        elements,
        potentials,
        eams,
        nvt_positions,
        nvt_atom_types,
        mean_length,
    )
    nvt = run_stage(
        MolecularDynamics,
        MDDriver,
        Units,
        nvt_material,
        output_dir / "nvt",
        MDDriver.NVT,
        args.cycles_thermalization,
        args.production_cycles,
        args.nvt_velocity_seed,
    )

    result = {
        "case": "6 at.% Zn alpha-brass at 300 K, 500 atoms",
        "cycles_thermalization": args.cycles_thermalization,
        "production_cycles": args.production_cycles,
        "placement_seed": args.placement_seed,
        "npt_velocity_seed": args.npt_velocity_seed,
        "nvt_velocity_seed": args.nvt_velocity_seed,
        "npt": stage_summary(npt, Units, include_elasticity=False),
        "nvt": stage_summary(nvt, Units, include_elasticity=True),
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_decimated_series(output_dir / "npt_series_decimated.csv", npt, Units, False)
    write_decimated_series(output_dir / "nvt_series_decimated.csv", nvt, Units, True)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
