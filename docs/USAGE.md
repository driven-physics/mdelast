# Usage

## Simulation workflow

1. Build `md_driver_module`.
2. Add the selected build directory and `src/python` to `sys.path`.
3. create a `Material`;
4. construct `MolecularDynamics` with keyword arguments;
5. call `.run()`.

The scripts under `examples/` implement this workflow without machine-specific
paths.

## Material definitions

For a single-component FCC material, provide a pair potential, atomic weight,
name, and lattice constant. `MolecularDynamics` generates a periodic FCC
supercell when `N_lattice_side` is supplied.

For an alloy or custom structure, provide:

- one pair potential for each triangular species pair;
- one atomic weight and EAM object per species;
- an explicit list of `Vector3d` positions;
- an integer species index for every atom;
- `(Lx, Ly, Lz)` as `box_length`.

`get_Ar_LJ()` creates a single-species Lennard-Jones argon material (no EAM,
no lattice constant — supply `a_lattice` and `N_lattice_side` to
`MolecularDynamics`).
`get_EAM_file("Cu01.eam.alloy.txt")` creates the included pure-Cu material.
`get_EAM_alloy_file("ZnCu.eam.alloy")` returns element metadata, pair
potentials, and EAM objects for an explicit alloy `Material`.

## `MolecularDynamics` parameters

The following keyword arguments are required:

- `material`: a `Material` instance
- `ensamble`: `MDDriver.NVE`, `MDDriver.NVT`, or `MDDriver.NPT`
- `temperature_initial`: initial velocity temperature in K
- `temperature_thermostat`: thermostat target in K
- `tau_thermostat`: thermostat time constant in seconds
- `pressure_barostat`: NPT pressure target in CGS pressure units
- `tau_barostat`: barostat time constant in seconds
- `dt_global`: integration time step in seconds
- `cycles_thermalization`: samples excluded from final averages
- `plot_freq`: production/plot interval used by the run loop
- `cycles_end`: total number of integration cycles
- `dump_cycles`: trajectory dump interval

Common optional arguments:

- `N_lattice_side`: FCC cells per direction for generated lattices
- `nchains_thermostat`, `nchains_barostat`: Nose-Hoover chain lengths
- `thermostat_coupling`: kinetic, configurational, or average coupling
- `neighbors`, `r_skin`, `cycles_neigbors`: neighbor-list controls
- `gr_stride`: post-thermalization radial-distribution sampling interval
- `initial_velocities`: explicit velocity vectors
- `random_seed`: deterministic velocity seed
- `plot_graphs`, `show_graphs`, `slurm_plt`: plotting controls
- `fatal_conservation`: energy-conservation guard

The misspellings `ensamble` and `cycles_neigbors` are stable historical API
names and must be used exactly as shown.

## Included examples

### Ar LJ gas at 300 K

```bash
python3 examples/ar_lj_gas/run.py \
  --build-dir build/Release \
  --output-dir outputs/ar_lj_gas \
  --velocity-seed 0
```

A 5 x 5 x 5 FCC supercell (500 atoms) of argon interacting via the standard
Lennard-Jones potential (epsilon = 0.0104 eV, sigma = 3.405 A). The lattice
constant is 8.0 A — well above the solid equilibrium — giving a gas-like
reduced density rho* ~ 0.3 at supercritical conditions. The simulation runs in
NVT at 300 K with a 5 fs timestep.

This is the fastest example (5,000 thermalization + 10,000 production steps by
default) and is suitable as a quick installation check.

### Pure Cu at 300 K

```bash
python3 examples/cu_300k/run.py \
  --build-dir build/Release \
  --output-dir outputs/cu_300k \
  --velocity-seed 0
```

This is a 5 x 5 x 5 FCC supercell (500 atoms) using
`Cu01.eam.alloy.txt` in NPT at effectively zero external pressure.

### 6 at.% Zn at 300 K

```bash
python3 examples/zn06_300k/run.py \
  --build-dir build/Release \
  --output-dir outputs/zn06_300k \
  --placement-seed 42 \
  --velocity-seed 0
```

Thirty of the 500 FCC sites are selected as Zn, giving exactly 6 at.% Zn. The
example uses `ZnCu.eam.alloy`, NPT, and the production lattice constant
3.6345 Angstrom.

All examples accept `--cycles-thermalization` and `--production-cycles`.
The Cu and Zn defaults (50,000 and 100,000) are the article-scale schedule;
the Ar LJ gas defaults (5,000 and 10,000) provide a quick demonstration.
`--python-dir` can select another copy of the Python wrapper for regression
testing; normal users should keep its package-relative default.

The `zn06_300k` example is an NPT-only demonstration. The manuscript
regression is the two-stage `sample_case/run.py` workflow. After its full run,
apply the shipped tolerances with:

```bash
python3 sample_case/validate.py sample_case_output/result.json
```

The command exits nonzero if metadata, box lengths, temperatures, or NVT
elastic constants fail their acceptance checks.

## Output files

Each output directory contains:

- `series.h5`: sampled data, cumulative statistics, and thermal results
- `dump/*.xyz`: OVITO-readable trajectory snapshots
- `result.json`: concise values and seeds emitted by the example
- `figs/`: optional diagnostic figures when plotting is enabled

The library retains a `save_pickle()` method for historical workflows, but the
included run scripts do not call it.

Important HDF5 groups include:

- `series/<observable>/data`, `mean`, `std`, and cumulative statistics
- component series such as `stress_tensor_1` and
  `elasticity_tensor_1_1`
- `thermal_results/` for final bulk, elastic, radial-distribution, and
  consistency quantities
