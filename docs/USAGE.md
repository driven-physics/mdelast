# Usage

## Workflow

1. Build `md_driver_module`.
2. Add the build directory and `src/python` to `sys.path`.
3. Create a `Material`.
4. Construct `MolecularDynamics` with keyword arguments.
5. Call `.run()`.

The scripts under `examples/` do this.

## Materials

For a single-component FCC material, give a pair potential, atomic weight,
name, and lattice constant. `MolecularDynamics` builds a periodic FCC supercell
when `N_lattice_side` is given.

For an alloy or another structure, give:

- one pair potential for each triangular species pair
- one atomic weight and EAM object per species
- a list of `Vector3d` positions
- an integer species index for every atom
- `(Lx, Ly, Lz)` as `box_length`

`get_Ar_LJ()` returns a single-species Lennard-Jones argon material. It has no
EAM and no lattice constant, so pass `a_lattice` and `N_lattice_side` to
`MolecularDynamics`.
`get_potential.get_EAM_file(file_name="Cu01.eam.alloy.txt")` returns the
included pure-Cu material.
`get_EAM_alloy_file("ZnCu.eam.alloy")` returns element metadata, pair
potentials, and EAM objects for an alloy `Material`.

Both loaders of `get_potential` take `interpolation="quintic"` (default) or `"bspline"`.
`"bspline"` selects the cubic B-spline tables of mdelast ≤ 1.1 and reproduces
version 1.1 results with the same build and inputs. See *Table interpolation* in
`API_REFERENCE.md`.

## `MolecularDynamics` parameters

Required keyword arguments:

- `material`: a `Material` instance
- `ensamble`: `MDDriver.NVE`, `MDDriver.NVT`, or `MDDriver.NPT`.
  In version 1.2 an NVE run integrates correctly but `run()` raises at the final
  averaging step. Drive `MDDriver.integrate()` directly, as benchmark S does.
  The NPT elastic-constant estimator assumes cubic symmetry. For a non-cubic
  crystal, use NVT at a separately established equilibrium cell shape.
- `temperature_initial`: initial velocity temperature in K
- `temperature_thermostat`: thermostat target in K
- `tau_thermostat`: thermostat time constant in seconds
- `pressure_barostat`: NPT pressure target in CGS pressure units
- `tau_barostat`: barostat time constant in seconds
- `dt_global`: time step in seconds
- `cycles_thermalization`: MD steps excluded from the averages
- `plot_freq`: averages, `series.h5` and plots are written every `plot_freq`
  steps after thermalization. Set it to the production length for one final output.
- `cycles_end`: nominal run length in MD steps. The current loop executes
  through `cycles_end + 1`.
- `dump_cycles`: trajectory dump interval
- `fatal_conservation`: the run stops when the relative error of the conserved
  energy exceeds this value (the scripts use a loose threshold of 1e3)

Optional arguments:

- `N_lattice_side`: FCC cells per direction for generated lattices
- `nchains_thermostat`, `nchains_barostat`: Nose-Hoover chain lengths (10 and 10)
- `thermostat_coupling`: `"kinetic"` (default), the scheme of the paper. The values
  `"configurational"` and `"average"` are experimental. They are not derived or
  validated in the paper and do not sample a known ensemble.
- `neighbors` (default False: all pairs), `r_skin` (list radius as a multiple of
  the cut-off, at least 1, and 1.2 in the paper), `cycles_neigbors` (rebuild
  interval in steps, default 15, and 20 in the paper)
- `gr_stride`: sampling interval of the radial distribution after thermalization
- `sample_stride`: compute and store the observables every N MD steps (default 1).
  On the other steps the Born elasticity tensor and the configurational-temperature
  Laplacian are not computed and no sample is stored. The stress-fluctuation term,
  the Born term, and T_conf are averaged over the same samples. The trajectory
  is unchanged. `cycles_thermalization` and `plot_freq` must be multiples of N, and
  N > 1 requires kinetic thermostat coupling.
- `initial_velocities`: velocity vectors
- `random_seed`: velocity seed
- `plot_graphs`, `show_graphs`, `slurm_plt`: plotting controls
- `work_dir`: output directory (default: the current directory)

Write `ensamble` and `cycles_neigbors` exactly as shown. The spelling is
historical and stable.

## Examples

### Ar LJ gas at 300 K

```bash
python3 examples/ar_lj_gas/run.py \
  --build-dir build/Release \
  --output-dir outputs/ar_lj_gas \
  --velocity-seed 0
```

A 5 x 5 x 5 FCC supercell (500 atoms) of Lennard-Jones argon
(epsilon = 0.0104 eV, sigma = 3.405 A) with a lattice constant of 8.0 A, which
gives a gas-like reduced density rho* ~ 0.3. NVT at 300 K, 5 fs timestep.
This is the fastest example (5,000 thermalization + 10,000 production steps by
default).

### Pure Cu at 300 K

```bash
python3 examples/cu_300k/run.py \
  --build-dir build/Release \
  --output-dir outputs/cu_300k \
  --velocity-seed 0
```

A 5 x 5 x 5 FCC supercell (500 atoms) with `Cu01.eam.alloy.txt`, NPT at zero
external pressure.

### 6 at.% Zn at 300 K

```bash
python3 examples/zn06_300k/run.py \
  --build-dir build/Release \
  --output-dir outputs/zn06_300k \
  --placement-seed 42 \
  --velocity-seed 0
```

Thirty of the 500 FCC sites are Zn, exactly 6 at.% Zn. The example uses
`ZnCu.eam.alloy`, NPT, and the lattice constant 3.6345 Angstrom.

All examples accept `--cycles-thermalization` and `--production-cycles`.
The Cu and Zn defaults (50,000 and 100,000) are the schedule of the paper. The
Ar defaults are 5,000 and 10,000. `--python-dir` selects another copy of the
Python layer for regression tests.

The `zn06_300k` example runs NPT only. The two-stage calculation of the paper
is `sample_case/run.py`. After its full run, check it with:

```bash
python3 sample_case/validate.py sample_case_output/result.json
```

The command exits nonzero if metadata, box lengths, temperatures, or NVT
elastic constants fail.

## Output files

Each output directory contains:

- `series.h5`: sampled data, cumulative statistics, and thermal results
- `dump/*.xyz`: trajectory snapshots, readable by OVITO
- `result.json`: summary values and seeds written by the example
- `figs/`: diagnostic figures when plotting is on

`save_pickle()` is kept for old workflows. The run scripts do not call it.

HDF5 groups:

- `series/<observable>/data`, `mean`, `std`, and cumulative statistics
- component series such as `stress_tensor_1` and
  `elasticity_tensor_1_1`
- `thermal_results/` for final bulk, elastic, radial-distribution, and
  consistency quantities
