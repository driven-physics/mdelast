# API reference

## Python simulation API

### `molecular_dynamics.Material`

Holds pair potentials, EAM density and embedding functions, atomic weights,
species names, and the lattice. Pair potentials, atomic weights, and non-None
EAM objects are stored as lists. For a generated FCC structure set `a_lattice`.
For an explicit structure set `lattice`, `atom_type`, and `box_length`.

### `molecular_dynamics.MolecularDynamics`

Runs a simulation. The constructor checks the inputs, creates the work and
dump directories, builds or takes a lattice, sets the velocities, and creates
the compiled driver. `run()` integrates, samples the observables, averages,
writes `series.h5`, and plots if asked. `save_pickle()` is kept for old
workflows. The run scripts do not call it.

The NPT elastic-constant estimator assumes cubic symmetry. For a non-cubic
crystal, use NVT at a separately established equilibrium cell shape.

Other methods:

- `build_fcc_lattice()`: build the periodic FCC supercell
- `set_initial_velocities()`: thermal velocities with zero total momentum
- `dump()`: write one XYZ trajectory frame
- `perform_thermal_average()`: final properties of the ensemble
- `save_h5()` and `save_pickle()`: write results
- `plot()`: call `PlotSeries`

`USAGE.md` lists the inputs and units.

### Material helpers

- `get_EAM_file(*, file_name, a_factor=1.0)`: single-element SETFL/FS loader, default
  tables only. Use `get_potential.get_EAM_file` to choose the interpolation.
- `get_Cu_Morse(a_factor=1.0)`: analytic Cu Morse material
- `get_Cu_Table_Morse()`: unavailable in version 1.2. The helper raises
  `AttributeError` before it constructs the table.
- `get_Ar_LJ(epsilon=None, sigma=None, rcut_fac=3.5)`: Ar Lennard-Jones material
- `get_deriv(x, y)`: finite differences at interval midpoints
- `get_derivatives(x, y)`: quintic-spline first and second derivatives

## Potential-file API

`get_potential.get_EAM_alloy_file(file_name, interpolation=None)` reads a
multi-element LAMMPS SETFL alloy file and returns:

1. element records `[symbol, atomic_number, atomic_weight, lattice_constant]`;
2. pair potentials in triangular pair ordering;
3. one `EAM` object per species.

Potential files are read from the top-level `potentials/` directory. The
single-element `get_EAM_file` and the Lennard-Jones helper are kept for
compatibility.

### Table interpolation

The `interpolation` argument of `get_EAM_alloy_file` and `get_EAM_file` selects
how the tabulated functions φ(r), ρ(r) and F(ρ) are evaluated:

| value | scheme |
|---|---|
| `"quintic"` (default) | one quintic Hermite polynomial per grid interval, built from f, f′, f″ at its two nodes. f′ and f″ are the exact derivatives of f and are continuous across nodes |
| `"bspline"` | three independent cubic B-splines through the node values of f, f′ and f″. This is the scheme of mdelast ≤ 1.1, kept to reproduce version 1.1 with the same build and inputs |

`interpolation=None` uses the module default `get_potential.DEFAULT_INTERPOLATION`
(`"quintic"`). Assigning to it changes the default for every later load in the
process. Any other value raises `ValueError`.

```python
elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy", interpolation="bspline")
material = get_potential.get_EAM_file(file_name="Cu01.eam.alloy.txt", interpolation="bspline")
```

The `"bspline"` scheme needs at least 5 grid points per table.

## Time-series API

- `ThermalResults`: namespace for the final properties
- `TimeSeries`: scalar samples, mean and standard deviation, cumulative
  statistics, Gaussian curve data, and HDF5 output
- `TimeSeriesVector(N)`: `N` component series
- `TimeSeriesMatrix(N)`: `N x N` component series
- `Series`: all sampled temperatures, energies, pressure/stress, box changes,
  Born tensors, and final thermal results

`write_h5(file_h5, group_name)` stores each field below the given group.
Vector and matrix components have one-based suffixes in the HDF5 names.

## Plotting API

`plot.PlotSeries` takes a `Series` and writes diagnostic plots. The plot
groups are:

- `temperature`
- `pressure`
- `elasticity`
- `volume_fluctuations`
- `conservation`
- `energies_cv`
- `partial_pair_correlation`

`plot_splitted_x` and `get_splitted_axes` make broken-axis figures.

## VACF and VDOS API

`vacf.py` post-processes trajectory dumps, independent of the run loop:

- `read_dump_xyz`: read positions, velocities, species, times, and box data
- `compute_vacf`: mass-weighted velocity autocorrelation
- `compute_vdos`: vibrational density of states, the Fourier transform of the VACF
- `compute_vdos_moments`: spectral moments
- `plot_vacf` and `plot_vdos`: figures
- `write_spectral_moments` and `save_results_npz`: write results

## Compiled `md_driver_module` API

### Geometry and units

- `Vector3d(x, y, z)`: mutable components plus `length()` and `length2()`
- `Box(h)`: cell matrix with `volume()` and `length()`
- `Units`: CGS conversion constants `k_boltz`, `Navogadro`, `ev`, `joule`,
  `angstrom`, `nm`, `bar`, `kbar`, `atm`, `Pa`, `GPa`, `mus`, `ns`, `ps`,
  and `fs`

### Potentials

- `Potential`: base class exposing `rcut` and `calculate(r)`
- `PotentialValue`: returned `V`, `dVdr`, and `d2Vdr2`
- `Morse(rcut, D, alpha, r0)`
- `LennardJones(rcut, epsilon, sigma)`
- `TablePotential(rcut, rmin_fatal, dr, V, dVdr, d2Vdr2[, interpolation])`:
  tabulated pair potential with two derivatives
- `EAM(rcut, rhomax_fatal, dr, drho, rho, drhodr, d2rhodr2, F, dFdrho, d2Fdrho2[, interpolation])`:
  tabulated density ρ(r) and embedding function F(ρ) with their derivatives

  For both, the optional last positional argument is the string `"quintic"`
  (default) or `"bspline"` (see *Table interpolation*). The `.interpolation`
  attribute holds the choice.
- `PotentialVector` and `EAMVector`: pybind11 vector containers

### Extended-system controls

- `Thermostat`: Nose-Hoover chain state
- `Barostat`: isotropic barostat state
- `BarostatFlex`: placeholder, no equations of motion
- `MDDriver.Ensamble`: `NVE`, `NVT`, `NPT`. `NPT_FLEX` is reserved and not
  implemented. Selecting it aborts the process.
- `MDDriver.ThermostatCoupling`: `THERMOSTAT_KINETIC` (default).
  `THERMOSTAT_CONFIGURATIONAL` and `THERMOSTAT_AVERAGE` are experimental (see `USAGE.md`).

### `MDDriver`

The driver has read-only properties for the state arrays, energies,
temperature estimators, pressure and stress, Born stress and elasticity
tensors, pair correlations, current time and cycle, atom types, and box. Its
methods are:

- `integrate()`
- `calculate_neighbors()`
- `calculate_pair_correlation()`
- `calculate_pair_correlations()`
- `calculate_partial_pair_correlation()`

One property is read-write: `sample_observables` (default True). When False, the
next `integrate()` skips the Born elasticity tensor and the Laplacian, which keep
their last values. `total_laplacian` and `total_force2` hold the two sums of the
configurational-temperature ratio from the last sampled step.

Use `MolecularDynamics` unless you need the driver itself. The `MDDriver`
constructor takes many positional arguments, all in internal units.
