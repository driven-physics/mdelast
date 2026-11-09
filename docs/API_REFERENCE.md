# API reference

## Python simulation API

### `molecular_dynamics.Material`

Container for pair potentials, EAM density/embedding functions, atomic
weights, species names, and lattice information. Pair potentials, atomic
weights, and non-None EAM objects are normalized to lists. For generated FCC structures set `a_lattice`; for an
explicit structure set `lattice`, `atom_type`, and `box_length`.

### `molecular_dynamics.MolecularDynamics`

The high-level simulation controller. Construction validates inputs, creates
the work and dump directories, builds or accepts a lattice, initializes
velocities, and creates the compiled driver. `run()` advances the simulation,
samples observables, performs final averages, writes `series.h5`, and
optionally plots. `save_pickle()` remains available for historical workflows
but is not called by the included run scripts.

The NPT elastic-constant estimator assumes cubic symmetry. For a non-cubic
crystal, use NVT at a separately established equilibrium cell shape.

Supporting methods:

- `build_fcc_lattice()`: construct the requested periodic FCC supercell
- `set_initial_velocities()`: initialize zero-momentum thermal velocities
- `dump()`: write one XYZ trajectory frame
- `perform_thermal_average()`: calculate ensemble-dependent final properties
- `save_h5()` and `save_pickle()`: serialize results
- `plot()`: invoke `PlotSeries`

Inputs and units are detailed in `USAGE.md`.

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

Potential files are resolved from the top-level `potentials/` directory. The
single-element `get_EAM_file` and Lennard-Jones helper are retained for
compatibility.

### Table interpolation

The tabulated functions φ(r), ρ(r) and F(ρ) are evaluated with one of two
schemes, chosen by the `interpolation` argument of `get_EAM_alloy_file` and
`get_EAM_file`:

| value | scheme |
|---|---|
| `"quintic"` (default) | one quintic Hermite polynomial per grid interval, built from f, f′, f″ at its two nodes; f′ and f″ are the exact derivatives of f and are continuous across nodes |
| `"bspline"` | three independent cubic B-splines through the node values of f, f′ and f″ — the scheme of mdelast ≤ 1.1, kept to reproduce version 1.1 with the same build and inputs |

`interpolation=None` uses the module default `get_potential.DEFAULT_INTERPOLATION`
(`"quintic"`); assigning it changes the default for every later load in the
process. Any other value raises `ValueError`.

```python
elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy", interpolation="bspline")
material = get_potential.get_EAM_file(file_name="Cu01.eam.alloy.txt", interpolation="bspline")
```

The `"bspline"` scheme needs at least 5 grid points per table.

## Time-series API

- `ThermalResults`: dynamic namespace for final aggregate properties
- `TimeSeries`: scalar samples, population mean/std, cumulative statistics,
  Gaussian curve data, and HDF5 serialization
- `TimeSeriesVector(N)`: `N` component series
- `TimeSeriesMatrix(N)`: `N x N` component series
- `Series`: all sampled temperatures, energies, pressure/stress, box changes,
  Born tensors, and final thermal results

`write_h5(file_h5, group_name)` stores each available field below the supplied
group. Vector and matrix components use one-based suffixes in HDF5 names.

## Plotting API

`plot.PlotSeries` consumes a `Series` and writes diagnostic plots. Public plot
groups are:

- `temperature`
- `pressure`
- `elasticity`
- `volume_fluctuations`
- `conservation`
- `energies_cv`
- `partial_pair_correlation`

`plot_splitted_x` and `get_splitted_axes` support broken-axis figures.

## VACF and VDOS API

The documented functions in `vacf.py` provide post-processing independent of
the run loop:

- `read_dump_xyz`: read positions, velocities, species, times, and box data
- `compute_vacf`: mass-weighted velocity autocorrelation
- `compute_vdos`: Fourier-transform VACF to a vibrational density of states
- `compute_vdos_moments`: spectral moments
- `plot_vacf` and `plot_vdos`: publication figures
- `write_spectral_moments` and `save_results_npz`: persist results

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
  (default) or `"bspline"` (see *Table interpolation*); the choice is readable
  as the `.interpolation` attribute.
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

The low-level driver exposes state arrays, energies, temperature estimators,
pressure/stress, Born stress and elasticity tensors, pair correlations,
current time/cycle, atom types, and box state as read-only properties. Its
public operations are:

- `integrate()`
- `calculate_neighbors()`
- `calculate_pair_correlation()`
- `calculate_pair_correlations()`
- `calculate_partial_pair_correlation()`

One property is read-write: `sample_observables` (default True). When False, the
next `integrate()` skips the Born elasticity tensor and the Laplacian, which keep
their last values. `total_laplacian` and `total_force2` hold the two sums of the
configurational-temperature ratio from the last sampled step.

Most users should construct `MolecularDynamics` instead of instantiating
`MDDriver` directly; the latter has a long, order-sensitive constructor and
expects all values already converted to internal units.
