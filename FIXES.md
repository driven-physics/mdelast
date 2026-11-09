# Fixes in this version

## Version 1.2 (not yet released)

### 4. Table interpolation of the potential

**Files:** `src/poly_table.hpp` (new), `src/potential.hpp`, `src/eam.hpp`,
`src/md_driver_module.cpp`, `src/python/get_potential.py`

**Change.** Up to version 1.1 every tabulated function (phi(r), rho(r), F(rho))
was evaluated with three independent cubic B-splines through the node values
of f, f' and f''; f' and f'' were therefore not the derivatives of the f
interpolant. The default is now one quintic Hermite polynomial per grid
interval, built from f, f' and f'' at its two nodes. f' and f'' are the exact
derivatives of f, f, f' and f'' are continuous across nodes, and one lookup
returns all three. The previous scheme remains available as
`interpolation="bspline"` (`docs/API_REFERENCE.md`, "Table interpolation") and
reproduces version 1.1 bit for bit.

**Effect.** On fixed configurations of the Zn06 benchmark cell the Born tensor
changes by 4e-8, the configurational temperature by 8e-8 and the forces by
5e-6 (relative). Forces now agree with -dU/dx to 4e-8 at a 1e-4 A step
(B-splines: 4e-7). In benchmark S the Born-constant finite-difference residual
falls from 0.013 to 1.0e-4 GPa and the Laplacian residual from 3.7e-6 to
2.9e-8. The fixed-cell values of benchmark S change by up to 3.6e-7
(relative); `benchmarks/expected/smoke.json` lists both sets. MD results
change only through the different chaotic trajectory. A step costs about half
as long (benchmark Z1: 2574 s against 4524 s with the B-spline tables, same
build and node).

### 5. Sampling interval of the observables

**Files:** `src/md_driver.cpp`, `src/md_driver.hpp`, `src/md_driver_module.cpp`,
`src/python/molecular_dynamics.py`

**Change.** New keyword `sample_stride` of `MolecularDynamics` (default 1):
the Born tensor, the configurational-temperature sums and the stored samples
are computed every `sample_stride` steps, while forces and the virial are
computed at every step. The trajectory does not depend on `sample_stride`.

**Effect.** None for the default `sample_stride = 1`, which is bit-identical to
version 1.1. With `sample_stride = 10` the standard errors of the elastic
constants are unchanged (benchmark state: 0.204, 0.155, 0.045 GPa against
0.204, 0.154, 0.045 GPa for C11, C12, C44).

### Reference values

`benchmarks/expected/zn06_300K.json` and `l12_cuzn3_5K.json` now hold the
values of version 1.2 reference runs (default options, same seeds); the
version 1.0 values are kept as `expected_previous` and the tolerances are
unchanged. The version 1.0 values carry the 3N degrees-of-freedom bias of entry
3: Z1 runs of version 1.1 already gave an NPT box ratio L/L0 0.77 tolerances
below the version 1.0 value.
`sample_case/expected_output` keeps the manuscript values of version 1.0;
version 1.2 passes all its checks.

## Version 1.1

### 3. Degrees of freedom of the thermostat and the kinetic temperature

**File:** `src/python/molecular_dynamics.py` (`self.Nf`)

**Problem.** The number of degrees of freedom was `Nf = 3N`, although the three
center-of-mass momentum components are set to zero at initialization and stay
zero. The thermostat held `<2K> = 3N kT` on `3N-3` active momenta, so the atoms
sat at `T (1 + 1/(N-1))`, 0.2% above the target at N = 500, while the reported
kinetic temperature, computed with the same `Nf`, read exactly `T`. The
configurational temperature, which counts no degrees of freedom, showed the
excess.

**Fix.** `Nf = 3(N-1)`, used consistently for the kinetic temperature, the
thermostat force and mass, the conserved energy, and the barostat
(`alpha = 1 + 3/Nf`, `W = Nf kT tau_P^2`).

**Effect.** Verified with 200 velocity-seed replica pairs of the sample case
(500 atoms, 6 at.% Zn, 300 K): `T_conf - T_K` falls by 0.204 +- 0.013% of T
(predicted 0.200%), and the NVT elastic constants change by +0.027% (C11),
+0.015% (C12), and +0.035% (C44), below their statistical uncertainties. The
reference outputs in `sample_case/expected_output` and `benchmarks/expected`
were produced with version 1.0; version 1.1 passes every check within the
tolerances.

## Version 1.0

Two bugs are fixed relative to the previous release. Neither changes any NVT result or any two-element (Cu–Zn) result.

## 1. NPT elastic constants: ideal-gas pressure in the fluctuating stress

**File:** `src/python/molecular_dynamics.py` (`perform_thermal_average`, the elasticity block)

**Problem.** In NPT the code evaluated

    C = <C^B> - <V>/kT cov(sigma^B, sigma^B) + <N/V> kT (dd + dd) + K dd

The Born stress `sigma^B` is configurational only. In NPT the ideal-gas pressure `N kT / V` fluctuates with the volume, so it has to stay inside the fluctuating stress. Using `sigma^B` together with the `+K` term double-counts it. The resulting bias is `+(2 p_id - p_id^2/K)` on C11, C22, C33, C12, C13 and C23. That is about +0.67 GPa at 300 K and about +1.7 GPa at 800 K for Cu–Zn. As a check, for an ideal gas the old formula gives C11 = 3P and C12 = P instead of 2P and 0.

**Fix.** In NPT the fluctuation term now uses

    sigma'_ab = sigma^B_ab - (N kT / V) delta_ab      (instantaneous V)

as in Lips & Maass, PRE 97, 053002 (2018), Eqs. 12 and 15–16. The ideal-gas and `+K` terms are unchanged. The covariance subtracts the production-window mean of `sigma'`.

**Not affected:** C44, C55, C66, every NVT result, and K from the volume fluctuations.

**Verified:** on a short 500-atom Zn06 300 K NPT run, the code's C11 and C12 match an independent evaluation of the corrected formula from the same `series.h5` to 1e-4 GPa. Relative to the old formula they shift by −0.63 GPa.

## 2. `eam.alloy` reader: pair potentials for three or more elements

**File:** `src/python/get_potential.py` (`get_EAM_alloy_file`)

**Problem.** A setfl file lists the pair blocks as 11, 21, 22, 31, 32, 33, …, but the reader walked them as 11, 12, 13, 22, 23, 33, …. With three or more elements, pair potentials therefore went into the wrong slots. For three elements the old slot map was `[0,1,3,2,4,5]`. For two elements both orders are the same, so Cu–Zn was never affected.

**Fix.** The reader now walks `(i, j) = (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), …`, which matches `pot_index = j(j+1)/2 + i` in `MDDriver`.

**Verified against LAMMPS `pair_style eam/alloy`** on a distorted 500-atom Ni–Al–Co cell (Mishin 2013 potential):

| Quantity | Difference from LAMMPS |
|---|---|
| Energy | 5e-12 eV/atom |
| Forces | 2e-7 eV/Å |
| Virial | 1e-6 GPa |
| Born matrix | 4e-4 GPa |

The Zn–Cu control is unchanged.
