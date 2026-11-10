# mdelast

<p align="center"><img src="docs/logo/mdelast-banner-dark.png" alt="mdelast: elastic constants from stress fluctuations" width="640"></p>

`mdelast` is a C++17/Python molecular-dynamics program. It calculates
finite-temperature elastic constants of multi-component embedded-atom-method
(EAM) alloys from equilibrium stress fluctuations.

This is the CPC program package, version 1.2, archived on Zenodo under
[10.5281/zenodo.21380002](https://doi.org/10.5281/zenodo.21380002) (all versions).
It holds the source, an MIT license, documentation, a sample case (the 6 at.% Zn
alpha-brass calculation of the paper) with its expected output, and benchmarks.
Repository: <https://github.com/driven-physics/mdelast>.

## Files

- `LICENSE` — MIT license for mdelast.
- `CITATION.cff` — citation metadata.
- `README.md` — this file.
- `FIXES.md` — changes by version.
- `requirements.txt` — pinned Python packages.
- `docs/`
  - `docs/API_REFERENCE.md` — Python and compiled-extension API.
  - `docs/BUILDING.md` — clusters, macOS, build problems.
  - `docs/USAGE.md` — parameters, examples, output files.
  - `docs/logo/` — project logo.
- `examples/` — 500-atom example inputs.
  - `examples/ar_lj_gas/run.py` in `examples/ar_lj_gas/` — argon Lennard-Jones, NVT.
  - `examples/cu_300k/run.py` in `examples/cu_300k/` — pure copper, NPT, 300 K.
  - `examples/zn06_300k/run.py` in `examples/zn06_300k/` — 6 at.% Zn alloy, NPT, 300 K.
- `potentials/` — tabulated EAM potentials.
  - `potentials/Cu01.eam.alloy.txt` — Mishin et al., pure Cu.
  - `potentials/ZnCu.eam.alloy` — Clement–Auger, Cu-Zn.
  - `potentials/Mishin_updated-Ni-Al-Co-2013.eam.alloy` — Purja Pun–Yamakov–Mishin, Ni-Al-Co, used by benchmark S.
- `THIRD_PARTY_LICENSES/` — licenses and notices of the potentials.
  - `THIRD_PARTY_LICENSES/Clement-Auger-LGPL-2.1.txt` — LGPL-2.1 license from the OpenKIM distribution of the Cu-Zn potential.
  - `THIRD_PARTY_LICENSES/Clement-Auger-NOTICE.md` — attribution, provenance, modification notice, checksum.
  - `THIRD_PARTY_LICENSES/Mishin-Cu-CDDL-1.0.txt` — CDDL-1.0 license from the OpenKIM distribution of the pure-Cu potential.
  - `THIRD_PARTY_LICENSES/Mishin-Cu-NOTICE.md` — attribution, provenance, checksum.
  - `THIRD_PARTY_LICENSES/PurjaPun-NiAlCo-NOTICE.md` — attribution, provenance, checksum.
- `sample_case/` — the test of the paper: zero-pressure NPT, then NVT.
  - `sample_case/run.py` — run script.
  - `sample_case/validate.py` — compares a finished run with the reference tolerances.
  - `sample_case/expected_output/` — reference output of the manuscript run.
    - `sample_case/expected_output/result.json` — summary, replica scatter, tolerances.
    - `sample_case/expected_output/npt_series_decimated.csv` — NPT temperature, pressure, box length, every 500th step.
    - `sample_case/expected_output/nvt_series_decimated.csv` — NVT temperature, pressure, running means of the elastic constants, every 500th step.
- `benchmarks/` — benchmarks S, Z1, and L1 with expected values and tolerances.
  - `benchmarks/README.md` — what each benchmark checks, commands, runtimes.
  - `benchmarks/verify.py` — compares a run with `benchmarks/expected/` and prints PASS or FAIL for every quantity.
  - `benchmarks/bench_lib.py` — shared helpers.
  - `benchmarks/update_expected.py` — takes new expected values of Z1 or L1 from a reference run and keeps the tolerances.
  - `benchmarks/make_expected_md.py` — renders `benchmarks/expected/*.json` as `EXPECTED.md`.
  - `benchmarks/smoke/run_smoke.py` — benchmark S: static finite-difference checks, short NVE, NVT, and NPT runs, and a three-element Ni-Al-Co cell against stored LAMMPS values.
  - `benchmarks/smoke/reference/` — the Ni-Al-Co cell, the LAMMPS input, and the LAMMPS forces and values.
  - `benchmarks/l12_cuzn3_5K/run.py` — benchmark L1: L1₂ CuZn₃ at 5 K, NPT then NVT.
  - `benchmarks/expected/` — expected values and tolerances (`smoke.json`, `zn06_300K.json`, `l12_cuzn3_5K.json`, and the readable `EXPECTED.md`).
  - `benchmarks/slurm_example.sh` — one-core Slurm array that runs and verifies S, Z1, and L1.
- `tests/` — unit tests and short compiled runs.
  - `tests/test_averages.py` — averaging helpers.
  - `tests/test_eam_reader.py` — SETFL reader with three and four elements
    (synthetic file, no external potential).
  - `tests/test_examples_smoke.py` — short runs of all examples and of the sample case.
  - `tests/test_reference_data.py` — reference CSV, metadata, potential checksums.
  - `tests/test_tables.py` — table interpolation option and the quintic tables.
  - `tests/test_validator.py` — validator pass, fail, and tolerance boundaries.
- `src/` — program source.
  - `src/CMakeLists.txt` — CMake build of the Python extension.
  - `src/barostat.cpp`, `src/barostat.hpp` — isotropic barostat.
  - `src/barostat_flex.cpp`, `src/barostat_flex.hpp` — placeholder data structure for a
    future flexible-cell barostat, not used.
  - `src/eam.hpp` — EAM density and embedding functions with their derivatives.
  - `src/md_driver.cpp`, `src/md_driver.hpp` — molecular-dynamics driver.
  - `src/md_driver_module.cpp` — pybind11 bindings.
  - `src/numerics.hpp` — numerical helpers.
  - `src/poly_table.hpp` — interpolation tables of the potential functions (quintic and B-spline).
  - `src/potential.hpp` — pair-potential classes.
  - `src/thermostat.cpp`, `src/thermostat.hpp` — Nose-Hoover-chain thermostat.
  - `src/units.hpp` — CGS unit constants.
  - `src/vector3d.hpp` — three-dimensional vector class.
  - `src/python/` — Python layer.
    - `src/python/averages.py` — averaging utilities.
    - `src/python/get_potential.py` — potential-file readers and material helpers.
    - `src/python/molecular_dynamics.py` — setup, run loop, analysis, output.
    - `src/python/plot.py` — diagnostic plots.
    - `src/python/series.py` — time series and HDF5 output.
    - `src/python/vacf.py` — velocity autocorrelation.
- `tools/` — developer checks of the potential tables.
  - `tools/test_poly_table.cpp` — C++ unit test of the quintic table class (build command in its header).
  - `tools/test_table_shape.py` — roots, extrema, inflections, continuity, and accuracy of the tables of a potential file, quintic against B-spline, with plots.
  - `tools/eval_tables.py` — helper of `test_table_shape.py`.

`build/` and `sample_case_output/` are generated and not part of the package.

## Source provenance

Version 1.0 of `src/` and `examples/`, the potential data, `requirements.txt`,
and `LICENSE` are from the deposited `mdelast-main.zip`
(MD5 `928812de6fe872e058e2fade2c41ded5`) at
<https://doi.org/10.5281/zenodo.21380003>. `FIXES.md` lists the later changes.
Version 1.1 changes one line of `src/python/molecular_dynamics.py`. Version 1.2
adds the quintic tables (default, the B-spline tables remain available) and the
`sample_stride` option, and removes `-ffloat-store` from the Release flags in
`src/CMakeLists.txt`.

## Requirements

The reference environment was Linux x86-64 with:

- GCC 14.2.0 with C++17 support
- CMake 3.13 or newer for the commands below (CMake 4 needs `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`)
- Python 3.12 and Python development headers
- pybind11 2.13.6
- Eigen 3.4.0 (minimum 3.3)
- Boost 1.84.0 headers
- NumPy 1.26.2, SciPy 1.11.4, h5py 3.10.0, and Matplotlib 3.8.2

Keep NumPy 1.x: `src/python/vacf.py` calls `np.trapz`, which NumPy 2 no longer
provides. Other C++17 compilers and Python 3 versions may work. Results can
differ slightly with compiler, CPU, and floating-point optimization.
`docs/BUILDING.md` covers clusters with modules and macOS with Apple clang.

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libeigen3-dev libboost-dev \
  python3-dev python3-venv
```

Fedora/RHEL:

```bash
sudo dnf install gcc-c++ cmake eigen3-devel boost-devel python3-devel
```

Python environment, from the package root:

```bash
unset PYTHONPATH
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --isolated --upgrade pip
python3 -m pip install --isolated pybind11==2.13.6 -r requirements.txt
```

Build:

```bash
cmake -S src -B build/Release \
  -DCMAKE_BUILD_TYPE=Release \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
cmake --build build/Release --parallel
```

Import check:

```bash
PYTHONPATH="build/Release:src/python" python3 -c \
  "import md_driver_module, molecular_dynamics; print(md_driver_module.__file__)"
```

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The suite has unit tests, checks of the reference data, and short runs of the
three examples and of both stages of the sample case. The short runs use
10 thermalization and 20 production steps. They test that the code runs, not
its results. For an extension in another directory:

```bash
MDELAST_BUILD_DIR=/path/to/build python3 -m unittest discover -s tests -v
```

## Benchmarks

`benchmarks/` holds three benchmarks for verifying an installation: S (static
checks and short runs, under 1 min), Z1 (the sample case below, about 19 min),
and L1 (L1₂ CuZn₃ at 5 K, about 17 min), each on one core.
`benchmarks/verify.py` checks a run against `benchmarks/expected/` and prints
PASS or FAIL for every quantity. `benchmarks/README.md` gives the commands.

## Sample case

The sample case is the workflow of the paper for a 500-atom random FCC solid
solution with 30 Zn atoms (6 at.% Zn) at 300 K. It uses the Clement–Auger
potential, a 2 fs timestep, a 10-link Nose-Hoover chain with a 50 fs relaxation
time, and a neighbor list rebuilt every 20 steps. The stages are:

1. 50,000 steps of zero-pressure NPT equilibration and 100,000 steps of NPT production
2. transfer of the NPT mean volume and final configuration to NVT
3. 50,000 steps of NVT equilibration and 100,000 steps of NVT production, in
   which the elastic constants are accumulated

Run it from the package root:

```bash
python3 sample_case/run.py \
  --build-dir build/Release \
  --output-dir sample_case_output
```

The output directory must not exist. The run writes `npt/` and `nvt/` with the
full `series.h5` files and trajectory dumps, and at the top level:

- `result.json` — NPT and NVT summary values
- `npt_series_decimated.csv` — compact NPT series
- `nvt_series_decimated.csv` — compact NVT series

The full output needs about 1.3 GB of disk space. The compact reference files
in this package keep every 500th sample, without averaging.

A quick run that only checks the build:

```bash
python3 sample_case/run.py \
  --build-dir build/Release \
  --output-dir sample_case_quick \
  --cycles-thermalization 10 \
  --production-cycles 20
```

## Reference values and tolerances

Compare `sample_case_output/result.json` of the full run (200 ps of production
per stage) with `sample_case/expected_output/result.json`. The manuscript
values, from version 1.0, are:

- NPT mean box length: 18.33463 A
- NVT C11: 153.0657 GPa
- NVT C12: 102.6290 GPa
- NVT C44: 71.8933 GPa
- NVT kinetic temperature: 299.9901 K
- NVT configurational temperature: 300.5768 K as the mean of the per-step
  samples (299.9750 K for the reciprocal-of-mean estimator, which the program
  also prints)

The nine placement-seed replicas of the manuscript gave standard deviations of
0.3643, 0.3648, and 0.1048 GPa for C11, C12, and C44, 0.1022 K for the kinetic
temperature, 0.5900 K for the configurational temperature, and 0.00084 A for
the box length. The tolerances in the reference JSON are about three times the
replica scatter for the elastic constants and the configurational temperature:
1.1 GPa for C11/C12, 0.32 GPa for C44, and 1.8 K. The kinetic temperature has
the same 1.8 K, and the box length a wider 0.01 A.

Check a full run:

```bash
python3 sample_case/validate.py sample_case_output/result.json
```

The validator checks the case metadata, cycle counts, all three random seeds,
NPT and NVT box lengths, NVT elastic constants, and both temperature estimators
in both stages. It prints every comparison and exits with status 0 only when
all pass. Pressure is printed but not checked, because the NPT target is zero.

The same run is benchmark Z1, checked against
`benchmarks/expected/zn06_300K.json`:

```bash
python3 benchmarks/verify.py zn06_300K sample_case_output
```

## Citation

R. Shikler, M. Krief, and Y. Ashkenazy, “Elastic constants of multi-component
EAM alloys from single-trajectory stress fluctuations,” *Computer Physics
Communications* (submitted).

Until the article has its final bibliographic details, also cite the archived
software, <https://doi.org/10.5281/zenodo.21380002> (all versions).
`CITATION.cff` gives the citation metadata.

## Runtime

On one core of an Intel Xeon Gold 6534, version 1.2 with its default options
needs about 3.8 ms per step for the 500-atom 6 at.% Zn case, so the full
NPT-plus-NVT workflow (300,000 steps) takes about 19 minutes. With
`--interpolation bspline` (the tables of versions up to 1.1) it takes about
38 minutes (7.5 ms per step).

## Potentials

`potentials/ZnCu.eam.alloy` is the Clement–Auger file of the NIST Interatomic
Potentials Repository, where it is named `CuZn.eam.alloy`:

<https://www.ctcms.nist.gov/potentials/entry/2022--Clement-A-Auger-T--Cu-Zn/2022--Clement-A--Cu-Zn--LAMMPS--ipr1.html>

NIST states that Antoine Clement provided the file on 3 April 2023. OpenKIM
distributes the model as
`EAM_Dynamo_ClementAuger_2022_CuZn__MO_875501152989_000` under LGPL-2.1. The
potential is covered by the license and modification notice in
`THIRD_PARTY_LICENSES/`, not by mdelast's MIT license. Cite:

A. Clement and T. Auger, “An EAM potential for alpha-brass copper-zinc alloys:
application to plasticity and fracture,” *Modelling Simul. Mater. Sci. Eng.*
31 (2023) 015004, <https://doi.org/10.1088/1361-651X/aca4ec>.

`potentials/Cu01.eam.alloy.txt` is the pure-Cu potential of Mishin et al.,
distributed by NIST and OpenKIM. The file is byte-for-byte identical to that of
OpenKIM model `MO_346334655118_006`, which is under CDDL-1.0. Its license and
notice are in `THIRD_PARTY_LICENSES/`. mdelast's MIT license does not apply to
it. NIST asks users to credit the repository as well as the paper:

<https://www.ctcms.nist.gov/potentials/entry/2001--Mishin-Y-Mehl-M-J-Papaconstantopoulos-D-A-et-al--Cu-1/>

Y. Mishin et al., “Structural stability and lattice defects in copper: Ab
initio, tight-binding, and embedded-atom calculations,” *Phys. Rev. B* 63
(2001) 224106, <https://doi.org/10.1103/PhysRevB.63.224106>.

`potentials/Mishin_updated-Ni-Al-Co-2013.eam.alloy` is the Ni-Al-Co potential of
Purja Pun, Yamakov, and Mishin, used by Part C of benchmark S. The file is
byte-for-byte identical to the parameter file of OpenKIM model
`MO_826591359508_001`, which is dedicated to the public domain. Its notice is
in `THIRD_PARTY_LICENSES/`. NIST distributes the same data:

<https://www.ctcms.nist.gov/potentials/entry/2015--Purja-Pun-G-P-Yamakov-V-Mishin-Y--Ni-Al-Co/>

G. P. Purja Pun, V. Yamakov, and Y. Mishin, “Interatomic potential for the
ternary Ni-Al-Co system and application to atomistic modeling of the B2-L1_0
martensitic transformation,” *Modelling Simul. Mater. Sci. Eng.* 23 (2015)
065006, <https://doi.org/10.1088/0965-0393/23/6/065006>.

## Units and names

The C++ driver works in CGS units. Use the constants of
`md_driver_module.Units` when you build inputs. The public names `ensamble` and
`cycles_neigbors` keep their historical spelling.
