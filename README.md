# mdelast

<p align="center"><img src="docs/logo/mdelast-banner-dark.png" alt="mdelast: elastic constants from stress fluctuations" width="640"></p>

`mdelast` is a C++17/Python molecular-dynamics program for calculating
finite-temperature elastic constants of multi-component embedded-atom-method
(EAM) alloys from equilibrium stress fluctuations.

This CPC program package, version 1.2, is archived on Zenodo under
[10.5281/zenodo.21380002](https://doi.org/10.5281/zenodo.21380002) (all versions).
It contains the program source, an MIT
license, user documentation, a complete input/output test based on the
manuscript's representative 6 at.% Zn alpha-brass calculation, and
verification benchmarks. Repository: <https://github.com/driven-physics/mdelast>.

## Package inventory

Every distributed file and directory is listed below.

- `LICENSE` — MIT license for mdelast.
- `CITATION.cff` — citation metadata for this version.
- `README.md` — package inventory, installation, execution, and validation instructions.
- `requirements.txt` — pinned Python runtime dependencies.
- `docs/` — user and developer documentation.
  - `docs/API_REFERENCE.md` — Python and compiled-extension API reference.
  - `docs/BUILDING.md` — detailed build requirements and troubleshooting.
  - `docs/USAGE.md` — simulation parameters, examples, and output format.
  - `docs/logo/` — project logo.
- `examples/` — standalone example inputs.
  - `examples/ar_lj_gas/` — 500-atom argon Lennard-Jones NVT example.
    - `examples/ar_lj_gas/run.py` — argon example run script.
  - `examples/cu_300k/` — 500-atom pure-copper NPT example at 300 K.
    - `examples/cu_300k/run.py` — pure-copper example run script.
  - `examples/zn06_300k/` — 500-atom, 6 at.% Zn alloy NPT example at 300 K.
    - `examples/zn06_300k/run.py` — alloy example run script.
- `potentials/` — tabulated EAM potentials used by the examples and sample case.
  - `potentials/Cu01.eam.alloy.txt` — Mishin et al. pure-Cu EAM potential.
  - `potentials/ZnCu.eam.alloy` — Clement–Auger Cu-Zn EAM potential.
- `THIRD_PARTY_LICENSES/` — licenses and notices for redistributed third-party data.
  - `THIRD_PARTY_LICENSES/Clement-Auger-LGPL-2.1.txt` — LGPL-2.1 license supplied with the OpenKIM distribution of the Cu-Zn potential.
  - `THIRD_PARTY_LICENSES/Clement-Auger-NOTICE.md` — potential attribution, provenance, modification notice, and checksum.
  - `THIRD_PARTY_LICENSES/Mishin-Cu-CDDL-1.0.txt` — CDDL-1.0 license supplied with the OpenKIM distribution of the pure-Cu potential.
  - `THIRD_PARTY_LICENSES/Mishin-Cu-NOTICE.md` — pure-Cu potential attribution, provenance, and checksum.
- `sample_case/` — comprehensive manuscript test: zero-pressure NPT followed by NVT elasticity accumulation.
  - `sample_case/run.py` — self-contained NPT-to-NVT input/run script.
  - `sample_case/validate.py` — automated comparison of a completed run against the reference tolerances.
  - `sample_case/expected_output/` — reference output from the manuscript run.
    - `sample_case/expected_output/result.json` — reference summary, replica scatter, and acceptance tolerances.
    - `sample_case/expected_output/npt_series_decimated.csv` — NPT temperature, pressure, and box-length series, every 500th step.
    - `sample_case/expected_output/nvt_series_decimated.csv` — NVT temperature, pressure, and elastic-constant running means, every 500th step.
- `benchmarks/` — verification benchmarks S, Z1, and L1 with expected values and tolerances.
  - `benchmarks/README.md` — what each benchmark checks, how to run and verify it, and runtimes.
  - `benchmarks/verify.py` — compares a run with `benchmarks/expected/` and prints PASS or FAIL for every quantity.
  - `benchmarks/bench_lib.py` — helpers shared by the benchmark scripts.
  - `benchmarks/update_expected.py` — replaces the expected values of Z1 or L1 by those of a new reference run, keeping the tolerances.
  - `benchmarks/make_expected_md.py` — renders `benchmarks/expected/*.json` as `EXPECTED.md`.
  - `benchmarks/smoke/run_smoke.py` — benchmark S: static finite-difference checks and short NVE, NVT, and NPT runs.
  - `benchmarks/l12_cuzn3_5K/run.py` — benchmark L1: L1₂ CuZn₃ at 5 K, NPT followed by NVT.
  - `benchmarks/expected/` — expected values and tolerances (`smoke.json`, `zn06_300K.json`, `l12_cuzn3_5K.json`, and the readable `EXPECTED.md`).
  - `benchmarks/slurm_example.sh` — one-core Slurm array that runs and verifies S, Z1, and L1.
- `tests/` — standard-library unit, integrity, and compiled smoke tests.
  - `tests/test_averages.py` — statistical-helper unit tests.
  - `tests/test_examples_smoke.py` — short compiled runs of all examples and the NPT-to-NVT sample.
  - `tests/test_reference_data.py` — reference CSV, metadata, and potential-checksum integrity tests.
  - `tests/test_tables.py` — table interpolation option and properties of the quintic tables.
  - `tests/test_validator.py` — validator pass/fail and tolerance-boundary tests.
- `src/` — complete program source.
  - `src/CMakeLists.txt` — CMake build definition for the Python extension.
  - `src/barostat.cpp` — isotropic barostat implementation.
  - `src/barostat.hpp` — isotropic barostat declarations.
  - `src/barostat_flex.cpp`, `src/barostat_flex.hpp` — placeholder data structure for a
    future flexible-cell barostat, not used.
  - `src/eam.hpp` — EAM density and embedding functions with their derivatives.
  - `src/md_driver.cpp` — molecular-dynamics driver implementation.
  - `src/md_driver.hpp` — molecular-dynamics driver declarations.
  - `src/md_driver_module.cpp` — pybind11 Python bindings.
  - `src/numerics.hpp` — numerical helper routines.
  - `src/poly_table.hpp` — interpolation tables of the potential functions (quintic and B-spline).
  - `src/potential.hpp` — pair-potential classes.
  - `src/thermostat.cpp` — Nose-Hoover-chain thermostat implementation.
  - `src/thermostat.hpp` — thermostat declarations.
  - `src/units.hpp` — CGS unit constants.
  - `src/vector3d.hpp` — three-dimensional vector class.
  - `src/python/` — high-level Python simulation and analysis source.
    - `src/python/averages.py` — statistical averaging utilities.
    - `src/python/get_potential.py` — potential-file parsers and material helpers.
    - `src/python/molecular_dynamics.py` — simulation setup, run loop, analysis, and output.
    - `src/python/plot.py` — diagnostic plotting routines.
    - `src/python/series.py` — time-series containers and HDF5 serialization.
    - `src/python/vacf.py` — velocity-autocorrelation analysis.

- `tools/` — developer checks of the potential tables.
  - `tools/test_poly_table.cpp` — C++ unit test of the quintic table class (build command in its header).
  - `tools/test_table_shape.py` — roots, extrema, inflections, continuity, and accuracy of the tables of a potential file, quintic against B-spline, with plots.
  - `tools/eval_tables.py` — helper of `test_table_shape.py`.

Generated directories such as `build/` and `sample_case_output/` are not part of
the distributed package.

## Source provenance

Version 1.0 of the program source in `src/` and `examples/`, the potential data,
`requirements.txt`, and `LICENSE` are from the deposited `mdelast-main.zip`
(MD5 `928812de6fe872e058e2fade2c41ded5`) at
<https://doi.org/10.5281/zenodo.21380003>. The CPC-specific README,
`sample_case/`, `benchmarks/`, third-party notices, and documentation
corrections do not modify the deposited program source. Later changes to the
source are listed in `FIXES.md`. Version 1.1 differs
from the deposited source in one line of `src/python/molecular_dynamics.py`
(`FIXES.md`, version 1.1). Version 1.2 adds the quintic table interpolation of
the potential (default; the previous B-spline tables remain available) and the
`sample_stride` option, and removes `-ffloat-store` from the Release flags in
`src/CMakeLists.txt` (`FIXES.md`, version 1.2).

## Requirements

The reference environment was Linux x86-64 with:

- GCC 14.2.0 with C++17 support;
- CMake 3.13 or newer for the commands below (CMake 4 needs `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`);
- Python 3.12 and Python development headers;
- pybind11 2.13.6;
- Eigen 3.4.0 (minimum 3.3);
- Boost 1.84.0 headers;
- NumPy 1.26.2, SciPy 1.11.4, h5py 3.10.0, and Matplotlib 3.8.2.

Keep NumPy 1.x: `src/python/vacf.py` calls `np.trapz`, which current NumPy 2
releases no longer provide.

Other C++17 compilers and Python 3 versions may work, but regression values can
vary slightly with compiler, CPU, and floating-point optimization.
`docs/BUILDING.md` gives the module setup for clusters and the configure
options for macOS with Apple clang.

On Ubuntu/Debian, install the system build dependencies with:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libeigen3-dev libboost-dev \
  python3-dev python3-venv
```

On Fedora/RHEL-family systems, use:

```bash
sudo dnf install gcc-c++ cmake eigen3-devel boost-devel python3-devel
```

From the package root, create the Python environment:

```bash
unset PYTHONPATH
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --isolated --upgrade pip
python3 -m pip install --isolated pybind11==2.13.6 -r requirements.txt
```

Configure and compile:

```bash
cmake -S src -B build/Release \
  -DCMAKE_BUILD_TYPE=Release \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
cmake --build build/Release --parallel
```

Confirm that the extension and Python layer import:

```bash
PYTHONPATH="build/Release:src/python" python3 -c \
  "import md_driver_module, molecular_dynamics; print(md_driver_module.__file__)"
```

## Automated tests

After building the extension, run the complete test suite:

```bash
python3 -m unittest discover -s tests -v
```

This runs unit and reference-integrity tests plus short end-to-end simulations
of all three examples and both stages of the manuscript sample. The smoke runs
use 10 thermalization and 20 production steps and test execution only; they are
not scientific regressions. If the extension is in a non-default location:

```bash
MDELAST_BUILD_DIR=/path/to/build python3 -m unittest discover -s tests -v
```

## Benchmarks

`benchmarks/` holds three benchmarks for verifying an installation: S (static
checks and short runs, under 1 min), Z1 (the sample case below, about 19 min),
and L1 (L1₂ CuZn₃ at 5 K, about 17 min), each on one core.
`benchmarks/verify.py` checks a run against `benchmarks/expected/` and prints
PASS or FAIL for every quantity. `benchmarks/README.md` gives the commands.

## Comprehensive sample case

The sample reproduces the manuscript workflow for a 500-atom random FCC solid
solution containing 30 Zn atoms (6 at.% Zn) at 300 K. It uses the Clement–Auger
EAM potential, a 2 fs timestep, a 10-link Nose-Hoover thermostat chain with
50 fs relaxation time, neighbor-list rebuilding every 20 steps, and:

1. 50,000-step zero-pressure NPT equilibration plus 100,000-step NPT production;
2. transfer of the NPT mean volume and final configuration to NVT;
3. 50,000-step NVT equilibration plus 100,000-step NVT production while
   accumulating the elastic constants.

From the package root, after activating the environment and building, run:

```bash
python3 sample_case/run.py \
  --build-dir build/Release \
  --output-dir sample_case_output
```

The output directory must not already exist. It will contain `npt/` and `nvt/`
subdirectories with full `series.h5` files and trajectory dumps, plus these
compact files at its root:

- `result.json` — NPT and NVT summary values;
- `npt_series_decimated.csv` — compact NPT series;
- `nvt_series_decimated.csv` — compact NVT series.

The full generated output requires approximately 1.3 GB of free disk space.
The compact reference files included in this package are decimated by retaining
every 500th sample; no averaging was applied during decimation.

For a fast build/run check only (not a scientifically meaningful comparison):

```bash
python3 sample_case/run.py \
  --build-dir build/Release \
  --output-dir sample_case_quick \
  --cycles-thermalization 10 \
  --production-cycles 20
```

## Reference values and tolerances

For the complete 200 ps production phases, compare
`sample_case_output/result.json` with
`sample_case/expected_output/result.json`. The manuscript reference values,
produced with version 1.0, are:

- NPT mean box length: 18.33463 A;
- NVT C11: 153.0657 GPa;
- NVT C12: 102.6290 GPa;
- NVT C44: 71.8933 GPa;
- NVT kinetic temperature: 299.9901 K;
- NVT configurational temperature: 300.5768 K when per-step temperature
  samples are averaged (299.9750 K for the reciprocal-of-mean estimator also
  emitted by the program).

The nine placement-seed replicas in the manuscript gave standard deviations of
0.3643, 0.3648, and 0.1048 GPa for C11, C12, and C44; 0.1022 K for kinetic
temperature; 0.5900 K for configurational temperature; and 0.00084 A for box
length. The acceptance tolerances in the reference JSON are about three times the
replica scatter for the elastic constants and the configurational temperature:
1.1 GPa for C11/C12, 0.32 GPa for C44, and 1.8 K. The kinetic temperature has
the same 1.8 K, and the box length a wider 0.01 A. Small
platform-dependent changes are expected; a systematic excess beyond these
tolerances warrants investigation.

Validate a completed full run automatically:

```bash
python3 sample_case/validate.py sample_case_output/result.json
```

The validator checks the case metadata, cycle counts, all three random seeds,
NPT/NVT box lengths, NVT elastic constants, and both temperature estimators in
both stages. It prints every comparison and exits with status 0 only when all
checks pass. Pressure is reported for diagnosis but is not used as an
acceptance criterion because the NPT target is effectively zero.

The same run can also be checked as benchmark Z1, against
`benchmarks/expected/zn06_300K.json`:

```bash
python3 benchmarks/verify.py zn06_300K sample_case_output
```

## Please cite

When using mdelast, please cite the associated paper:

R. Shikler, M. Krief, and Y. Ashkenazy, “Elastic constants of multi-component
EAM alloys from single-trajectory stress fluctuations,” *Computer Physics
Communications* (submitted).

Until the article receives its final bibliographic details, also cite the
archived software, <https://doi.org/10.5281/zenodo.21380002> (all versions).
`CITATION.cff` gives the citation metadata.

## Runtime estimate

On one core of an Intel Xeon Gold 6534, version 1.2 with its default options
needs about 3.8 ms per step for the 500-atom 6 at.% Zn case, so the complete
NPT-plus-NVT workflow (300,000 steps) takes about 19 minutes. With
`--interpolation bspline` (the tables of versions up to 1.1) it takes about
38 minutes (7.5 ms per step). Actual time depends on hardware and compiler.

## Potential provenance

`potentials/ZnCu.eam.alloy` is the Clement–Auger file distributed by the NIST
Interatomic Potentials Repository (the repository file is named
`CuZn.eam.alloy`):

<https://www.ctcms.nist.gov/potentials/entry/2022--Clement-A-Auger-T--Cu-Zn/2022--Clement-A--Cu-Zn--LAMMPS--ipr1.html>

NIST states that the file was provided by Antoine Clement on 3 April 2023.
OpenKIM distributes the model as
`EAM_Dynamo_ClementAuger_2022_CuZn__MO_875501152989_000` under LGPL-2.1. The
potential is therefore covered by the separate license and modification notice
in `THIRD_PARTY_LICENSES/`, not by mdelast's root MIT license. It is included
because it is required by the test and was part of the cited Zenodo software
deposit. Cite:

A. Clement and T. Auger, “An EAM potential for alpha-brass copper-zinc alloys:
application to plasticity and fracture,” *Modelling Simul. Mater. Sci. Eng.*
31 (2023) 015004, <https://doi.org/10.1088/1361-651X/aca4ec>.

`potentials/Cu01.eam.alloy.txt` is the Mishin et al. pure-Cu potential
distributed by NIST and OpenKIM. The bundled file is byte-for-byte identical to
OpenKIM model `MO_346334655118_006`, which is distributed under CDDL-1.0. Its
separate license and attribution notice are included in
`THIRD_PARTY_LICENSES/`; mdelast's MIT license does not apply to it. NIST asks
users to credit the repository as well as the original paper:

<https://www.ctcms.nist.gov/potentials/entry/2001--Mishin-Y-Mehl-M-J-Papaconstantopoulos-D-A-et-al--Cu-1/>

Y. Mishin et al., “Structural stability and lattice defects in copper: Ab
initio, tight-binding, and embedded-atom calculations,” *Phys. Rev. B* 63
(2001) 224106, <https://doi.org/10.1103/PhysRevB.63.224106>.

## Units and compatibility

The C++ driver uses CGS internally. Use constants from
`md_driver_module.Units` when constructing inputs. The historical public names
`ensamble` and `cycles_neigbors` are intentionally retained for compatibility.
