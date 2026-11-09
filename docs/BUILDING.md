# Building

## Requirements

The extension requires:

- a C++17 compiler
- CMake 3.0 or newer (CMake 4 needs `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` at
  configure time, because `src/CMakeLists.txt` declares version 3.0)
- Python 3 and its development headers
- pybind11
- Eigen3 3.3 or newer
- Boost headers

The Release environment used GCC 14.2.0, pybind11 2.13.6, Eigen 3.4.0,
Boost 1.84.0, and Python 3.12. The Python runtime additionally requires the
packages in `requirements.txt`.

On a module-based cluster, an equivalent environment is:

```bash
module load gcc/14.2.0
module load python/3.12.1
module load pybind11/2.13.6
module load eigen/3.4.0
module load boost/1.84.0
module load cmake/3.31.7
```

Compute nodes may have no system CMake even when the login node has one, so
load CMake as a module too. The pybind11 module sets `pybind11_DIR`; pass it
to CMake with `-Dpybind11_DIR="$pybind11_DIR"`. On the reference cluster the
`python/3.12.1` module provides the packages of `requirements.txt` through
its `PYTHONPATH`; clear `PYTHONPATH` before loading modules, not after, and
set `PYTHONNOUSERSITE=1` so that packages in `~/.local` do not interfere.
Where the Python module lacks these packages, install them in a virtual
environment created with the module Python. Package names and module commands
vary by system. CMake must be able to locate the package configuration files
for pybind11 and Eigen3.

## macOS with Apple clang

The Release flags contain two GCC-only options (`-Wno-maybe-uninitialized`
and `-rdynamic` at compile time), and `-Werror` turns clang's warnings about
them into errors. Boost is added as a regular
include directory, so `-Wdouble-promotion` warnings inside the Boost 1.84
headers also stop the build. Without changing `src/`, configure with:

```bash
cmake -S src -B build/Release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_CXX_FLAGS="-isystem <prefix>/include -Wno-unused-command-line-argument -Wno-unknown-warning-option" \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
```

Here `<prefix>/include` is the directory that holds the Boost headers, for
example the `include/` directory of a conda environment or `/opt/homebrew/include`.
A clang build on arm64 is not bit-identical to the GCC x86-64 reference.

## Release build

From the top-level release directory:

```bash
cmake -S src -B build/Release -DCMAKE_BUILD_TYPE=Release
make -j -C build/Release
```

The resulting extension is named similar to:

```text
build/Release/md_driver_module.cpython-312-x86_64-linux-gnu.so
```

The exact suffix is selected by Python. The Release configuration uses `-O2`
with strict warning and aliasing checks. Version 1.2 removed `-ffloat-store`
from these flags (`FIXES.md`, entry 6).

## Debug build

```bash
cmake -S src -B build/Debug -DCMAKE_BUILD_TYPE=Debug
make -j -C build/Debug
```

Pass `--build-dir build/Debug` to an example to use this module. Regression
comparisons should use builds made with the same compiler and build type;
floating-point evaluation can differ across compilers, CPUs, and optimization
levels.

## Import check

```bash
PYTHONPATH="build/Release:src/python" python3 -c \
  "import md_driver_module, molecular_dynamics; print(md_driver_module.__file__)"
```

## Quick run check

After a successful import, the Ar LJ gas example is the fastest way to verify
the full simulation pipeline works end-to-end:

```bash
python3 examples/ar_lj_gas/run.py --build-dir build/Release
```

## Automated test suite

After building, run:

```bash
python3 -m unittest discover -s tests -v
```

The suite includes statistical-helper tests, reference-data and potential
checksum checks, validator tests, and short compiled smoke runs of every
example plus the two-stage sample case. Set `MDELAST_BUILD_DIR` when testing an
extension outside `build/Release`.

## Common problems

- `Could not find pybind11`: load/install pybind11 and set `pybind11_DIR` to
  the directory containing `pybind11Config.cmake`.
- `Could not find Eigen3`: load/install Eigen3 and set `Eigen3_DIR`.
- `ModuleNotFoundError: md_driver_module`: build the extension and pass the
  correct build directory; do not point at `src/`.
- A module imports but immediately fails: ensure the Python interpreter used
  at runtime matches the interpreter and ABI used by CMake.
- Potential-file errors: run from this release and keep `potentials/` adjacent
  to `src/`; the loaders resolve that package-relative directory.
