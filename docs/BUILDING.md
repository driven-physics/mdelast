# Building

`README.md` lists the requirements and the standard build commands. This file
covers clusters, macOS, the `make` builds, and common problems.

CMake 4 needs `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` at configure time, because
`src/CMakeLists.txt` declares version 3.0.

## Cluster with modules

```bash
module load gcc/14.2.0
module load python/3.12.1
module load pybind11/2.13.6
module load eigen/3.4.0
module load boost/1.84.0
module load cmake/3.31.7
```

Load CMake as a module too, since compute nodes may have no system CMake. The
pybind11 module sets `pybind11_DIR`. Pass it to CMake with
`-Dpybind11_DIR="$pybind11_DIR"`. On the reference cluster the `python/3.12.1`
module provides the packages of `requirements.txt` through its `PYTHONPATH`.
Clear `PYTHONPATH` before loading modules, not after, and set
`PYTHONNOUSERSITE=1` so that packages in `~/.local` do not interfere. If the
Python module lacks the packages, install them in a virtual environment made
with the module Python. Module names vary by system. CMake must find the
package configuration files of pybind11 and Eigen3.

## macOS with Apple clang

The Release flags contain two GCC-only options (`-Wno-maybe-uninitialized`
and `-rdynamic` at compile time), and `-Werror` turns clang's warnings about
them into errors. Boost is added as a regular include directory, so
`-Wdouble-promotion` warnings inside the Boost 1.84 headers also stop the
build. Without changing `src/`, configure with:

```bash
cmake -S src -B build/Release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_CXX_FLAGS="-isystem <prefix>/include -Wno-unused-command-line-argument -Wno-unknown-warning-option" \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
```

`<prefix>/include` is the directory with the Boost headers, for example the
`include/` directory of a conda environment or `/opt/homebrew/include`.
A clang build on arm64 is not bit-identical to the GCC x86-64 reference.

## Release build

From the top-level directory:

```bash
cmake -S src -B build/Release -DCMAKE_BUILD_TYPE=Release
make -j -C build/Release
```

The extension has a name like:

```text
build/Release/md_driver_module.cpython-312-x86_64-linux-gnu.so
```

Python sets the suffix. The Release build uses `-O2` with strict warning and
aliasing checks. Version 1.2 removed `-ffloat-store` from these flags
(`FIXES.md`, entry 6).

## Debug build

```bash
cmake -S src -B build/Debug -DCMAKE_BUILD_TYPE=Debug
make -j -C build/Debug
```

Pass `--build-dir build/Debug` to an example to use it. Compare results only
between builds with the same compiler and build type. Floating-point results
can differ across compilers, CPUs, and optimization levels.

## Quick run

The import check and the test suite are in `README.md`. The Ar LJ gas example
is the fastest full run:

```bash
python3 examples/ar_lj_gas/run.py --build-dir build/Release
```

## Common problems

- `Could not find pybind11`: load or install pybind11 and set `pybind11_DIR` to
  the directory with `pybind11Config.cmake`.
- `Could not find Eigen3`: load or install Eigen3 and set `Eigen3_DIR`.
- `ModuleNotFoundError: md_driver_module`: build the extension and pass the
  build directory, not `src/`.
- The module imports but fails at once: the Python interpreter at run time
  must match the interpreter and ABI that CMake used.
- Potential-file errors: keep `potentials/` next to `src/`. The loaders read
  that directory.
