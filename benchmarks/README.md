# mdelast benchmarks: inputs, expected outputs and verification

This directory documents how to verify an mdelast installation independently and how to reproduce the reference calculations of the paper
(R. Shikler, M. Krief, Y. Ashkenazy, *Elastic constants of multi-component EAM alloys from single-trajectory stress fluctuations*, Comput. Phys. Commun.).

All paths below are relative to the repository root.
It needs only the files already in the repository (`src/`, `potentials/ZnCu.eam.alloy`, `sample_case/run.py`).

| ID | Benchmark | What it verifies | Input | Runtime, 1 core |
|---|---|---|---|---|
| **S** | Smoke test: static finite-difference checks plus short NVE/NVT/NPT runs | Forces = −∇U; Born stress = (1/V) ∂U/∂η; Born constants C^B_ab = (1/V0) ∂²U/(∂η_a∂η_b) at zero strain; configurational-temperature Laplacian = Σ∇²U; net force = 0; absolute energy, pressure, Born constants, Laplacian and T_conf of a fixed cell; NVE energy error scales as dt²; extended-Hamiltonian conservation and temperatures in NVT and NPT | `benchmarks/smoke/run_smoke.py` | under 1 min |
| **Z1** | 6 at.% Zn α-brass, 300 K, NPT → NVT (the workflow of Sec. 5.1, Figs. 1–5) | Full production protocol: zero-pressure volume, C11, C12, C44, their Born, stress-fluctuation and kinetic parts, T_K, both T_conf estimators, pressure, conserved-energy error and drift | `sample_case/run.py` (already in the repository) | about 19 min |
| **L1** | Ordered L1₂ CuZn₃, 5 K, NPT → NVT (Table 1 of the manuscript) | Same quantities as Z1 near the athermal limit, plus agreement with the independent 0 K reference values of the potential (NIST repository) | `benchmarks/l12_cuzn3_5K/run.py` | about 17 min |

Runtimes were measured with version 1.2 and its default options on one core of an Intel Xeon Gold 6534 (about 3.8 ms per MD step for 500 atoms; see "Runtime" below).
Z1 and L1 each run 2 × 150,000 MD steps.

`benchmarks/verify.py` compares any of these outputs with the expected values in `benchmarks/expected/*.json`, prints `PASS`/`FAIL` for every quantity, and exits with status 0 only if nothing fails.
It needs only NumPy and h5py, so outputs can be checked on another machine.
`benchmarks/expected/EXPECTED.md` lists all expected values and tolerances in human-readable form.

## 1. Build

Build the compiled extension as described in the repository `README.md` (same requirements and commands), e.g.

```bash
cmake -S src -B build/Release -DCMAKE_BUILD_TYPE=Release \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
cmake --build build/Release --parallel
```

All scripts take `--build-dir` (default `build/Release`) and write only into the `--output-dir` you give them; the output directory must not exist yet.

## 2. Run and verify

```bash
# S: smoke test (<1 min)
python3 benchmarks/smoke/run_smoke.py --build-dir build/Release --output-dir benchmark_output/smoke
python3 benchmarks/verify.py smoke benchmark_output/smoke

# Z1: 6 at.% Zn, 300 K (~19 min); this is the repository's sample case
python3 sample_case/run.py --build-dir build/Release --output-dir benchmark_output/zn06_300K
python3 benchmarks/verify.py zn06_300K benchmark_output/zn06_300K

# L1: L1_2 CuZn3, 5 K (~17 min)
python3 benchmarks/l12_cuzn3_5K/run.py --build-dir build/Release --output-dir benchmark_output/l12_cuzn3_5K
python3 benchmarks/verify.py l12_cuzn3_5K benchmark_output/l12_cuzn3_5K


# everything that exists under benchmark_output/
python3 benchmarks/verify.py all benchmark_output
```

Each run is serial; Z1 and L1 are independent and can run at the same time on two cores.

All three run scripts take `--interpolation quintic` (default) or `--interpolation bspline`, the table interpolation of the potential (`docs/API_REFERENCE.md`, "Table interpolation"); `bspline` reproduces version 1.1 with the same build and inputs.
The choice is stored in `result.json`.
The fixed-cell values of S depend on it, and `verify.py` uses the matching expected values (an output without the field is taken as `bspline`).
Z1 and L1 are statistical checks and apply to both choices.
`benchmarks/slurm_example.sh` shows a one-core Slurm array that runs S, Z1 and L1 and verifies them.

A typical line of `verify.py` output is

```
  PASS  nvt.C11 = 153.24 GPa, expected 153.239 +/- 1 (diff +0.000, 0.00 of tol)
```

`INFO` lines are reported but not checked (see "Ensemble notes").
`verify.py` also fails if the run length differs from the reference schedule, because the tolerances hold only for the full 50,000 + 100,000-step runs.

## 3. Expected values (excerpt)

Cubic averages, GPa; the tolerance follows ±. Full tables, including the NPT stages, are in `expected/EXPECTED.md`.

| | C11 | C12 | C44 | C^B_11 | C^B_12 | C^B_44 | T_K [K] | T_conf [K] |
|---|---|---|---|---|---|---|---|---|
| Z1, NVT stage | 153.24 ± 1.0 | 102.63 ± 0.88 | 71.95 ± 0.42 | 163.396 ± 0.069 | 108.182 ± 0.050 | 76.071 ± 0.033 | 300.06 ± 0.83 | 299.48 ± 3.1 |
| L1, NVT stage | 97.625 ± 0.015 | 60.818 ± 0.023 | 42.840 ± 0.007 | 97.724 ± 0.008 | 60.829 ± 0.005 | 42.872 ± 0.005 | 4.999 ± 0.011 | 5.083 ± 0.2 |

Z1 also checks the NPT box ratio L/L0 = 1.008903 ± 0.000017 (box length 18.3343 Å), and L1 checks C(5 K) against the 0 K values 97.638, 60.838, 42.885 GPa within 0.2%.
T_conf in this table is the reciprocal-of-mean estimator 1/⟨1/T_conf⟩; the per-step mean ⟨T_conf⟩ is checked too.

Benchmark S checks finite-difference consistency at the levels listed in `expected/EXPECTED.md` (forces to 2×10⁻⁵ eV/Å, Born stress to 10⁻⁵ GPa, Born constants to 0.05 GPa, Laplacian to 2×10⁻⁵ relative) and reproduces the energy, pressure, Born constants, Laplacian and T_conf of its fixed cell to a relative 10⁻⁷.
With the default quintic tables the Born-constant residual is about 1.0×10⁻⁴ GPa and the Laplacian residual 3×10⁻⁸ (relative), because the tabulated first and second derivatives are the exact derivatives of the interpolated function.
With `--interpolation bspline` the Born-constant residual is about 0.013 GPa (10⁻⁴ relative) and does not change with the finite-difference step between 10⁻⁴ and 10⁻⁶; it reflects the separately interpolated tables of the potential and its first and second derivatives, not truncation error.

## 4. How the tolerances were set

For every checked quantity q of Z1 and L1:

* **σ_run(q) = max(SE_run, s_vel)**.
  * SE_run is the standard error of the reference run: block bootstrap of the full nonlinear estimator for C_ij and the fluctuation terms (64 blocks of 3.1 ps, 2000 resamples), and the blocking standard error for per-step averages (Born terms, temperatures, pressure, box ratio).
  * s_vel is the standard deviation of q over nine runs that differ only in the initial velocities (same atomic arrangement). It includes the propagation of the NPT volume into the NVT stage, which dominates for the Born terms.
* **tol(q) = 4 × √2 × σ_run(q)**. The factor √2 holds because a re-run is an independent sample compared with one reference sample. The factor 4 keeps the false-alarm probability of the ~30 checks of one benchmark below 0.2% for Gaussian errors.
* **Conserved energy**: one-sided bounds of three times the largest value seen in the reference run and its velocity replicas.

A run on another compiler or CPU follows a different chaotic trajectory, so it is checked statistically.
The expected values of Z1 and L1 come from version 1.2 reference runs with the default options and the original seeds (`benchmarks/update_expected.py`); SE_run, s_vel and the tolerances are those of the version 1.0 reference set, whose values are kept as `expected_previous` in the JSON files.
Runs of versions 1.1 and 1.2, with either table interpolation and with `sample_stride` 1 or 10, pass every check.
Version 1.0 outputs fail two Z1 checks, the NPT box ratio L/L0 and the kinetic term (by 1.1 tolerances), because version 1.0 counted 3N degrees of freedom and ran 0.2% hotter (`FIXES.md`, entry 3).

The tolerances are deliberately tight enough to detect a different random Zn arrangement. Running Z1 with another placement seed fails the Born-term and volume checks, because the configurational scatter of C^B_11 between arrangements (0.28 GPa) exceeds the tolerance (0.07 GPa).
Keep the seeds in the scripts: placement seed 42 and velocity seeds 0/0 for Z1; velocity seeds 0/500 for L1. The arrangement comes from NumPy's legacy `RandomState`, whose stream NumPy keeps fixed across versions.

The fixed-configuration values of S are deterministic and use relative tolerances of 10⁻⁷, with one set of values per table interpolation. Its finite-difference bounds are about five times the residuals observed with the B-spline tables, which the quintic tables meet with a wide margin.

## 5. Relation to the paper

| Paper item | Benchmark | Compare |
|---|---|---|
| Sec. 5.1, Figs. 1–5 (workflow, T_K/T_conf, NPT volume, Born convergence, conserved energy) | Z1 | `npt/series.h5`, `nvt/series.h5`: same seeds as the figures (which were made with version 1.0) |
| Sec. 5.3 / Fig. 11 (uncertainty) | Z1 | the per-run standard errors and replica scatter that set the tolerances (`expected/zn06_300K.json`: `se_run`, `s_vel`) |
| Table 1, CuZn₃ L1₂ | L1 | `nvt.C11/C12/C44`; the table quotes the mean of nine velocity replicas (97.634, 60.824, 42.841 GPa), L1 is replica 0 |
| Sec. 5.5 (T_K vs T_conf) and Appendix (Laplacian) | S, Z1, L1 | S checks the Laplacian against finite differences; Z1 and L1 check both T_conf estimators |
| Sec. 5.6 (runtime) | all | wall time per step printed by S for each short run (`dynamics.*.ms_per_step`) |

For Z1 and L1, `verify.py` reads `npt/series.h5` and `nvt/series.h5`. Their tolerances apply to the reference state points, seeds, and 150,000-step schedule in each stage.
The quantity extraction is in `benchmarks/bench_lib.py`, if you want to tabulate other state points.

## 6. Ensemble notes and known limitations

* **NPT C11, C12 and K** are printed as `INFO` and are not checked.
  * The paper recommends the NVT stage for the elastic constants (Sec. 5.3), and the NPT stress-fluctuation estimates of C11 and C12 converge slowly (per-run standard error about 11 GPa at 300 K).
  * NPT C44, the Born terms and the volume are checked.
  * This version of the code includes the corrected NPT estimator (the ideal-gas pressure N k_B T / V stays inside the fluctuating stress; `FIXES.md`, item 1). With the corrected formula the version 1.0 reference trajectories gave Z1 146.765 and 96.278 GPa (0.689 GPa below the previous estimator) and L1 46.051 and 9.246 GPa (0.014 GPa below). The INFO values now printed are those of the version 1.2 reference runs (`expected/EXPECTED.md`). NPT C44, K and every NVT quantity are unchanged by the fix.
* **Shear components.** mdelast uses an orthorhombic cell, so benchmark S cannot apply shear strains. Shear components are checked through C44 and C^B_44 by Z1 and L1.
* **NVE runs.** In the current code, `MolecularDynamics.run()` stops at the final averaging step for the NVE ensemble. Benchmark S therefore drives the compiled integrator directly for its NVE check.
* **Disk space.** Z1 and L1 each write about 1.3 GB (full per-step `series.h5` and trajectory dumps).

## 7. Runtime

On one core of an Intel Xeon Gold 6534 (version 1.2, GCC 14.2, Release build), for both table interpolations:

| Benchmark | MD steps | Wall time, quintic (default) | Wall time, bspline |
|---|---|---|---|
| S | 7,000 MD steps + ~1,600 static evaluations | 35 s | 71 s |
| Z1 | 300,000 | 1,141 s (3.8 ms/step) | 2,256 s (7.5 ms/step) |
| L1 | 300,000 | 1,032 s (3.4 ms/step) | 2,001 s (6.7 ms/step) |

With `sample_stride = 10` (observables every tenth step) the quintic Z1 and L1 runs took 729 s and 676 s.
These times are for the default Release build of `src/CMakeLists.txt` (`-O2`), with other single-core jobs running on the same node.
Versions up to 1.1 also compiled with `-ffloat-store`, which made a step 2.2 times slower without changing any result (`FIXES.md`, entry 6).
Compiling with `-O3` instead of `-O2` shortens Z1 by about 6% (1,073 s) and gives bit-identical results; the reference times and values are those of the default build.

## 8. Files

```
benchmarks/
  README.md                     this file
  verify.py                     PASS/FAIL comparison for all benchmarks
  update_expected.py            new Z1/L1 expected values from a reference run (tolerances kept)
  make_expected_md.py           renders expected/*.json as expected/EXPECTED.md
  bench_lib.py                  reading series.h5, comparison helpers
  slurm_example.sh              one-core Slurm array: S, Z1, L1 + verification
  expected/
    smoke.json zn06_300K.json l12_cuzn3_5K.json   expected values and tolerances
    EXPECTED.md                 the same, as tables
  smoke/run_smoke.py            benchmark S
  l12_cuzn3_5K/run.py           benchmark L1
```

Benchmark Z1 uses the repository's `sample_case/run.py` unchanged; its compact reference output is in `sample_case/expected_output/`.
