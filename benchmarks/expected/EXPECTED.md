# Expected benchmark values

Generated from `expected/*.json` (the files `verify.py` reads) by `make_expected_md.py`. All stresses and elastic constants in GPa, temperatures in K. Cubic averages: C11 = (C11+C22+C33)/3, C12 = (C12+C13+C23)/3, C44 = (C44+C55+C66)/3.

Tolerance rule for Z1 and L1: tol = 4 x sqrt(2) x max(per-run SE, velocity-replica std); conserved-energy checks are upper bounds of 3 x the largest value in the reference set. See README.md, Section 'How the tolerances were set'.

## Benchmark S: smoke test (`smoke/run_smoke.py`)

The fixed-cell values depend on the table interpolation of the potential; `verify.py` uses the column that matches the `interpolation` field of `result.json` (outputs of versions <= 1.1 have no such field and used the B-spline tables). The other bounds apply to both.

| check | unit | expected / bound, quintic tables | expected / bound, bspline tables | observed, quintic | observed, bspline |
|---|---|---|---|---|---|
| `static.net_force_rel` |  | <= 1e-12 | <= 1e-12 | 1.3e-15 | 2.93e-15 |
| `static.force_fd_max_abs_err_eV_A` | eV/A | <= 2e-05 | <= 2e-05 | 3.64e-06 | 3.68e-06 |
| `static.born_stress_fd_max_abs_err_GPa` | GPa | <= 1e-05 | <= 1e-05 | 1.55e-06 | 1.21e-06 |
| `static.born_C_fd_max_abs_err_GPa` | GPa | <= 0.05 | <= 0.05 | 0.000105 | 0.01288 |
| `static.laplacian_fd_rel_err` |  | <= 2e-05 | <= 2e-05 | 2.89e-08 | 3.69e-06 |
| `static.static_epot_eV_per_atom` | eV | -3.49917801 +/- 3.5e-07 | -3.49917801 +/- 3.5e-07 | = expected | = expected |
| `static.static_pressure_GPa` | GPa | 5.131070279 +/- 5.1e-07 | 5.131072142 +/- 5.1e-07 | = expected | = expected |
| `static.static_CB11_GPa` | GPa | 197.3963943 +/- 2e-05 | 197.396401 +/- 2e-05 | = expected | = expected |
| `static.static_CB12_GPa` | GPa | 126.5757292 +/- 1.3e-05 | 126.5757319 +/- 1.3e-05 | = expected | = expected |
| `static.static_CB44_GPa` | GPa | 92.13965308 +/- 9.2e-06 | 92.13965579 +/- 9.2e-06 | = expected | = expected |
| `static.laplacian_eV_A2` | eV/A^2 | 6567.770289 +/- 0.00066 | 6567.770385 +/- 0.00066 | = expected | = expected |
| `static.static_T_conf_K` | K | 1000.120485 +/- 0.0001 | 1000.120487 +/- 0.0001 | = expected | = expected |
| `ternary.epot_absdiff_eV_per_atom` | eV | <= 1e-10 | <= 1e-10 | 1.39e-12 | 4.72e-12 |
| `ternary.force_max_absdiff_eV_A` | eV/A | <= 2e-06 | <= 2e-06 | 9.25e-08 | 1.71e-07 |
| `ternary.virial_max_absdiff_GPa` | GPa | <= 1e-05 | <= 1e-05 | 1.09e-06 | 9.37e-07 |
| `ternary.born_max_absdiff_GPa` | GPa | <= 0.004 | <= 0.004 | 0.000424 | 0.000425 |
| `dynamics.nve_dt_ratio` |  | >= 3.5 and <= 4.5 | >= 3.5 and <= 4.5 | 3.998 | 3.998 |
| `dynamics.nve_dt2.econs_err_max` |  | <= 5e-05 | <= 5e-05 | 1.06e-05 | 1.06e-05 |
| `dynamics.nve_dt1.econs_err_max` |  | <= 1.5e-05 | <= 1.5e-05 | 2.65e-06 | 2.65e-06 |
| `dynamics.nvt.econs_err_max` |  | <= 5e-05 | <= 5e-05 | 1.36e-05 | 1.37e-05 |
| `dynamics.npt.econs_err_max` |  | <= 5e-05 | <= 5e-05 | 1.23e-05 | 1.22e-05 |
| `dynamics.nvt.T_kin` | K | >= 290 and <= 310 | >= 290 and <= 310 | 303 | 300 |
| `dynamics.nvt.T_conf_recip` | K | >= 285 and <= 315 | >= 285 and <= 315 | 304.9 | 302.5 |
| `dynamics.npt.T_kin` | K | >= 290 and <= 310 | >= 290 and <= 310 | 301.2 | 299.7 |
| `dynamics.npt.T_conf_recip` | K | >= 285 and <= 315 | >= 285 and <= 315 | 301.7 | 300.4 |

## Benchmark Z1: 6 at.% Zn, 300 K, NPT -> NVT (`sample_case/run.py`)

**NPT stage** (50,000 + 100,000 steps; reference: version 1.2 reference run (default options: quintic table interpolation, sample_stride 1; same seeds as the version 1.0 reference); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std | version 1.0 reference |
|---|---|---|---|---|---|---|
| L/L0 (mean box length / initial) |  | 1.0089 | +/- 1.7e-05 | 2.99e-06 | 2.89e-06 | 1.00892 |
| C44 | GPa | 72.0448 | +/- 0.32 | 0.052 | 0.055 | 71.8665 |
| C^B_11 (Born) | GPa | 163.42 | +/- 0.056 | 0.0054 | 0.0098 | 163.394 |
| C^B_12 (Born) | GPa | 108.202 | +/- 0.045 | 0.0045 | 0.0079 | 108.188 |
| C^B_44 (Born) | GPa | 76.0777 | +/- 0.025 | 0.0032 | 0.0044 | 76.0626 |
| T_K | K | 299.902 | +/- 0.79 | 0.14 | 0.14 | 300.172 |
| T_conf, 1/<1/T_conf> | K | 299.605 | +/- 2.8 | 0.33 | 0.49 | 300.632 |
| T_conf, <T_conf> | K | 300.193 | +/- 2.8 | 0.33 | 0.49 | 301.219 |
| mean pressure | GPa | -2.38e-05 | +/- 0.00097 | 0.00017 | 3.75e-05 | 9.16e-06 |
| C11 | GPa | 151.748 | not checked (INFO) | 11 | - | 146.765 |
| C12 | GPa | 101.321 | not checked (INFO) | 11 | - | 96.278 |
| NPT K term | GPa | 119.23 | not checked (INFO) | - | - | 116.736 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 4.5e-05 | one-sided | reference 1.32e-05 | max of set 1.47e-05 | - |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 2.4e-05 | one-sided | reference -2.29e-06 | max of set 7.93e-06 | - |

**NVT stage** (50,000 + 100,000 steps; reference: version 1.2 reference run (default options: quintic table interpolation, sample_stride 1; same seeds as the version 1.0 reference); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std | version 1.0 reference |
|---|---|---|---|---|---|---|
| C11 | GPa | 153.239 | +/- 1 | 0.18 | 0.11 | 153.066 |
| C12 | GPa | 102.63 | +/- 0.88 | 0.15 | 0.12 | 102.629 |
| C44 | GPa | 71.9473 | +/- 0.42 | 0.058 | 0.074 | 71.8933 |
| C^B_11 (Born) | GPa | 163.396 | +/- 0.069 | 0.0065 | 0.012 | 163.37 |
| C^B_12 (Born) | GPa | 108.182 | +/- 0.05 | 0.0054 | 0.0088 | 108.171 |
| C^B_44 (Born) | GPa | 76.0713 | +/- 0.033 | 0.0038 | 0.0058 | 76.053 |
| stress-fluctuation term, 11 | GPa | 10.8296 | +/- 1 | 0.18 | 0.11 | 10.9762 |
| stress-fluctuation term, 12 | GPa | 5.55184 | +/- 0.88 | 0.15 | 0.13 | 5.54157 |
| stress-fluctuation term, 44 | GPa | 4.46008 | +/- 0.44 | 0.058 | 0.076 | 4.49571 |
| kinetic (ideal-gas) term, 11 | GPa | 0.672065 | +/- 3.3e-05 | 0 | 5.78e-06 | 0.672028 |
| T_K | K | 300.064 | +/- 0.83 | 0.13 | 0.15 | 299.99 |
| T_conf, 1/<1/T_conf> | K | 299.484 | +/- 3.1 | 0.4 | 0.54 | 299.975 |
| T_conf, <T_conf> | K | 300.075 | +/- 3.1 | 0.4 | 0.53 | 300.577 |
| mean pressure | GPa | -0.00138464 | +/- 0.009 | 0.000998 | 0.0016 | -0.00381331 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 3.2e-05 | one-sided | reference 8.28e-06 | max of set 1.05e-05 | - |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 6.6e-06 | one-sided | reference 1.25e-06 | max of set 2.18e-06 | - |

## Benchmark L1: L1_2 CuZn3, 5 K, NPT -> NVT (`l12_cuzn3_5K/run.py`)

**NPT stage** (50,000 + 100,000 steps; reference: version 1.2 reference run (default options: quintic table interpolation, sample_stride 1; same seeds as the version 1.0 reference); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std | version 1.0 reference |
|---|---|---|---|---|---|---|
| L/L0 (mean box length / initial) |  | 1.00005 | +/- 6.7e-06 | 2.5e-07 | 1.17e-06 | 1.00005 |
| C44 | GPa | 42.8434 | +/- 0.0084 | 0.000635 | 0.0015 | 42.8412 |
| C^B_11 (Born) | GPa | 97.7266 | +/- 0.0089 | 0.000243 | 0.0016 | 97.7264 |
| C^B_12 (Born) | GPa | 60.8317 | +/- 0.0028 | 6.93e-05 | 0.00048 | 60.8313 |
| C^B_44 (Born) | GPa | 42.875 | +/- 0.0034 | 9.08e-05 | 0.000593 | 42.8745 |
| T_K | K | 4.99777 | +/- 0.015 | 0.0025 | 0.0025 | 5.00036 |
| T_conf, 1/<1/T_conf> | K | 5.17158 | +/- 0.23 | 0.0082 | 0.04 | 5.18936 |
| T_conf, <T_conf> | K | 5.18048 | +/- 0.23 | 0.0082 | 0.04 | 5.19944 |
| mean pressure | GPa | 3.73e-06 | +/- 2.3e-05 | 3.99e-06 | 1.63e-06 | 4.44e-06 |
| C11 | GPa | 50.3521 | not checked (INFO) | 7.5 | - | 46.051 |
| C12 | GPa | 13.5399 | not checked (INFO) | 7.5 | - | 9.2458 |
| NPT K term | GPa | 53.1902 | not checked (INFO) | - | - | 51.7296 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 3.9e-06 | one-sided | reference 1.29e-06 | max of set 1.29e-06 | - |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 3.4e-07 | one-sided | reference 4.98e-08 | max of set 1.11e-07 | - |

**NVT stage** (50,000 + 100,000 steps; reference: version 1.2 reference run (default options: quintic table interpolation, sample_stride 1; same seeds as the version 1.0 reference); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std | version 1.0 reference |
|---|---|---|---|---|---|---|
| C11 | GPa | 97.6246 | +/- 0.015 | 0.0025 | 0.0025 | 97.6332 |
| C12 | GPa | 60.8175 | +/- 0.023 | 0.0021 | 0.0039 | 60.819 |
| C44 | GPa | 42.8403 | +/- 0.007 | 0.000588 | 0.0012 | 42.8395 |
| C^B_11 (Born) | GPa | 97.7243 | +/- 0.0078 | 0.00017 | 0.0014 | 97.7271 |
| C^B_12 (Born) | GPa | 60.8287 | +/- 0.005 | 0.000135 | 0.000877 | 60.8289 |
| C^B_44 (Born) | GPa | 42.8719 | +/- 0.0049 | 0.00015 | 0.000861 | 42.872 |
| stress-fluctuation term, 11 | GPa | 0.109962 | +/- 0.019 | 0.0026 | 0.0033 | 0.104134 |
| stress-fluctuation term, 12 | GPa | 0.0112295 | +/- 0.025 | 0.0021 | 0.0042 | 0.00993313 |
| stress-fluctuation term, 44 | GPa | 0.0366504 | +/- 0.0068 | 0.000579 | 0.0012 | 0.037578 |
| kinetic (ideal-gas) term, 11 | GPa | 0.0102107 | +/- 2.1e-07 | 0 | 3.58e-08 | 0.0102107 |
| T_K | K | 4.99919 | +/- 0.011 | 0.0019 | 0.001 | 5.0012 |
| T_conf, 1/<1/T_conf> | K | 5.08325 | +/- 0.2 | 0.014 | 0.035 | 4.99842 |
| T_conf, <T_conf> | K | 5.09146 | +/- 0.2 | 0.014 | 0.034 | 5.00839 |
| mean pressure | GPa | -0.00105246 | +/- 0.0023 | 5.56e-05 | 0.000403 | -0.001117 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 1.8e-07 | one-sided | reference 4.13e-08 | max of set 5.99e-08 | - |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 9.2e-08 | one-sided | reference -3.06e-08 | max of set 3.06e-08 | - |

Additional check: |C(5 K) - C(0 K)| / C(0 K) <= 0.002 against the 0 K reference values C11 = 97.638, C12 = 60.838, C44 = 42.885 GPa (NIST Interatomic Potentials Repository, Clement-Auger Cu-Zn potential (manuscript Table 1)). Table 1 lists the nine-replica means C11 = 97.634, C12 = 60.824, C44 = 42.841 GPa.
