# Expected benchmark values

Generated from `expected/*.json` (the files `verify.py` reads). All stresses and elastic constants in GPa, temperatures in K. Cubic averages: C11 = (C11+C22+C33)/3, C12 = (C12+C13+C23)/3, C44 = (C44+C55+C66)/3.

Tolerance rule for Z1 and L1: tol = 4 x sqrt(2) x max(per-run SE, velocity-replica std); conserved-energy checks are upper bounds of 3 x the largest value in the reference set. See README.md, Section 'How the tolerances were set'.

## Benchmark S: smoke test (`smoke/run_smoke.py`)

| check | unit | expected / bound | observed (public code) |
|---|---|---|---|
| `static.net_force_rel` |  | <= 1e-12 | 2.93e-15 |
| `static.force_fd_max_abs_err_eV_A` | eV/A | <= 2e-05 | 3.68e-06 |
| `static.born_stress_fd_max_abs_err_GPa` | GPa | <= 1e-05 | 1.21e-06 |
| `static.born_C_fd_max_abs_err_GPa` | GPa | <= 0.05 | 0.01288 |
| `static.laplacian_fd_rel_err` |  | <= 2e-05 | 3.69e-06 |
| `static.static_epot_eV_per_atom` | eV | -3.49917801 +/- 3.5e-07 | = expected |
| `static.static_pressure_GPa` | GPa | 5.131072142 +/- 5.1e-07 | = expected |
| `static.static_CB11_GPa` | GPa | 197.396401 +/- 2e-05 | = expected |
| `static.static_CB12_GPa` | GPa | 126.5757319 +/- 1.3e-05 | = expected |
| `static.static_CB44_GPa` | GPa | 92.13965579 +/- 9.2e-06 | = expected |
| `static.laplacian_eV_A2` | eV/A^2 | 6567.770385 +/- 0.00066 | = expected |
| `static.static_T_conf_K` | K | 1000.120487 +/- 0.0001 | = expected |
| `dynamics.nve_dt_ratio` |  | >= 3.5 and <= 4.5 | 3.998 |
| `dynamics.nve_dt2.econs_err_max` |  | <= 5e-05 | 1.06e-05 |
| `dynamics.nve_dt1.econs_err_max` |  | <= 1.5e-05 | 2.65e-06 |
| `dynamics.nvt.econs_err_max` |  | <= 5e-05 | 1.37e-05 |
| `dynamics.npt.econs_err_max` |  | <= 5e-05 | 1.22e-05 |
| `dynamics.nvt.T_kin` | K | >= 290 and <= 310 | 300 |
| `dynamics.nvt.T_conf_recip` | K | >= 285 and <= 315 | 302.5 |
| `dynamics.npt.T_kin` | K | >= 290 and <= 310 | 299.7 |
| `dynamics.npt.T_conf_recip` | K | >= 285 and <= 315 | 300.4 |

## Benchmark Z1: 6 at.% Zn, 300 K, NPT -> NVT (`sample_case/run.py`)

**NPT stage** (50,000 + 100,000 steps; reference: the run of manuscript Sec. 5.1 (Figs. 1-5); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std |
|---|---|---|---|---|---|
| L/L0 (mean box length / initial) |  | 1.00892 | +/- 1.7e-05 | 2.99e-06 | 2.89e-06 |
| C44 | GPa | 71.8665 | +/- 0.32 | 0.052 | 0.055 |
| C^B_11 (Born) | GPa | 163.394 | +/- 0.056 | 0.0054 | 0.0098 |
| C^B_12 (Born) | GPa | 108.188 | +/- 0.045 | 0.0045 | 0.0079 |
| C^B_44 (Born) | GPa | 76.0626 | +/- 0.025 | 0.0032 | 0.0044 |
| T_K | K | 300.172 | +/- 0.79 | 0.14 | 0.14 |
| T_conf, 1/<1/T_conf> | K | 300.632 | +/- 2.8 | 0.33 | 0.49 |
| T_conf, <T_conf> | K | 301.219 | +/- 2.8 | 0.33 | 0.49 |
| mean pressure | GPa | 9.16e-06 | +/- 0.00097 | 0.00017 | 3.75e-05 |
| C11 | GPa | 146.765 | not checked (decision (b)) | 11 | - |
| C12 | GPa | 96.278 | not checked (decision (b)) | 11 | - |
| NPT K term | GPa | 116.736 | not checked (decision (b)) | - | - |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 4.5e-05 | one-sided | reference 1.32e-05 | max of set 1.47e-05 |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 2.4e-05 | one-sided | reference -2.29e-06 | max of set 7.93e-06 |

**NVT stage** (50,000 + 100,000 steps; reference: the run of manuscript Sec. 5.1 (Figs. 1-5); 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std |
|---|---|---|---|---|---|
| C11 | GPa | 153.066 | +/- 1 | 0.18 | 0.11 |
| C12 | GPa | 102.629 | +/- 0.88 | 0.15 | 0.12 |
| C44 | GPa | 71.8933 | +/- 0.42 | 0.058 | 0.074 |
| C^B_11 (Born) | GPa | 163.37 | +/- 0.069 | 0.0065 | 0.012 |
| C^B_12 (Born) | GPa | 108.171 | +/- 0.05 | 0.0054 | 0.0088 |
| C^B_44 (Born) | GPa | 76.053 | +/- 0.033 | 0.0038 | 0.0058 |
| stress-fluctuation term, 11 | GPa | 10.9762 | +/- 1 | 0.18 | 0.11 |
| stress-fluctuation term, 12 | GPa | 5.54157 | +/- 0.88 | 0.15 | 0.13 |
| stress-fluctuation term, 44 | GPa | 4.49571 | +/- 0.44 | 0.058 | 0.076 |
| kinetic (ideal-gas) term, 11 | GPa | 0.672028 | +/- 3.3e-05 | 0 | 5.78e-06 |
| T_K | K | 299.99 | +/- 0.83 | 0.13 | 0.15 |
| T_conf, 1/<1/T_conf> | K | 299.975 | +/- 3.1 | 0.4 | 0.54 |
| T_conf, <T_conf> | K | 300.577 | +/- 3.1 | 0.4 | 0.53 |
| mean pressure | GPa | -0.00381331 | +/- 0.009 | 0.000998 | 0.0016 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 3.2e-05 | one-sided | reference 8.28e-06 | max of set 1.05e-05 |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 6.6e-06 | one-sided | reference 1.25e-06 | max of set 2.18e-06 |

## Benchmark L1: L1_2 CuZn3, 5 K, NPT -> NVT (`l12_cuzn3_5K/run.py`)

**NPT stage** (50,000 + 100,000 steps; reference: replica 0 of the set behind Table tab:elastic-avg; 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std |
|---|---|---|---|---|---|
| L/L0 (mean box length / initial) |  | 1.00005 | +/- 6.7e-06 | 2.5e-07 | 1.17e-06 |
| C44 | GPa | 42.8412 | +/- 0.0084 | 0.000635 | 0.0015 |
| C^B_11 (Born) | GPa | 97.7264 | +/- 0.0089 | 0.000243 | 0.0016 |
| C^B_12 (Born) | GPa | 60.8313 | +/- 0.0028 | 6.93e-05 | 0.00048 |
| C^B_44 (Born) | GPa | 42.8745 | +/- 0.0034 | 9.08e-05 | 0.000593 |
| T_K | K | 5.00036 | +/- 0.015 | 0.0025 | 0.0025 |
| T_conf, 1/<1/T_conf> | K | 5.18936 | +/- 0.23 | 0.0082 | 0.04 |
| T_conf, <T_conf> | K | 5.19944 | +/- 0.23 | 0.0082 | 0.04 |
| mean pressure | GPa | 4.44e-06 | +/- 2.3e-05 | 3.99e-06 | 1.63e-06 |
| C11 | GPa | 46.051 | not checked (decision (b)) | 7.5 | - |
| C12 | GPa | 9.2458 | not checked (decision (b)) | 7.5 | - |
| NPT K term | GPa | 51.7296 | not checked (decision (b)) | - | - |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 3.9e-06 | one-sided | reference 1.29e-06 | max of set 1.29e-06 |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 3.4e-07 | one-sided | reference 4.98e-08 | max of set 1.11e-07 |

**NVT stage** (50,000 + 100,000 steps; reference: replica 0 of the set behind Table tab:elastic-avg; 9 velocity replicas)

| quantity | unit | expected | tolerance | per-run SE | velocity-replica std |
|---|---|---|---|---|---|
| C11 | GPa | 97.6332 | +/- 0.015 | 0.0025 | 0.0025 |
| C12 | GPa | 60.819 | +/- 0.023 | 0.0021 | 0.0039 |
| C44 | GPa | 42.8395 | +/- 0.007 | 0.000588 | 0.0012 |
| C^B_11 (Born) | GPa | 97.7271 | +/- 0.0078 | 0.00017 | 0.0014 |
| C^B_12 (Born) | GPa | 60.8289 | +/- 0.005 | 0.000135 | 0.000877 |
| C^B_44 (Born) | GPa | 42.872 | +/- 0.0049 | 0.00015 | 0.000861 |
| stress-fluctuation term, 11 | GPa | 0.104134 | +/- 0.019 | 0.0026 | 0.0033 |
| stress-fluctuation term, 12 | GPa | 0.00993313 | +/- 0.025 | 0.0021 | 0.0042 |
| stress-fluctuation term, 44 | GPa | 0.037578 | +/- 0.0068 | 0.000579 | 0.0012 |
| kinetic (ideal-gas) term, 11 | GPa | 0.0102107 | +/- 2.1e-07 | 0 | 3.58e-08 |
| T_K | K | 5.0012 | +/- 0.011 | 0.0019 | 0.001 |
| T_conf, 1/<1/T_conf> | K | 4.99842 | +/- 0.2 | 0.014 | 0.035 |
| T_conf, <T_conf> | K | 5.00839 | +/- 0.2 | 0.014 | 0.034 |
| mean pressure | GPa | -0.001117 | +/- 0.0023 | 5.56e-05 | 0.000403 |
| max abs(E_cons(t)/E_cons(0) - 1) |  | <= 1.8e-07 | one-sided | reference 4.13e-08 | max of set 5.99e-08 |
| abs(linear drift of E_cons/E_cons(0)) over production | 1/ns | <= 9.2e-08 | one-sided | reference -3.06e-08 | max of set 3.06e-08 |

Additional check: |C(5 K) - C(0 K)| / C(0 K) <= 0.002 against the 0 K reference values {'C11': 97.638, 'C12': 60.838, 'C44': 42.885} (NIST Interatomic Potentials Repository, Clement-Auger Cu-Zn potential (manuscript Table tab:elastic-avg)). Table tab:elastic-avg lists the nine-replica means {'C11': 97.634, 'C12': 60.824, 'C44': 42.841}.
