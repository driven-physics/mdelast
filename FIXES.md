# Fixes in this version

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
