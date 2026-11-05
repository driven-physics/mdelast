#!/usr/bin/env python3
"""Benchmark S (smoke test, a few minutes on one core).

Part A -- static finite-difference checks on one fixed, distorted 256-atom
Cu-Zn cell (no dynamics, deterministic up to floating-point round-off):
  A1  net force vanishes;
  A2  analytic forces F = -dU/dr against central differences of the energy;
  A3  Born (configurational) stress against (1/V) dU/deta_aa;
  A4  Born elastic constants C^B_ab (a,b = xx,yy,zz) against central
      differences of the analytic Born stress with respect to the Lagrangian
      strain (orthorhombic cell);
  A5  configurational-temperature Laplacian sum_i lap_i U against central
      differences of the analytic forces (every atom, every direction);
  A6  absolute reference values at this configuration (energy, pressure,
      Born constants, Laplacian, T_conf) for regression across platforms.
Part B -- short dynamics (500 atoms, 6 at.% Zn, 300 K):
  B1  NVE: conserved-energy error at dt = 1 and 2 fs and its ~dt^2 scaling;
  B2  NVT and NPT: extended-Hamiltonian error and temperatures stay sane.

Shear components (44, 55, 66) cannot be strained in mdelast's orthorhombic
cell; they enter benchmarks Z1 and L1 through C44.

Usage (from the repository root, after building):
    python3 benchmarks/smoke/run_smoke.py --build-dir build/Release \
        --output-dir benchmark_output/smoke
    python3 benchmarks/verify.py smoke benchmark_output/smoke
"""

from __future__ import annotations

import argparse
import contextlib
import ctypes
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo-root", type=Path, default=HERE.parents[1],
                   help="mdelast repository root (contains src/python and potentials/)")
    p.add_argument("--build-dir", type=Path, default=None,
                   help="directory with md_driver_module*.so (default: <repo>/build/Release)")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--skip-dynamics", action="store_true", help="run Part A only")
    return p.parse_args()


_LIBC = ctypes.CDLL(None)


@contextlib.contextmanager
def quiet_stdout():
    """Silence the C++ printf output of the driver constructor (fd 1)."""
    sys.stdout.flush()
    _LIBC.fflush(None)
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    try:
        yield
    finally:
        _LIBC.fflush(None)
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


# ---------------------------------------------------------------------------
# Part A: static checks
# ---------------------------------------------------------------------------
A0_ANG = 3.6345          # production lattice constant of 6 at.% Zn [A]
NSIDE_STATIC = 4         # 256 atoms
N_ZN_STATIC = 15         # ~6 at.% Zn
DISPLACEMENT_ANG = 0.08  # Gaussian random displacement (breaks symmetry)
SEED_STATIC = 2026
H_POS_ANG = 1e-4         # finite-difference step for positions [A]
H_STRAIN = 1e-4          # finite-difference step for Lagrangian strain (energy)
H_STRAIN_BORN = 1e-5     # finite-difference step for Lagrangian strain (Born stress)
N_FORCE_ATOMS = 16       # atoms used in the force check (all 3 directions)


def fcc_positions(nside, a):
    basis = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    pos = np.array([(np.array([i, j, k]) + 0.25 + b) * a
                    for i in range(nside) for j in range(nside)
                    for k in range(nside) for b in basis])
    return pos - nside * a / 2.0


class StaticCell:
    """Build MDDriver objects for one configuration (cgs units inside)."""

    def __init__(self, mdm, get_potential):
        self.m = mdm
        self.U = mdm.Units
        elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy")
        self.elements = [e[0] for e in elements]
        self.mass = [e[2] / self.U.Navogadro for e in elements]
        self.pot = mdm.PotentialVector(potentials)
        self.eam = mdm.EAMVector(eams)
        rs = np.random.RandomState(SEED_STATIC)
        pos = fcc_positions(NSIDE_STATIC, A0_ANG)
        pos = pos + rs.normal(scale=DISPLACEMENT_ANG, size=pos.shape)
        self.L = NSIDE_STATIC * A0_ANG
        self.pos = np.mod(pos + self.L / 2, self.L) - self.L / 2          # [A]
        types = np.full(len(pos), self.elements.index("Cu"))
        types[rs.choice(len(pos), size=N_ZN_STATIC, replace=False)] = self.elements.index("Zn")
        self.types = [int(t) for t in types]
        self.N = len(pos)

    def driver(self, pos_ang=None, lengths_ang=None):
        m, U = self.m, self.U
        pos = self.pos if pos_ang is None else pos_ang
        lengths = np.full(3, self.L) if lengths_ang is None else np.asarray(lengths_ang)
        pos = np.mod(pos + lengths / 2, lengths) - lengths / 2
        P = [m.Vector3d(*(p * U.angstrom)) for p in pos]
        V = [m.Vector3d(0.0, 0.0, 0.0) for _ in pos]
        th = m.Thermostat(300.0, [1e-30], [0.0], 3.0 * self.N)
        bth = m.Thermostat(300.0, [1e-30], [0.0], 1.0)
        baro = m.Barostat(0.0, 1e-30, 0.0, bth)
        box = m.Box(np.diag(lengths * U.angstrom))
        drg = A0_ANG / 500 * U.angstrom
        rg = [(k + 0.5) * drg for k in range(10)]
        with quiet_stdout():
            d = m.MDDriver(P, V, self.pot, self.eam, m.MDDriver.NVE,
                           m.MDDriver.THERMOSTAT_KINETIC, th, baro, box, rg,
                           3.0 * self.N, self.mass, self.types, 1e-15, False, 1.0, 1e3)
        return d

    def energy(self, d):
        return d.epot_atom * self.N  # erg

    @staticmethod
    def forces(d):
        return np.array([[f.x, f.y, f.z] for f in d.forces])  # dyn


def run_static(mdm, get_potential):
    U = mdm.Units
    cell = StaticCell(mdm, get_potential)
    d0 = cell.driver()
    N, L = cell.N, cell.L
    V = (L * U.angstrom) ** 3
    E0 = cell.energy(d0)
    F0 = cell.forces(d0)
    res = {"N": N, "n_Zn": int(sum(1 for t in cell.types if cell.elements[t] == "Zn")),
           "box_A": L}
    ev_a = U.ev / U.angstrom

    # A1 net force
    res["net_force_rel"] = float(np.abs(F0.sum(0)).max() / np.abs(F0).max())

    # A2 forces by central differences of the energy
    rs = np.random.RandomState(1)
    atoms = rs.choice(N, size=N_FORCE_ATOMS, replace=False)
    fd, an = [], []
    for i in atoms:
        for a in range(3):
            e = []
            for s in (+1, -1):
                p = cell.pos.copy()
                p[i, a] += s * H_POS_ANG
                e.append(cell.energy(cell.driver(p)))
            fd.append(-(e[0] - e[1]) / (2 * H_POS_ANG * U.angstrom))
            an.append(F0[i, a])
    fd, an = np.array(fd), np.array(an)
    res["force_fd_max_abs_err_eV_A"] = float(np.abs(fd - an).max() / ev_a)
    res["force_fd_max_rel_err"] = float(np.abs(fd - an).max() / np.abs(an).max())
    res["force_max_abs_eV_A"] = float(np.abs(F0).max() / ev_a)

    # A3 Born stress = (1/V) dU/deta_a, central differences of the energy
    def strained(eta):
        lam = np.sqrt(1.0 + 2.0 * np.asarray(eta))   # exact Lagrangian strain eta_a
        return cell.driver(cell.pos * lam, L * lam), lam

    born_stress_an = np.array(d0.born_stress_tensor).ravel()[:3]
    dstress = []
    for a in range(3):
        ep = np.zeros(3)
        ep[a] = H_STRAIN
        dstress.append((cell.energy(strained(ep)[0]) - cell.energy(strained(-ep)[0]))
                       / (2 * H_STRAIN) / V)
    dstress = np.array(dstress)
    res["born_stress_GPa"] = (born_stress_an / 1e10).tolist()
    res["born_stress_fd_GPa"] = (dstress / 1e10).tolist()
    res["born_stress_fd_max_abs_err_GPa"] = float(np.abs(dstress - born_stress_an).max() / 1e10)

    # A4 Born constants C^B_ab = (1/V0) d T_a / d eta_b, where
    #    T_a(eta) = V(eta) sigma^B_a(eta) / lambda_a^2 = dU/deta_a (reference frame),
    #    i.e. central differences of the analytic Born stress.
    def t_ref(eta):
        d, lam = strained(eta)
        vol = V * np.prod(lam)
        return vol * np.array(d.born_stress_tensor).ravel()[:3] / lam ** 2

    CB_an = np.array(d0.born_elasticity_tensor)
    CB_an = np.where(CB_an == 0.0, CB_an.T, CB_an)       # upper triangle is stored
    M = np.zeros((3, 3))
    for b in range(3):
        ep = np.zeros(3)
        ep[b] = H_STRAIN_BORN
        M[:, b] = (t_ref(ep) - t_ref(-ep)) / (2 * H_STRAIN_BORN) / V
    M = 0.5 * (M + M.T)
    comps = {}
    for a, b in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
        comps[f"C{a + 1}{b + 1}"] = (float(CB_an[a, b] / 1e10), float(M[a, b] / 1e10))
    res["born_C_analytic_vs_fd_GPa"] = comps
    res["born_C_fd_max_abs_err_GPa"] = float(max(abs(x - y) for x, y in comps.values()))
    res["born_C_fd_max_rel_err"] = float(max(abs(x - y) / abs(x) for x, y in comps.values()))

    # A5 Laplacian from central differences of the forces
    lap_fd = 0.0
    hc = H_POS_ANG * U.angstrom
    for i in range(N):
        for a in range(3):
            f = []
            for s in (+1, -1):
                p = cell.pos.copy()
                p[i, a] += s * H_POS_ANG
                f.append(cell.forces(cell.driver(p))[i, a])
            lap_fd += -(f[0] - f[1]) / (2 * hc)
    lap_an = d0.total_laplacian
    res["laplacian_eV_A2"] = float(lap_an / (U.ev / U.angstrom ** 2))
    res["laplacian_fd_eV_A2"] = float(lap_fd / (U.ev / U.angstrom ** 2))
    res["laplacian_fd_rel_err"] = float(abs(lap_fd - lap_an) / abs(lap_an))

    # A6 absolute reference values
    CB = CB_an / 1e10
    res["static_epot_eV_per_atom"] = float(d0.epot_atom / U.ev)
    res["static_pressure_GPa"] = float(d0.pressure / 1e10)
    res["static_CB11_GPa"] = float(np.mean([CB[0, 0], CB[1, 1], CB[2, 2]]))
    res["static_CB12_GPa"] = float(np.mean([CB[0, 1], CB[0, 2], CB[1, 2]]))
    res["static_CB44_GPa"] = float(np.mean([CB[3, 3], CB[4, 4], CB[5, 5]]))
    res["static_T_conf_K"] = float(1.0 / d0.temperature_conf_inverse)
    return res


# ---------------------------------------------------------------------------
# Part B: short dynamics through the documented Python API
# ---------------------------------------------------------------------------
NSIDE_DYN = 5           # 500 atoms
N_ZN_DYN = 30           # 6 at.% Zn
SEED_PLACEMENT = 42     # same Zn arrangement as sample_case/run.py
T_DYN = 300.0


def run_dynamics(mdm, md_mod, get_potential, out_dir, ensemble, dt_fs, n_th, n_prod):
    U = mdm.Units
    elements, potentials, eams = get_potential.get_EAM_alloy_file("ZnCu.eam.alloy")
    pos = fcc_positions(NSIDE_DYN, A0_ANG)
    L = NSIDE_DYN * A0_ANG
    rs = np.random.RandomState(SEED_PLACEMENT)
    types = [1] * len(pos)
    for i in rs.choice(len(pos), size=N_ZN_DYN, replace=False):
        types[int(i)] = 0
    mat = md_mod.Material(potential=potentials, atomic_weight=[e[2] for e in elements],
                          atom_name="-".join(e[0] for e in elements), eam=eams,
                          lattice=[mdm.Vector3d(*(p * U.angstrom)) for p in pos],
                          atom_type=types, box_length=(L * U.angstrom,) * 3)
    t0 = time.time()
    md = md_mod.MolecularDynamics(
        work_dir=str(out_dir), ensamble=ensemble, material=mat,
        dt_global=dt_fs * 1e-15, tau_thermostat=5e-14, tau_barostat=5e-13,
        temperature_initial=T_DYN, temperature_thermostat=T_DYN,
        pressure_barostat=1e-10 * U.GPa, dump_cycles=10 ** 9, slurm_plt=True,
        show_graphs=False, plot_graphs=False, fatal_conservation=1e3,
        cycles_thermalization=n_th, plot_freq=n_prod, cycles_end=n_th + n_prod,
        neighbors=True, r_skin=1.2, cycles_neigbors=20, gr_stride=100, random_seed=0,
    )
    if ensemble == mdm.MDDriver.NVE:
        # MolecularDynamics.run() fails at the final NVE averaging step in the
        # current code (np.polyfit on TimeSeries objects), so the NVE check
        # drives the compiled integrator directly, with the same neighbour-list
        # schedule as run().
        d = md.driver
        econs, tk, tci = [], [], []
        for step in range(n_th + n_prod + 1):
            if step % 20 == 0:
                d.calculate_neighbors()
            d.integrate()
            econs.append(d.etot_cons)
            tk.append(d.temperature_kin)
            tci.append(d.temperature_conf_inverse)
        n0 = n_th
        T_kin = float(np.mean(tk[n0:]))
        T_conf = float(1.0 / np.mean(tci[n0:]))
    else:
        md.run()
        s = md.series
        econs = s.etot_cons.data
        T_kin = float(s.temperature_kin.mean)
        T_conf = float(1.0 / s.temperature_conf_inverse.mean)
    wall = time.time() - t0
    econs = np.asarray(econs, dtype=float)
    rel = econs / econs[0] - 1.0
    half = rel[len(rel) // 2:]
    return {
        "dt_fs": dt_fs, "steps": n_th + n_prod,
        "econs_err_max": float(np.abs(rel).max()),
        "econs_err_rms_second_half": float(np.sqrt(np.mean((half - half.mean()) ** 2))),
        "T_kin": T_kin,
        "T_conf_recip": T_conf,
        "wall_s": wall, "ms_per_step": 1e3 * wall / (n_th + n_prod),
    }


def main():
    args = parse_args()
    repo = args.repo_root.resolve()
    build = (args.build_dir or repo / "build" / "Release").resolve()
    out = args.output_dir.resolve()
    if out.exists():
        raise SystemExit(f"output directory already exists: {out}")
    out.mkdir(parents=True)
    sys.path.insert(0, str(build))
    sys.path.insert(0, str(repo / "src" / "python"))
    import md_driver_module as mdm
    import get_potential
    import molecular_dynamics as md_mod

    result = {"benchmark": "smoke", "build_dir": str(build)}
    t0 = time.time()
    result["static"] = run_static(mdm, get_potential)
    result["static"]["wall_s"] = time.time() - t0
    print(json.dumps(result["static"], indent=1), flush=True)

    if not args.skip_dynamics:
        dyn = {}
        # B1: NVE, same physical time (2 ps) at two timesteps
        dyn["nve_dt2"] = run_dynamics(mdm, md_mod, get_potential, out / "nve_dt2",
                                      mdm.MDDriver.NVE, 2.0, 500, 500)
        dyn["nve_dt1"] = run_dynamics(mdm, md_mod, get_potential, out / "nve_dt1",
                                      mdm.MDDriver.NVE, 1.0, 1000, 1000)
        dyn["nve_dt_ratio"] = dyn["nve_dt2"]["econs_err_rms_second_half"] / \
            dyn["nve_dt1"]["econs_err_rms_second_half"]
        # B2: thermostatted / barostatted, production settings, 4 ps each
        dyn["nvt"] = run_dynamics(mdm, md_mod, get_potential, out / "nvt",
                                  mdm.MDDriver.NVT, 2.0, 1000, 1000)
        dyn["npt"] = run_dynamics(mdm, md_mod, get_potential, out / "npt",
                                  mdm.MDDriver.NPT, 2.0, 1000, 1000)
        result["dynamics"] = dyn
        print(json.dumps(dyn, indent=1), flush=True)

    result["wall_s_total"] = time.time() - t0
    (out / "result.json").write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {out / 'result.json'}")


if __name__ == "__main__":
    main()
