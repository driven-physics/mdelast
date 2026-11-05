"""High-level Python interface for the molecular-dynamics engine.

The compiled :mod:`md_driver_module` extension performs integration and force
evaluation.  This module prepares materials and lattices, manages simulations,
collects time series, computes thermodynamic averages, and writes HDF5 output.
Internal quantities use the CGS unit constants exposed by ``Units``.
"""

import pickle
import h5py
import numpy as np
import socket
import time
import os
from os import path
import sys
from scipy import integrate
from scipy.interpolate import UnivariateSpline

try:
    import md_driver_module
    from md_driver_module import (
        Box,
        Vector3d,
        MDDriver,
        Thermostat,
        Barostat,
        BarostatFlex,
        Units,
        Potential,
        EAM,
        PotentialVector,
        EAMVector,
    )
except:
    pass

from averages import *
from series import *


import logging

if len(logging.root.handlers) == 0:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("MD")


class Material:
    """Describe a single- or multi-species material used by a simulation.

    Parameters are accepted either as scalars or per-species sequences.
    ``lattice`` and ``box_length`` may be supplied for an explicit structure;
    otherwise :class:`MolecularDynamics` constructs an FCC lattice from
    ``a_lattice`` and ``N_lattice_side``.
    """

    def __init__(
        self,
        potential,
        atomic_weight,
        atom_name,
        eam=None,
        a_lattice=None,
        lattice=None,
        atom_type=None,
        box_length=None,
    ):
        if isinstance(potential, (list, tuple, np.ndarray)):
            self.potential = list(potential)
        else:
            self.potential = [potential]

        if isinstance(atomic_weight, (list, tuple, np.ndarray)):
            self.atomic_weight = list(atomic_weight)
        else:
            self.atomic_weight = [atomic_weight]
        self.atom_name = atom_name
        # Preserve "no EAM" as None (many code paths key off this),
        # otherwise treat EAM as a per-type list.
        if eam is None:
            self.eam = None
        elif isinstance(eam, (list, tuple, np.ndarray)):
            self.eam = list(eam)
        else:
            self.eam = [eam]

        self.a_lattice = a_lattice
        self.lattice = lattice
        if lattice is None:
            self.atom_type = [0]
            self.box_length = None
        else:
            assert (atom_type is not None) and (len(atom_type) == len(lattice))
            self.atom_type = atom_type
            assert box_length is not None
            self.box_length = box_length

    def __str__(self):
        return (
            f"Material: '{self.atom_name}', "
            f"A={self.atomic_weight}, "
            f"potential: '{self.potential}'"
            f"EAM: '{self.eam}'"
        )


def get_deriv(x, y):
    """Return midpoint coordinates and first finite-difference derivatives."""

    return 0.5 * (x[1:] + x[:-1]), np.diff(y) / np.diff(x)


def get_derivatives(x, y):
    """Evaluate first and second derivatives of an exact quintic spline."""

    spl = UnivariateSpline(x=x, y=y, k=5, s=0)
    dydx = spl.derivative(n=1)(x)
    d2ydx2 = spl.derivative(n=2)(x)
    return dydx, d2ydx2


def get_EAM_file(*, file_name, a_factor=1.0):
    """Load a supported single-element EAM potential from ``potentials/``.

    ``file_name`` must use the ``eam.alloy.txt`` (SETFL) or ``eam.fs.txt``
    suffix expected by the historical parser. ``a_factor`` scales only the
    lattice constant stored in the returned :class:`Material`.
    """

    # SETFL potential files taken from:
    # https://www.ctcms.nist.gov/potentials/entry/2001--Mishin-Y-Mehl-M-J-Papaconstantopoulos-D-A-et-al--Cu-1/2001--Mishin-Y--Cu-1--LAMMPS--ipr1.html
    # FS Mendelev file taken from
    # https://www.ctcms.nist.gov/potentials/entry/2013--Mendelev-M-I-King-A-H--Cu/2012--Mendelev-M-I--Cu--LAMMPS--ipr1.html
    # format documentation given in https://lammps.sandia.gov/doc/pair_eam.html

    pot_dir = path.join(path.dirname(path.abspath(__file__)), "../../potentials")
    file_path = path.join(pot_dir, file_name)
    logger.info(f"Reading file potential... '{file_path}'")

    # Read the five-line header.
    with open(file_path, "r") as f:
        line = f.readlines()[5].split("\n")[0].split()[:-1]
        atomic_number, atomic_weight, a_lattice = [float(s) for s in line]
        a_lattice *= Units.angstrom
        logger.info(
            f"Z={atomic_number:g}, ATW={atomic_weight:g}, a={a_lattice/Units.angstrom:g}[A]"
        )

    Nrho, drho, Nr, dr, rcut = np.loadtxt(file_path, skiprows=4, max_rows=1)[:]
    Nrho = int(Nrho)
    Nr = int(Nr)

    r_grid = np.linspace(0.0, dr * (Nr - 1.0), Nr)  # angstrom
    rho_grid = np.linspace(0.0, drho * (Nrho - 1.0), Nrho)  # angstrom**-3

    # read the rest of the file
    data = np.loadtxt(file_path, skiprows=6)

    if file_name[-13:] == "eam.alloy.txt":
        # 'SETFL' format (see https://lammps.sandia.gov/doc/pair_eam.html)
        assert len(data.shape) == 1 and data.shape[0] == Nrho + Nr + Nr
        F = data[:Nrho]  # embedding function F(rho) (Nrho values)
        rho = data[Nrho : Nrho + Nr]  # density function rho(r) (Nr values)
        phir = data[
            Nrho + Nr :
        ]  # phi function are listed in setfl files directly as r*phi (in units of eV-Angstroms)
    elif file_name[-10:] == "eam.fs.txt":
        # 'FS' format (see https://lammps.sandia.gov/doc/pair_eam.html)
        assert data.shape[1] == 5
        assert data.shape[0] * 5 == Nrho + Nr + Nr
        assert Nrho % 5 == 0
        assert Nr % 5 == 0
        F = []
        rho = []
        phir = []

        for i in range(Nrho // 5):
            F += list(data[i])
        F = np.array(F)
        for i in range(Nrho // 5, Nrho // 5 + Nr // 5):
            rho += list(data[i])
        assert set(rho[0:167]) == {
            1e12
        }  # Mendelev_Cu2_2012.eam.fs.txt fills those values with 1e12
        rho = np.array(rho)
        rho[0:167] = rho[167]

        for i in range(Nrho // 5 + Nr // 5, Nrho // 5 + 2 * Nr // 5):
            phir += list(data[i])
        assert set(phir[0:167]) == {
            1e12
        }  # Mendelev_Cu2_2012.eam.fs.txt fills those values with 1e12
        phir = np.array(phir)
        phir[0:167] = phir[167]
    else:
        logger.fatal(f"undefined potential format '{file_name}'")
        sys.exit(1)

    # Convert tabulated quantities to CGS.
    F *= Units.ev
    phi = phir / r_grid * Units.ev
    phi[0] = phi[1]

    dr *= Units.angstrom
    rcut *= Units.angstrom
    r_grid *= Units.angstrom

    drho /= Units.angstrom**3
    rho /= Units.angstrom**3
    rho_grid /= Units.angstrom**3

    rmin_fatal = 1e-2 * Units.angstrom  # below this radius the interpolation is unreliable
    dphidr, d2phidr2 = get_derivatives(x=r_grid, y=phi)
    drhodr, d2rhodr2 = get_derivatives(x=r_grid, y=rho)
    dFdrho, d2Fdrho = get_derivatives(x=rho_grid, y=F)

    logger.info(
        f"a_lattice={a_lattice*a_factor/Units.angstrom:g}[A] "
        f"rcut={rcut/Units.angstrom:g}[A] "
        f"r_grid[-1]={r_grid[-1]/Units.angstrom:g}[A]"
    )

    return Material(
        atomic_weight=atomic_weight,
        atom_name="Cu",
        potential=md_driver_module.TablePotential(
            r_grid[-1], rmin_fatal, dr, phi, dphidr, d2phidr2
        ),
        eam=md_driver_module.EAM(
            r_grid[-1], rho_grid[-1], dr, drho, rho, drhodr, d2rhodr2, F, dFdrho, d2Fdrho
        ),
        a_lattice=a_lattice * a_factor,
    )


def get_Cu_Table_Morse():
    """Return copper represented by a tabulated Morse pair potential."""

    morse = get_Cu_Morse().potential
    rcut = morse.rcut
    Ninterp = 10000
    r_grid = np.linspace(0.0, rcut, Ninterp)
    dr = r_grid[1] - r_grid[0]
    phi_morse = np.array([morse.calculate(r).V for r in r_grid])
    dphi_morse = np.array([morse.calculate(r).dVdr for r in r_grid])
    d2phi_morse = np.array([morse.calculate(r).d2Vdr2 for r in r_grid])

    return Material(
        potential=md_driver_module.TablePotential(
            r_grid[-1], 0.0, dr, phi_morse, dphi_morse, d2phi_morse
        ),
        atomic_weight=63.546,
        atom_name="Cu",
        a_lattice=3.615 * Units.angstrom,
    )


def get_Cu_Morse(a_factor=1.0):
    """Return copper represented by the built-in Morse pair potential."""

    # Morse parameters for copper.
    from md_driver_module import Morse

    D = 0.3429 * Units.ev
    alpha = 1.3588 / Units.angstrom
    r0 = 2.866 * Units.angstrom
    rcut = 10.0 / alpha + r0
    a_lattice = 3.615 * Units.angstrom
    return Material(
        potential=md_driver_module.Morse(rcut, D, alpha, r0),
        atomic_weight=63.546,
        atom_name="Cu",
        a_lattice=a_factor * a_lattice,
    )


def get_Ar_LJ(*, epsilon=None, sigma=None, rcut_fac=3.5):
    """Return argon represented by a Lennard-Jones pair potential."""

    # Standard Lennard-Jones parameters for argon.
    if epsilon is None:
        epsilon = 0.0104 * Units.ev
    if sigma is None:
        sigma = 3.405 * Units.angstrom
    rcut = rcut_fac * sigma
    return Material(
        potential=md_driver_module.LennardJones(rcut, epsilon, sigma),
        atomic_weight=39.948,
        atom_name="Ar",
    )


class MolecularDynamics:
    """Configure and run an MD trajectory.

    Parameters
    ----------
    work_dir : str
        Directory receiving ``series.h5``, ``series.pickle``, trajectory dumps,
        and optional figures.
    **kwargs
        Required historical keys are ``material``, ``ensamble``,
        ``temperature_initial``, ``temperature_thermostat``, ``tau_thermostat``,
        ``pressure_barostat``, ``tau_barostat``, ``dt_global``,
        ``cycles_thermalization``, ``plot_freq``, ``cycles_end``, and
        ``dump_cycles``. Optional controls include neighbor-list settings,
        thermostat chains, plotting, explicit initial velocities, and
        ``random_seed``.

    Notes
    -----
    ``ensamble`` and ``cycles_neigbors`` are intentionally retained misspelled
    because they are part of the published input API.
    """

    def __init__(self, *, work_dir=".", **kwargs):
        np.random.seed(kwargs.get("random_seed", 0))

        self.runtime_begin = time.time()
        self.hostname = socket.gethostname()

        keys_must = {
            "material",
            "temperature_initial",
            "temperature_thermostat",
            "tau_thermostat",
            "pressure_barostat",
            "tau_barostat",
            "ensamble",
            "dump_cycles",
            "cycles_thermalization",
            "plot_freq",
            "cycles_end",
            "dt_global",
        }
        self.work_dir = work_dir
        for key in keys_must:
            assert key in kwargs, f"missing: '{key}'"

        keys_default = {
            "nchains_thermostat",
            "nchains_barostat",
            "fatal_conservation",
            "slurm_plt",
            "show_graphs",
            "plot_graphs",
            "N_lattice_side",
            "a_lattice",
            "neighbors",
            "r_skin",
            "cycles_neigbors",
            "gr_stride",
            "initial_velocities",
            "thermostat_coupling",
            "random_seed",
        }

        assert (
            kwargs.keys() <= keys_must | keys_default
        ), f"extra keys given: {(kwargs.keys() - (keys_must | keys_default))}"

        # -- Default parameters ----
        self.nchains_thermostat = 10
        self.nchains_barostat = 10
        self.slurm_plt = False
        self.show_graphs = False
        self.plot_graphs = True
        self.neighbors = False
        self.r_skin = 1000
        self.cycles_neigbors = 15
        # g(r) sampling stride (in MD steps) *after* thermalization.
        # Default of 1 preserves the old "sample every post-thermalization step" behavior.
        self.gr_stride = 1
        self.thermostat_coupling = MDDriver.THERMOSTAT_KINETIC
        # --------------------------

        for key, value in kwargs.items():
            logger.info(f"Setting parameter '{key}' : {value}")
            if key == "ensamble" and type(value) is str:
                assert value in {"NVE", "NVT", "NPT"}
                value = (
                    MDDriver.NVE
                    if value == "NVE"
                    else MDDriver.NVT
                    if value == "NVT"
                    else MDDriver.NPT
                )
            if key == "thermostat_coupling" and type(value) is str:
                thermostat_coupling_map = {
                    "kinetic": MDDriver.THERMOSTAT_KINETIC,
                    "configurational": MDDriver.THERMOSTAT_CONFIGURATIONAL,
                    "average": MDDriver.THERMOSTAT_AVERAGE,
                }
                assert (
                    value in thermostat_coupling_map
                ), f"thermostat_coupling must be one of {set(thermostat_coupling_map)}"
                value = thermostat_coupling_map[value]
            setattr(self, key, value)

        self.work_dir = path.abspath(self.work_dir)
        logger.info(f"creating output directory: '{self.work_dir}'")
        os.makedirs(self.work_dir, exist_ok=True)
        self.dump_dir = path.join(self.work_dir, "dump")
        os.makedirs(self.dump_dir, exist_ok=True)
        self.pickle_path = path.join(self.work_dir, "series.pickle")
        assert not path.isfile(self.pickle_path), f"'{self.pickle_path}' already exists"

        self.h5_path = path.join(self.work_dir, "series.h5")
        assert not path.isfile(self.h5_path), f"'{self.h5_path}' already exists"

        import matplotlib

        if self.slurm_plt:
            matplotlib.use("Agg")
        from plot import PlotSeries

        self.plot_series = PlotSeries(
            cycles_thermalization=self.cycles_thermalization,
            ensamble="NVE"
            if self.ensamble == MDDriver.NVE
            else "NVT"
            if self.ensamble == MDDriver.NVT
            else "NPT",
            temperature_thermostat=self.temperature_thermostat,
            pressure_barostat=self.pressure_barostat,
            dir_figs=path.join(self.work_dir, "figs"),
            show_graphs=self.show_graphs,
            atom_names=self.material.atom_name.split("-"),
            eam=not (self.material.eam is None and len(self.material.atomic_weight) == 1),
        )

        if self.material.lattice is not None:
            self.initial_positions = self.material.lattice
            self.Lx, self.Ly, self.Lz = self.material.box_length
            self.box = Box(np.diag((self.Lx, self.Ly, self.Lz)))
            self.volume = self.box.volume()
            self.Ntypes = self.material.atom_type
            self.Nparticles = len(self.initial_positions)
            self.a_lattice = np.power(self.volume / self.Nparticles, 1 / 3)
            self.mass = (
                np.array([self.material.atomic_weight[a] for a in self.material.atom_type])
                / Units.Navogadro
            )
            self.xyz_lattice_str = f"a= {self.a_lattice/Units.angstrom} n= {self.Lx/self.a_lattice} {self.Ly/self.a_lattice} {self.Lz/self.a_lattice}"

        elif self.material.a_lattice is not None or "a_lattice" in kwargs.keys():
            if self.material.a_lattice is not None:
                self.a_lattice = self.material.a_lattice
            else:
                self.a_lattice = kwargs["a_lattice"]

            if hasattr(self, "N_lattice_side"):
                self.NsideX, self.NsideY, self.NsideZ = [self.N_lattice_side] * 3
            else:
                assert (
                    hasattr(self, "NsideX") and hasattr(self, "NsideY") and hasattr(self, "NsideZ")
                )

            self.Lx, self.Ly, self.Lz = (
                self.NsideX * self.a_lattice,
                self.NsideY * self.a_lattice,
                self.NsideZ * self.a_lattice,
            )
            self.box = Box(np.diag((self.Lx, self.Ly, self.Lz)))
            self.volume = self.box.volume()
            self.build_fcc_lattice()

            self.Ntypes = [0] * self.Nparticles
            self.mass = (
                np.array([self.material.atomic_weight[0]] * self.Nparticles) / Units.Navogadro
            )

        else:
            assert hasattr(self, "a_lattice")

        self.dump_counter = 0

        self.ndensity = self.Nparticles / self.volume
        self.Nf = 3.0 * self.Nparticles  # 3. * (self.Nparticles - 1.)

        # timesteps
        self.thermalized = False

        # g(r)
        drg = self.a_lattice / 500.0
        Lg = min([self.Lx, self.Ly, self.Lz]) / 2.0
        ng = int(Lg / drg)
        self.rg = np.array([(k + 0.5) * drg for k in range(ng)])
        if getattr(self, "initial_velocities", None) is None:
            self.set_initial_velocities()
        else:
            if len(self.initial_velocities) != self.Nparticles:
                raise ValueError(
                    f"initial_velocities length {len(self.initial_velocities)} "
                    f"does not match Nparticles {self.Nparticles}"
                )

        # Nose-Hoover-chain thermostat
        kT = Units.k_boltz * self.temperature_thermostat
        q_thermostat = np.ones(self.nchains_thermostat) * self.tau_thermostat**2.0 * kT
        q_thermostat[0] *= self.Nf

        sigma_p_eta = np.sqrt(q_thermostat * kT)
        p_eta_initial = np.random.normal(0.0, sigma_p_eta)
        thermostat = Thermostat(self.temperature_thermostat, q_thermostat, p_eta_initial, self.Nf)

        # Barostat
        q_barostat = np.ones(self.nchains_barostat) * self.tau_barostat**2.0 * kT
        sigma_p_eta_baro = np.sqrt(q_barostat * kT)
        p_eta_baro_initial = np.random.normal(0.0, sigma_p_eta_baro)
        thermostat_baro = Thermostat(
            self.temperature_thermostat, q_barostat, p_eta_baro_initial, 1.0
        )

        w_barostat = self.Nf * self.tau_barostat**2.0 * kT
        p_epsilon = np.random.normal(0.0, np.sqrt(w_barostat * kT))
        barostat = Barostat(self.pressure_barostat, w_barostat, p_epsilon, thermostat_baro)

        self.driver = MDDriver(
            self.initial_positions,
            self.initial_velocities,
            PotentialVector(self.material.potential),
            EAMVector([] if self.material.eam is None else self.material.eam),
            self.ensamble,
            self.thermostat_coupling,
            thermostat,
            barostat,
            self.box,
            self.rg,
            self.Nf,
            np.array(self.material.atomic_weight) / Units.Navogadro,
            self.Ntypes,
            self.dt_global,
            self.neighbors,
            self.r_skin,
            self.fatal_conservation,
        )

        self.series = Series()
        self.gr_sum = np.zeros(len(self.driver.gr))
        self.gr_samples = 0
        self.Ntypes_unique = len(set(self.Ntypes))
        if self.Ntypes_unique > 1:
            self.gr_partial_sum = np.zeros(
                (len(self.driver.gr_partial), len(self.driver.gr_partial[0]))
            )
            self.gr_partial_samples = 0
        else:
            self.gr_partial_sum = None
            self.gr_partial_samples = 0

    def build_fcc_lattice(self):
        # Centered FCC basis.

        basis = np.zeros((4, 3))
        basis[0, 0] = 0.0
        basis[0, 1] = 0.0
        basis[0, 2] = 0.0
        basis[1, 0] = 0.5
        basis[1, 1] = 0.5
        basis[1, 2] = 0.0
        basis[2, 0] = 0.5
        basis[2, 1] = 0.0
        basis[2, 2] = 0.5
        basis[3, 0] = 0.0
        basis[3, 1] = 0.5
        basis[3, 2] = 0.5

        Nbasis = 4
        offset = np.ones(3) * 0.25

        self.initial_positions = []
        for i in range(self.NsideX):
            for j in range(self.NsideY):
                for k in range(self.NsideZ):
                    for m in range(Nbasis):
                        self.initial_positions.append(
                            Vector3d(
                                -self.Lx / 2.0 + (i + offset[0] + basis[m, 0]) * self.a_lattice,
                                -self.Ly / 2.0 + (j + offset[1] + basis[m, 1]) * self.a_lattice,
                                -self.Lz / 2.0 + (k + offset[2] + basis[m, 2]) * self.a_lattice,
                            )
                        )

        self.xyz_lattice_str = (
            f"a= {self.a_lattice/Units.angstrom} n= {self.NsideX} {self.NsideX} {self.NsideX}"
        )
        self.Nparticles = len(self.initial_positions)

    def set_initial_velocities(self):
        sigmav = np.sqrt(Units.k_boltz * self.temperature_initial / self.mass)
        vx = np.random.normal(0.0, sigmav)
        vy = np.random.normal(0.0, sigmav)
        vz = np.random.normal(0.0, sigmav)

        vx_tot = np.sum(vx * self.mass) / np.sum(self.mass)
        vy_tot = np.sum(vy * self.mass) / np.sum(self.mass)
        vz_tot = np.sum(vz * self.mass) / np.sum(self.mass)

        vx -= vx_tot
        vy -= vy_tot
        vz -= vz_tot

        ekin = (
            (3.0 * self.Nparticles - 3.0) * Units.k_boltz * self.temperature_initial / 2.0
        )  # Remove three degrees of freedom for NVE momentum conservation
        scale_factor = np.sqrt(ekin / (0.5 * np.sum(self.mass * (vx**2 + vy**2 + vz**2))))

        vx *= scale_factor
        vy *= scale_factor
        vz *= scale_factor

        assert (
            abs(ekin / (0.5 * np.sum(self.mass * (vx**2 + vy**2 + vz**2))) - 1.0) < 1e-8
        ), abs(ekin / (0.5 * np.sum(self.mass * (vx**2 + vy**2 + vz**2))) - 1.0)
        self.initial_velocities = [Vector3d(vxi, vyi, vzi) for vxi, vyi, vzi in zip(vx, vy, vz)]

    def dump(self):
        filename = path.join(self.dump_dir, f"{self.dump_counter}.xyz")
        logger.info(f"Dump to file: '{filename}'")
        str_file = (
            f"{self.Nparticles}\n"
            f"{self.xyz_lattice_str} "
            f"Properties=species:S:1:pos:R:3:vel:R:3 Time={self.driver.time/Units.ps}\n"
            + "\n".join(
                [
                    f'{self.material.atom_name.split("-")[a]} '
                    f"{p.x/Units.angstrom} {p.y/Units.angstrom} {p.z/Units.angstrom} "
                    f"{v.x/(Units.angstrom/Units.fs)} {v.y/(Units.angstrom/Units.fs)} {v.z/(Units.angstrom/Units.fs)}"
                    for p, v, a in zip(self.driver.positions, self.driver.velocities, self.Ntypes)
                ]
            )
        )

        print(str_file, file=open(filename, "w"))
        self.dump_counter += 1

    def perform_thermal_average(self):
        logger.info("performing thermodynamic averages")
        assert self.thermalized and self.cycles_thermalization < len(self.series.time)

        self.series.thermal_results.time = self.driver.time
        self.series.thermal_results.cycles_thermalization = self.cycles_thermalization
        self.series.temperature_kin.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.temperature_conf_inverse.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.total_laplacian.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.total_force2.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.pressure.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.stress_tensor.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.born_stress_tensor.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )

        self.series.epot_atom.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.ekin_atom.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.etot_atom.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.etot_cons.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.etot_cons_err.set_averages(0, cycles_thermalization_cum=0)

        self.series.V_V0.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.L_L0.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        logger.info(f"MEAN L/L0={self.series.L_L0.mean}")

        V0 = self.driver.box_initial.volume()
        Vmean = self.series.V_V0.mean * V0
        particle_density = TimeSeries()
        particle_density.data = self.Nparticles / (np.array(self.series.V_V0.data) * V0)
        particle_density.set_averages(
            self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
        )
        self.series.particle_density = particle_density

        # gr
        if self.gr_samples <= 0:
            logger.warning("No g(r) samples were collected; setting g(r) to zeros")
            self.series.thermal_results.gr = np.zeros_like(self.gr_sum)
        else:
            self.series.thermal_results.gr = self.gr_sum / self.gr_samples
        if self.gr_partial_sum is not None:
            if self.gr_partial_samples <= 0:
                logger.warning(
                    "No partial g(r) samples were collected; setting g_{a,b}(r) to zeros"
                )
                self.series.thermal_results.gr_partial = np.zeros_like(self.gr_partial_sum)
            else:
                self.series.thermal_results.gr_partial = (
                    self.gr_partial_sum / self.gr_partial_samples
                )

        # temperature
        if self.ensamble in {MDDriver.NVT, MDDriver.NPT}:
            self.series.thermal_results.tk_tnh_err = abs(
                self.series.temperature_kin.mean / self.driver.thermostat.temperature - 1.0
            )
            self.series.thermal_results.tc_tnh_err = abs(
                1.0 / self.series.temperature_conf_inverse.mean / self.driver.thermostat.temperature
                - 1.0
            )

        # pressure
        self.series.thermal_results.pressure_ideal = (
            Units.k_boltz * self.series.pressure.mean * self.ndensity
        )
        self.series.thermal_results.pressure_baro_err = abs(
            self.series.pressure.mean / self.driver.barostat.pressure - 1.0
        )

        # Pressure from the radial distribution function
        if self.material.eam is None and len(self.material.atomic_weight) == 1:
            prefactor = -2.0 / 3.0 * np.pi * self.Nparticles * self.ndensity
            wvec = np.array([self.material.potential[0].calculate(r).dVdr * r for r in self.rg])
            pgint = prefactor * integrate.simpson(
                self.rg**2 * wvec * self.series.thermal_results.gr, self.rg
            )
            self.series.thermal_results.pressure_g = (
                self.Nparticles * Units.k_boltz * self.series.temperature_kin.mean + pgint
            ) / self.volume
            self.series.thermal_results.pressure_g_err = abs(
                self.series.thermal_results.pressure_g / self.series.pressure.mean - 1.0
            )

        kTth = Units.k_boltz * self.temperature_thermostat
        # compressibility and thermal expansion
        self.series.thermal_results.bulk_modulus = 0.0
        if self.driver.ensamble == MDDriver.NPT:
            # isothermal compressibility
            compressibility = TimeSeries()
            compressibility.data = (self.series.V_V0.std_cum * V0) ** 2 / (
                kTth * self.series.V_V0.mean_cum * V0
            )
            self.series.compressibility = compressibility
            self.series.thermal_results.compressibility = (self.series.V_V0.std * V0) ** 2 / (
                Vmean * kTth
            )
            logger.info(
                f"isothermal compressibility:   {self.series.thermal_results.compressibility*Units.GPa} {compressibility.data[-1]*Units.GPa}[1/GPa]"
            )
            err = abs(self.series.thermal_results.compressibility / compressibility.data[-1] - 1.0)

            bulk_modulus = TimeSeries()
            bulk_modulus.data = 1.0 / compressibility.data
            self.series.bulk_modulus = bulk_modulus
            self.series.thermal_results.bulk_modulus = (
                1.0 / self.series.thermal_results.compressibility
            )
            logger.info(
                f"bulk_modulus:   {self.series.thermal_results.bulk_modulus/Units.GPa} {bulk_modulus.data[-1]/Units.GPa}[1/GPa]"
            )
            assert (
                abs(self.series.thermal_results.bulk_modulus / bulk_modulus.data[-1] - 1.0) < 1e-4
            )

            # Volumetric thermal expansion (Allen-Tildesley, eq. 2.94)
            dVvec = np.array(self.series.V_V0.data) * V0 - Vmean
            dEvec = (
                np.array(self.series.etot_atom.data) - self.series.etot_atom.mean
            ) * self.Nparticles

            volume_thermal_expansion = TimeSeries()
            volume_thermal_expansion.data = (
                dVvec
                * (dEvec + self.pressure_barostat * dVvec)
                / (Vmean * Units.k_boltz * self.temperature_thermostat**2)
            )
            volume_thermal_expansion.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )
            self.series.volume_thermal_expansion = volume_thermal_expansion
            self.series.thermal_results.volume_thermal_expansion = volume_thermal_expansion.mean
            logger.info(
                f"volume thermal expansion: {self.series.thermal_results.volume_thermal_expansion:g} 1/K"
            )

        #  born stress and elasticity tensors
        if self.driver.ensamble in {MDDriver.NVT, MDDriver.NPT}:
            fluc_stress = self.series.born_stress_tensor
            fluc_stress_tensor = np.array(fluc_stress.get_time_array(), dtype=float, copy=True)
            born_stress_tensor = self.series.born_stress_tensor.get_time_array()
            if self.driver.ensamble == MDDriver.NPT:
                # The Born stress is configurational only. In NPT the ideal-gas pressure
                # N k_B T / V fluctuates with V, so it must stay inside the fluctuating
                # stress: sigma' = sigma^B - (N k_B T / V) delta (Lips & Maass, PRE 97,
                # 053002, Eqs. 12, 15-16). Using sigma^B here together with "+K" biased
                # C11, C12 (C13, C22, C23, C33) by +(2 p_id - p_id^2 / K).
                fluc_stress_tensor[:3, :] -= (
                    np.asarray(self.series.particle_density.data) * kTth
                )[None, :]
            # production-window means of the fluctuating stress
            fluc_stress_mean = fluc_stress_tensor[:, self.cycles_thermalization :].mean(axis=1)

            self.series.stress_fluc = TimeSeriesMatrix(6)
            self.series.born_stress_fluc = TimeSeriesMatrix(6)
            self.series.elasticity_tensor = TimeSeriesMatrix(6)
            self.series.elasticity_ideal_gas = TimeSeriesMatrix(6)

            born_elasticity_tensor = self.series.born_elasticity_tensor.get_time_array()
            self.series.born_elasticity_tensor_11 = TimeSeries()
            self.series.born_elasticity_tensor_12 = TimeSeries()
            self.series.born_elasticity_tensor_44 = TimeSeries()

            self.series.born_elasticity_tensor_11.data = (
                born_elasticity_tensor[0, 0, :]
                + born_elasticity_tensor[1, 1, :]
                + born_elasticity_tensor[2, 2, :]
            ) / 3.0
            self.series.born_elasticity_tensor_12.data = (
                born_elasticity_tensor[0, 1, :]
                + born_elasticity_tensor[0, 2, :]
                + born_elasticity_tensor[1, 2, :]
            ) / 3.0
            self.series.born_elasticity_tensor_44.data = (
                born_elasticity_tensor[3, 3, :]
                + born_elasticity_tensor[4, 4, :]
                + born_elasticity_tensor[5, 5, :]
            ) / 3.0

            self.series.born_elasticity_tensor_11.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )
            self.series.born_elasticity_tensor_12.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )
            self.series.born_elasticity_tensor_44.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )

            self.series.thermal_results.elasticity_tensor = np.zeros((6, 6))

            for i in range(6):
                for j in range(6):
                    # born term
                    self.series.born_elasticity_tensor.data[i][j].data[:] = born_elasticity_tensor[
                        i, j, :
                    ]
                    self.series.born_elasticity_tensor.data[i][j].set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )

                    # stress fluctuations term

                    V_all_times = TimeSeries()
                    V_all_times.data = np.copy(self.series.V_V0.data) * V0
                    V_all_times.set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )
                    self.series.stress_fluc.data[i][j].data[:] = list(
                        V_all_times.mean
                        / kTth
                        * (
                            fluc_stress_tensor[i, :] * fluc_stress_tensor[j, :]
                            - fluc_stress_mean[i] * fluc_stress_mean[j]
                        )
                    )
                    self.series.stress_fluc.data[i][j].set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )

                    self.series.born_stress_fluc.data[i][j].data[:] = list(
                        V_all_times.mean
                        / kTth
                        * (
                            born_stress_tensor[i, :] * born_stress_tensor[j, :]
                            - self.series.born_stress_tensor.data[i].mean
                            * self.series.born_stress_tensor.data[j].mean
                        )
                    )
                    self.series.born_stress_fluc.data[i][j].set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )

                    C_kinetic = np.copy(self.series.particle_density.data) * kTth
                    if i == j:
                        if i < 3:
                            C_kinetic *= 2
                    else:
                        C_kinetic *= 0.0
                    self.series.elasticity_ideal_gas.data[i][j].data[:] = list(C_kinetic)
                    self.series.elasticity_ideal_gas.data[i][j].set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )

                    bulk_modulus_term = 0.0
                    if i < 3 and j < 3 and self.driver.ensamble == MDDriver.NPT:
                        bulk_modulus_term = self.series.thermal_results.bulk_modulus

                    self.series.elasticity_tensor.data[i][j].data[:] = list(
                        np.array(self.series.born_elasticity_tensor.data[i][j].data[:])
                        - np.array(self.series.stress_fluc.data[i][j].data[:])
                        + np.array(self.series.elasticity_ideal_gas.data[i][j].data[:])
                        + bulk_modulus_term
                    )

                    self.series.elasticity_tensor.data[i][j].set_averages(
                        self.cycles_thermalization,
                        cycles_thermalization_cum=self.cycles_thermalization,
                    )
                    self.series.thermal_results.elasticity_tensor[i, j] = (
                        self.series.born_elasticity_tensor.data[i][j].mean
                        - self.series.stress_fluc.data[i][j].mean
                        + self.series.elasticity_ideal_gas.data[i][j].mean
                        + bulk_modulus_term
                    )

                    if i <= j:
                        logger.info(
                            f"C{i+1}{j+1} ELASTIC CONSTANT (cycles {self.cycles_thermalization} {self.driver.cycle})"
                        )
                        logger.info(
                            f"C{i+1}{j+1} BORN TERM                 {self.series.born_elasticity_tensor.data[i][j].mean/Units.GPa}"
                        )
                        logger.info(
                            f"C{i+1}{j+1} STRESS CORRELATION TERM   {-self.series.stress_fluc.data[i][j].mean/Units.GPa}"
                        )
                        logger.info(
                            f"C{i+1}{j+1} IDEAL GAS TERMS           {self.series.elasticity_ideal_gas.data[i][j].mean/Units.GPa}"
                        )
                        logger.info(
                            f"C{i+1}{j+1} Bulk Modulus term         {bulk_modulus_term/Units.GPa}"
                        )
                        logger.info(
                            f"C{i+1}{j+1} TOTAL:                    {self.series.thermal_results.elasticity_tensor[i,j]/Units.GPa}"
                        )
                        logger.info(
                            f"C{i+1}{j+1} TOTAL Mean:               {self.series.elasticity_tensor.data[i][j].mean/Units.GPa}"
                        )
                        logger.info("------------------------------------------")

            # component averages to get C11 C12 C44
            self.series.thermal_results.elasticity_tensor_11 = (
                self.series.thermal_results.elasticity_tensor[0, 0]
                + self.series.thermal_results.elasticity_tensor[1, 1]
                + self.series.thermal_results.elasticity_tensor[2, 2]
            ) / 3.0
            self.series.thermal_results.elasticity_tensor_12 = (
                self.series.thermal_results.elasticity_tensor[0, 1]
                + self.series.thermal_results.elasticity_tensor[0, 2]
                + self.series.thermal_results.elasticity_tensor[1, 2]
            ) / 3.0
            self.series.thermal_results.elasticity_tensor_44 = (
                self.series.thermal_results.elasticity_tensor[3, 3]
                + self.series.thermal_results.elasticity_tensor[4, 4]
                + self.series.thermal_results.elasticity_tensor[5, 5]
            ) / 3.0

            self.series.elasticity_tensor_11 = TimeSeries()
            self.series.elasticity_tensor_12 = TimeSeries()
            self.series.elasticity_tensor_44 = TimeSeries()

            self.series.elasticity_tensor_11.data = (
                np.array(self.series.elasticity_tensor.data[0][0].data[:])
                + np.array(self.series.elasticity_tensor.data[1][1].data[:])
                + np.array(self.series.elasticity_tensor.data[2][2].data[:])
            ) / 3.0
            self.series.elasticity_tensor_12.data = (
                np.array(self.series.elasticity_tensor.data[0][1].data[:])
                + np.array(self.series.elasticity_tensor.data[0][2].data[:])
                + np.array(self.series.elasticity_tensor.data[1][2].data[:])
            ) / 3.0
            self.series.elasticity_tensor_44.data = (
                np.array(self.series.elasticity_tensor.data[3][3].data[:])
                + np.array(self.series.elasticity_tensor.data[4][4].data[:])
                + np.array(self.series.elasticity_tensor.data[5][5].data[:])
            ) / 3.0

            self.series.elasticity_tensor_11.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )
            self.series.elasticity_tensor_12.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )
            self.series.elasticity_tensor_44.set_averages(
                self.cycles_thermalization, cycles_thermalization_cum=self.cycles_thermalization
            )

            bulk_modulus_term = (
                self.series.thermal_results.bulk_modulus
                if self.driver.ensamble == MDDriver.NPT
                else 0.0
            )
            logger.info(f"C11->(C11+C22+C33)/3 Tensor Averages summary")
            logger.info(
                f"C11->(C11+C22+C33)/3 BORN TERM                 {self.series.born_elasticity_tensor_11.mean/Units.GPa}"
            )
            logger.info(
                f"C11->(C11+C22+C33)/3 STRESS CORRELATION TERM   {-(self.series.stress_fluc.data[0][0].mean+self.series.stress_fluc.data[1][1].mean+self.series.stress_fluc.data[2][2].mean)/3./Units.GPa}"
            )
            logger.info(
                f"C11->(C11+C22+C33)/3 KINETIC TERMS             {self.series.elasticity_ideal_gas.data[0][0].mean/Units.GPa}"
            )
            if self.driver.ensamble == MDDriver.NPT:
                logger.info(
                    f"C11->(C11+C22+C33)/3 Bulk Modulus term         {bulk_modulus_term/Units.GPa}"
                )
            logger.info(
                f"C11->(C11+C22+C33)/3 TOTAL:                    {self.series.thermal_results.elasticity_tensor_11/Units.GPa}"
            )
            logger.info("------------------------------------------")

            logger.info(f"C12->(C12+C13+C23)/3 Tensor Averages summary")
            logger.info(
                f"C12->(C12+C13+C23)/3 BORN TERM                 {self.series.born_elasticity_tensor_12.mean/Units.GPa}"
            )
            logger.info(
                f"C12->(C12+C13+C23)/3 STRESS CORRELATION TERM   {-(self.series.stress_fluc.data[0][1].mean+self.series.stress_fluc.data[0][2].mean+self.series.stress_fluc.data[1][2].mean)/3./Units.GPa}"
            )
            logger.info(
                f"C12->(C12+C13+C23)/3 KINETIC TERMS             {self.series.elasticity_ideal_gas.data[0][1].mean/Units.GPa}"
            )
            if self.driver.ensamble == MDDriver.NPT:
                logger.info(
                    f"C12->(C12+C13+C23)/3 Bulk Modulus term         {bulk_modulus_term/Units.GPa}"
                )
            logger.info(
                f"C12->(C12+C13+C23)/3 TOTAL:                    {self.series.thermal_results.elasticity_tensor_12/Units.GPa}"
            )
            logger.info("------------------------------------------")

            logger.info(f"C44->(C44+C55+C66)/3 Tensor Averages summary")
            logger.info(
                f"C44->(C44+C55+C66)/3 BORN TERM                 {self.series.born_elasticity_tensor_44.mean/Units.GPa}"
            )
            logger.info(
                f"C44->(C44+C55+C66)/3 STRESS CORRELATION TERM   {-(self.series.stress_fluc.data[3][3].mean+self.series.stress_fluc.data[4][4].mean+self.series.stress_fluc.data[5][5].mean)/3./Units.GPa}"
            )
            logger.info(
                f"C44->(C44+C55+C66)/3 KINETIC TERMS             {self.series.elasticity_ideal_gas.data[3][3].mean/Units.GPa}"
            )
            logger.info(
                f"C44->(C44+C55+C66)/3 TOTAL:                    {self.series.thermal_results.elasticity_tensor_44/Units.GPa}"
            )
            logger.info("------------------------------------------")

        # etot/epot/ekin

        # Energy from the radial distribution function
        if self.material.eam is None and len(self.material.atomic_weight) == 1:
            prefactor = 2.0 * np.pi * self.ndensity
            vvec = np.array([self.material.potential[0].calculate(r).V for r in self.rg])
            egint = prefactor * integrate.simpson(
                self.rg**2 * vvec * self.series.thermal_results.gr, self.rg
            )
            self.series.thermal_results.etot_atom_g = (
                1.5 * Units.k_boltz * self.series.temperature_kin.mean + egint
            )
            self.series.thermal_results.etot_atom_g_err = abs(
                self.series.thermal_results.etot_atom_g / self.series.etot_atom.mean - 1.0
            )

            logger.info(
                f"eg={self.series.thermal_results.etot_atom_g:g}, e={self.series.etot_atom.mean:g}, err={self.series.thermal_results.etot_atom_g_err:g}"
            )

        # energy conservation
        if self.driver.ensamble == MDDriver.NVE:
            self.series.thermal_results.epot_ekin_std_err = abs(
                self.series.epot_atom.std / self.series.ekin_atom.std - 1.0
            )
            self.series.thermal_results.epot_ekin_slope_fit = np.polyfit(
                self.series.ekin_atom, self.series.epot_atom, deg=1.0
            )[0]

        # Cv
        self.series.thermal_results.cv_id = 1.5 * Units.k_boltz

        if self.driver.ensamble == MDDriver.NVE:
            get_cv_nve = (
                lambda sigma, tmp: 1.5
                * Units.k_boltz
                / (1.0 - sigma**2.0 / (1.5 * self.Nparticles * (Units.k_boltz * tmp) ** 2))
            )
            self.series.thermal_results.cv_ekin = get_cv_nve(
                self.series.ekin_atom.std * self.Nparticles, self.series.temperature_kin.mean
            )
            self.series.thermal_results.cv_epot = get_cv_nve(
                self.series.epot_atom.std * self.Nparticles, self.series.temperature_kin.mean
            )
            self.series.thermal_results.cv_epot_ekin_err = abs(
                self.series.thermal_results.cv_epot / self.series.thermal_results.cv_ekin - 1.0
            )

        elif self.driver.ensamble == MDDriver.NVT:
            self.series.thermal_results.ekin_atom_std_exact = (
                np.sqrt(
                    1.5 * self.Nparticles * (Units.k_boltz * self.series.temperature_kin.mean) ** 2
                )
                / self.Nparticles
            )
            self.series.thermal_results.ekin_atom_std_err = abs(
                self.series.ekin_atom.std / self.series.thermal_results.ekin_atom_std_exact - 1.0
            )

            self.series.thermal_results.cv_etot = (
                (self.series.etot_atom.std * self.Nparticles) ** 2
                / (Units.k_boltz * self.series.temperature_kin.mean**2)
                / self.Nparticles
            )
            self.series.thermal_results.cv_epot = (
                self.series.epot_atom.std * self.Nparticles
            ) ** 2 / (
                Units.k_boltz * self.series.temperature_kin.mean**2
            ) / self.Nparticles + 1.5 * Units.k_boltz
            self.series.thermal_results.cv_etot_epot_err = abs(
                self.series.thermal_results.cv_epot / self.series.thermal_results.cv_etot - 1.0
            )
            logger.info(
                f"cv_id={self.series.thermal_results.cv_id:g} cv_epot={self.series.thermal_results.cv_epot:g} cv_id={self.series.thermal_results.cv_id:g}"
            )
            logger.info(
                f"dK {self.series.ekin_atom.std} {self.series.thermal_results.ekin_atom_std_exact}"
            )
        elif self.driver.ensamble == MDDriver.NPT:
            pass
        else:
            logger.fatal("bad ensamble")
            quit()

        self.series.thermal_results.rg = self.rg

    def save_pickle(self):
        logger.info(f"Pickle to file: '{self.pickle_path }'")

        pickle_full = self.pickle_path.replace(".pickle", ".pickle.full")
        with open(pickle_full, "wb") as file:
            pickle.dump(self.series, file)

    def save_h5(self):
        logger.info(f"Writing h5 file: '{self.h5_path }'")
        with h5py.File(self.h5_path, "w") as fh5:
            for key, value in self.series.thermal_results.__dict__.items():
                fh5[f"thermal_results/{key}"] = value

            fh5[f"series/time"] = self.series.time

            for key, value in self.series.__dict__.items():
                if type(value) in {TimeSeries, TimeSeriesVector, TimeSeriesMatrix}:
                    logger.info(f"writing {key!r} {type(value)}")
                    value.write_h5(file_h5=fh5, group_name=f"series/{key}")

        logger.info(f"Finished writing to h5 file: {self.h5_path!r}")

    def plot(self):
        assert self.cycles_thermalization is not None
        assert (
            self.driver.time == self.series.thermal_results.time
        )  # make sure we are plotting at the same time the averages were calculated
        if not self.plot_graphs:
            if self.driver.ensamble == MDDriver.NVE:
                logger.info("plotting NVE temperature/pressure graphs...")
                self.plot_series.temperature(series=self.series)
                self.plot_series.pressure(series=self.series)
            return

        logger.info("plotting thermodynamic graphs...")
        self.plot_series.temperature(series=self.series)
        self.plot_series.pressure(series=self.series)
        self.plot_series.elasticity(series=self.series)
        self.plot_series.volume_fluctuations(series=self.series)
        self.plot_series.energies_cv(series=self.series)
        self.plot_series.partial_pair_correlation(series=self.series)
        self.plot_series.conservation(series=self.series)

    def run(self):
        while self.driver.cycle <= self.cycles_end:
            if self.driver.cycle % self.dump_cycles == 0:
                self.dump()

            if self.driver.cycle % self.cycles_neigbors == 0 and self.neighbors:
                self.driver.calculate_neighbors()

            self.driver.integrate()

            # ---- store results for averages calculation
            self.series.time.append(self.driver.time)
            self.series.temperature_kin.append(self.driver.temperature_kin)
            self.series.temperature_conf_inverse.append(self.driver.temperature_conf_inverse)
            self.series.total_laplacian.append(self.driver.total_laplacian)
            self.series.total_force2.append(self.driver.total_force2)
            self.series.pressure.append(self.driver.pressure)
            self.series.stress_tensor.append(np.copy(self.driver.stress_tensor))
            self.series.etot_atom.append(self.driver.etot_atom)
            self.series.epot_atom.append(self.driver.epot_atom)
            self.series.epot_eam_atom.append(self.driver.epot_eam_atom)
            self.series.ekin_atom.append(self.driver.ekin_atom)
            self.series.etot_cons.append(self.driver.etot_cons)
            self.series.etot_cons_err.append(self.driver.etot_cons_err)
            self.series.born_stress_tensor.append(np.copy(self.driver.born_stress_tensor))
            self.series.born_elasticity_tensor.append(np.copy(self.driver.born_elasticity_tensor))

            V_V0 = self.driver.box.volume() / self.driver.box_initial.volume()
            self.series.V_V0.append(V_V0)
            self.series.L_L0.append(V_V0 ** (1.0 / 3.0))
            self.thermalized = self.driver.cycle > self.cycles_thermalization

            if self.thermalized:
                if self.gr_stride < 1:
                    raise ValueError(f"gr_stride must be >= 1, got {self.gr_stride}")
                if (self.driver.cycle - self.cycles_thermalization) % self.gr_stride == 0:
                    # g(r) is expensive; compute it explicitly only when we actually sample it.
                    # (The C++ driver no longer computes it automatically each step.)
                    if self.gr_partial_sum is not None:
                        self.driver.calculate_pair_correlations()
                        self.gr_sum += np.array(self.driver.gr)
                        self.gr_samples += 1
                        self.gr_partial_sum += np.array(self.driver.gr_partial)
                        self.gr_partial_samples += 1
                    else:
                        self.driver.calculate_pair_correlation()
                        self.gr_sum += np.array(self.driver.gr)
                        self.gr_samples += 1

            if (
                self.driver.cycle - self.cycles_thermalization
            ) % self.plot_freq == 0 and self.thermalized:
                self.perform_thermal_average()
                self.save_h5()
                self.plot()

        self.runtime = time.time() - self.runtime_begin
        logger.info(f"hostname: {self.hostname}")
        logger.info(f"rutime was: {self.runtime}")

        return self
