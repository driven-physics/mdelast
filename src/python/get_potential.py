"""Readers for tabulated embedded-atom method (EAM) potential files.

Potential files are resolved relative to the package's top-level
``potentials`` directory. Arrays are converted to the CGS convention used by
the compiled driver before pybind11 potential objects are constructed.
"""

import numpy as np
from os import path
import sys
from scipy.interpolate import UnivariateSpline
from molecular_dynamics import Material

try:
    import md_driver_module
    from md_driver_module import Units
except:
    pass

import logging

if len(logging.root.handlers) == 0:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("MD")


def get_deriv(x, y):
    """Return midpoint coordinates and first finite-difference derivatives."""

    return 0.5 * (x[1:] + x[:-1]), np.diff(y) / np.diff(x)


def get_derivatives(x, y):
    """Evaluate first and second derivatives of an exact quintic spline."""

    spl = UnivariateSpline(x=x, y=y, k=5, s=0)
    dydx = spl.derivative(n=1)(x)
    d2ydx2 = spl.derivative(n=2)(x)
    return dydx, d2ydx2


def get_EAM_alloy_file(file_name):
    """Load a multi-element LAMMPS SETFL alloy potential.

    Returns
    -------
    elements : list
        Per-element records ``[symbol, atomic_number, atomic_weight,
        lattice_constant_angstrom]``.
    potentials : list
        Pair potentials in the triangular ordering expected by ``MDDriver``.
    eams : list
        One embedding/density object per element.
    """

    pot_dir = path.join(path.dirname(path.abspath(__file__)), "../../potentials")
    file_path = path.join(pot_dir, file_name)
    logger.info(f"Reading file potential... '{file_path}'")

    with open(file_path, "r") as f:
        lines = f.readlines()

    header_lines = 3

    element_info = lines[header_lines].split()

    num_elements = int(element_info[0])
    elements = [[element_info[i]] for i in range(1, num_elements + 1)]

    param_line = lines[header_lines + 1].split()
    Nrho = int(param_line[0])
    drho = float(param_line[1]) / Units.angstrom**3
    Nr = int(param_line[2])
    dr = float(param_line[3]) * Units.angstrom
    rcut = float(param_line[4]) * Units.angstrom

    r_grid = np.linspace(0.0, dr * (Nr - 1.0), Nr)  # angstrom
    rho_grid = np.linspace(0.0, drho * (Nrho - 1.0), Nrho)

    index = 5
    eams = []
    for i in range(num_elements):
        atomic_number, atomic_weight, a_lattice = [
            float(s) for s in lines[index].split("\n")[0].split()[:-1]
        ]
        elements[i].extend([atomic_number, atomic_weight, a_lattice])

        a_lattice *= Units.angstrom
        logger.info(
            f"{elements[i][0]}: Z={atomic_number:g}, ATW={atomic_weight:g}, a={a_lattice/Units.angstrom:g}[A]"
        )
        index += 1
        data = np.loadtxt(file_path, skiprows=index, max_rows=Nrho + Nr)
        F = np.array(data[:Nrho]) * Units.ev  # embedding function F(rho) (Nrho values)
        rho = np.array(data[Nrho : Nrho + Nr]) / Units.angstrom**3
        drhodr, d2rhodr2 = get_derivatives(x=r_grid, y=rho)
        dFdrho, d2Fdrho = get_derivatives(x=rho_grid, y=F)
        eams.append(
            md_driver_module.EAM(
                r_grid[-1], rho_grid[-1], dr, drho, rho, drhodr, d2rhodr2, F, dFdrho, d2Fdrho
            )
        )
        index += Nrho + Nr
    potentials = [0] * int(num_elements * (num_elements + 1) / 2)
    i, j = 0, 0
    for _ in range(int(num_elements * (num_elements + 1) / 2)):
        pot_index = int(j * (j + 1) / 2 + i)
        phir = np.loadtxt(file_path, skiprows=index, max_rows=Nr)
        phi = phir / r_grid * Units.ev * Units.angstrom
        phi[0] = phi[1]

        rmin_fatal = 1e-2 * Units.angstrom
        dphidr, d2phidr2 = get_derivatives(x=r_grid, y=phi)

        potentials[pot_index] = md_driver_module.TablePotential(
            r_grid[-1], rmin_fatal, dr, phi, dphidr, d2phidr2
        )

        index += Nr
        # setfl stores the pair blocks as 11, 21, 22, 31, 32, 33, ...
        # i.e. (i, j) = (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ... with i <= j,
        # matching pot_index = j*(j+1)/2 + i used by MDDriver.
        if i == j:
            j += 1
            i = 0
        else:
            i += 1

    return elements, potentials, eams


def get_EAM_file(*, file_name, a_factor=1.0):
    """Load a supported single-element EAM file as a ``Material``."""

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


def get_Ar_LJ(*, epsilon=None, sigma=None, rcut_fac=3.5):
    """Return an argon Lennard-Jones ``Material``."""

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
