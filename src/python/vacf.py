"""
Velocity autocorrelation (VACF) and vibrational density of states (VDOS).

Mass-weighted total and per-species partial VACFs follow the normalization used
for binary alloys (total + Cu + Zn partials).
"""

from __future__ import annotations

import os
import re
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

try:
    from md_driver_module import Units, Vector3d
except ImportError:
    Units = None
    Vector3d = None

# THz -> cm^-1 (conventional spectroscopy factor)
THZ_TO_CM1 = 33.356411

FIGURE_RC_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "figure.figsize": (8.0, 5.0),
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.0,
    "lines.markersize": 4,
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "legend.edgecolor": "0.8",
}

VACF_COLORS = {
    "total": "#1f77b4",
    "Cu": "#ff7f0e",
    "Zn": "#2ca02c",
}


def _require_units():
    if Units is None or Vector3d is None:
        raise ImportError("md_driver_module is not available; add the build directory to sys.path")


def read_dump_xyz(
    filepath: str,
    *,
    element_order: Optional[Mapping[str, int]] = None,
) -> Dict:
    """
    Read an extended XYZ dump written by MolecularDynamics.dump().

    Parameters
    ----------
    filepath
        Path to .xyz file.
    element_order
        Map species name -> integer type index (e.g. {'Zn': 0, 'Cu': 1}).
        If None, species are assigned 0..K-1 in order of first appearance.

    Returns
    -------
    dict with keys: species, positions, velocities, box_lengths, atom_types,
    time_ps, n_atoms, a_angstrom, n_lattice
    """
    _require_units()
    vel_scale = Units.angstrom / Units.fs

    with open(filepath, "r", encoding="utf-8") as f:
        n_atoms = int(f.readline().strip())
        header = f.readline().strip()

        m_a = re.search(r"a=\s*([\d.eE+-]+)", header)
        m_n = re.search(r"n=\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", header)
        m_t = re.search(r"Time=([\d.eE+-]+)", header)
        if m_a is None or m_n is None:
            raise ValueError(f"Could not parse lattice header in {filepath!r}: {header!r}")

        a_angstrom = float(m_a.group(1))
        n_lattice = tuple(float(x) for x in m_n.groups())
        time_ps = float(m_t.group(1)) if m_t else 0.0

        species: list[str] = []
        positions = []
        velocities = []

        for _ in range(n_atoms):
            parts = f.readline().split()
            if len(parts) < 7:
                raise ValueError(f"Expected 7 fields per atom line in {filepath!r}")
            sp, x, y, z, vx, vy, vz = parts[0], *map(float, parts[1:7])
            species.append(sp)
            positions.append(
                Vector3d(
                    x * Units.angstrom,
                    y * Units.angstrom,
                    z * Units.angstrom,
                )
            )
            velocities.append(
                Vector3d(
                    vx * vel_scale,
                    vy * vel_scale,
                    vz * vel_scale,
                )
            )

    if element_order is None:
        unique = []
        for sp in species:
            if sp not in unique:
                unique.append(sp)
        element_order = {sp: i for i, sp in enumerate(unique)}

    unknown = sorted(set(species) - set(element_order))
    if unknown:
        raise ValueError(f"Species {unknown} not in element_order {list(element_order)}")

    atom_types = [element_order[sp] for sp in species]
    Lx = a_angstrom * n_lattice[0] * Units.angstrom
    Ly = a_angstrom * n_lattice[1] * Units.angstrom
    Lz = a_angstrom * n_lattice[2] * Units.angstrom

    return {
        "species": species,
        "positions": positions,
        "velocities": velocities,
        "box_lengths": (Lx, Ly, Lz),
        "atom_types": atom_types,
        "time_ps": time_ps,
        "n_atoms": n_atoms,
        "a_angstrom": a_angstrom,
        "n_lattice": n_lattice,
    }


def _mass_weighted_vacf_segment(
    vel: np.ndarray,
    masses: np.ndarray,
    mask: Optional[np.ndarray],
    max_lag: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    vel: (T, N, 3), masses: (N,), mask: (N,) bool or None.
    Returns time_lags (max_lag+1,) in snapshot index units, normalized vacf.
    """
    if mask is not None:
        vel = vel[:, mask, :]
        masses = masses[mask]

    if vel.shape[1] == 0:
        raise ValueError("No atoms in species mask")

    n_snap = vel.shape[0]
    max_lag = min(max_lag, n_snap - 1)
    numer = np.zeros(max_lag + 1, dtype=np.float64)
    denom = 0.0

    for lag in range(max_lag + 1):
        n_orig = n_snap - lag
        if n_orig <= 0:
            break
        v0 = vel[:n_orig]
        vt = vel[lag : lag + n_orig]
        dot = np.sum(v0 * vt, axis=2)
        weighted = np.sum(dot * masses[np.newaxis, :], axis=1)
        numer[lag] = np.mean(weighted)
        if lag == 0:
            denom = numer[0]

    if denom <= 0.0:
        raise ValueError("VACF denominator is zero (no kinetic energy in sample)")

    vacf = numer / denom
    lags = np.arange(max_lag + 1, dtype=np.float64)
    return lags, vacf


def compute_vacf(
    velocities: np.ndarray,
    masses: np.ndarray,
    species: Sequence[str],
    *,
    dt_ps: float,
    max_correlation_time_ps: float,
    species_names: Iterable[str] = ("Cu", "Zn"),
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Mass-weighted VACF (total + partial per species).

    Parameters
    ----------
    velocities
        (n_snapshots, n_atoms, 3) in consistent velocity units (e.g. CGS cm/s).
    masses
        Per-atom masses, same length as n_atoms.
    species
        Length n_atoms, element symbol per atom ('Cu', 'Zn', ...).
    dt_ps
        Time between velocity snapshots [ps].
    max_correlation_time_ps
        Maximum lag time for the correlation [ps].

    Returns
    -------
    dict label -> (time_ps, C_v(t)/C_v(0))
    """
    if velocities.ndim != 3 or velocities.shape[2] != 3:
        raise ValueError("velocities must have shape (n_snapshots, n_atoms, 3)")

    masses = np.asarray(masses, dtype=np.float64)
    species_arr = np.asarray(species)
    if masses.shape[0] != velocities.shape[1]:
        raise ValueError("masses length must match n_atoms")
    if species_arr.shape[0] != velocities.shape[1]:
        raise ValueError("species length must match n_atoms")

    max_lag = int(round(max_correlation_time_ps / dt_ps))
    max_lag = max(0, min(max_lag, velocities.shape[0] - 1))

    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    lags_idx, vacf_total = _mass_weighted_vacf_segment(velocities, masses, None, max_lag)
    time_ps = lags_idx * dt_ps
    out["total"] = (time_ps, vacf_total)

    for name in species_names:
        mask = species_arr == name
        if not np.any(mask):
            continue
        _, vacf_part = _mass_weighted_vacf_segment(velocities, masses, mask, max_lag)
        out[name] = (time_ps, vacf_part)

    return out


def compute_vdos_moments(
    freq_thz: np.ndarray,
    g: np.ndarray,
    *,
    min_freq_thz: float = 0.0,
    max_freq_thz: float = 10.0,
) -> Dict[str, float]:
    """
    Spectral moments from VDOS: <omega^n> and <omega^4>/<omega^2>^2.

    Moments use trapezoidal integration with g normalized to unit area.
    freq_thz is in THz; reported <omega^2> in (THz)^2, etc.
    Integration is restricted to [min_freq_thz, max_freq_thz].
    """
    freq = np.asarray(freq_thz, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    mask = (freq >= min_freq_thz) & (freq <= max_freq_thz)
    freq = freq[mask]
    g = np.maximum(g[mask], 0.0)
    if freq.size < 2:
        raise ValueError("Need at least two frequency points for moments")

    norm = np.trapz(g, freq)
    if norm <= 0.0:
        raise ValueError("VDOS integral is zero")

    g_norm = g / norm
    m2 = float(np.trapz(g_norm * freq**2, freq))
    m4 = float(np.trapz(g_norm * freq**4, freq))
    if m2 <= 0.0:
        raise ValueError("<omega^2> is zero")

    return {
        "omega2": m2,
        "omega4": m4,
        "omega4_over_omega2_sq": m4 / (m2 * m2),
    }


def compute_vdos(
    time_ps: np.ndarray,
    vacf: np.ndarray,
    *,
    apply_hann: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fourier transform of the VACF -> vibrational density of states.

    Returns
    -------
    freq_thz, g_omega (arbitrary units, nonnegative)
    """
    if len(time_ps) < 4:
        raise ValueError("Need at least 4 VACF points for VDOS FFT")

    dt_ps = float(np.mean(np.diff(time_ps)))
    if not np.allclose(np.diff(time_ps), dt_ps, rtol=1e-6, atol=1e-12):
        raise ValueError("time_ps must be uniformly spaced")

    c = np.asarray(vacf, dtype=np.float64).copy()
    if apply_hann:
        window = np.hanning(len(c))
        c *= window

    spectrum = np.fft.rfft(c)
    g = np.real(spectrum)
    g = np.maximum(g, 0.0)
    freq_thz = np.fft.rfftfreq(len(c), d=dt_ps)
    return freq_thz, g


def _apply_figure_style():
    mpl.rcParams.update(FIGURE_RC_PARAMS)


def _save_figure(fig: plt.Figure, out_dir: str, basename: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600, bbox_inches="tight")


def plot_vacf(
    vacf_results: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    out_dir: str,
    *,
    basename: str = "vacf",
    max_time_ps: float = 10.0,
    show: bool = False,
) -> plt.Figure:
    """Plot normalized VACF vs time [ps]."""
    _apply_figure_style()
    fig, ax = plt.subplots()

    for label, (t_ps, cv) in vacf_results.items():
        color = VACF_COLORS.get(label, None)
        mask = t_ps <= max_time_ps
        ax.plot(t_ps[mask], cv[mask], label=label, color=color)

    ax.set_xlabel(r"$t$ [ps]")
    ax.set_ylabel(r"$C_v(t)\,/\,C_v(0)$")
    ax.set_xlim(0, max_time_ps)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(bottom=max(0, ylo) if ylo > -0.5 else ylo)

    _save_figure(fig, out_dir, basename)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_vdos(
    vdos_results: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    out_dir: str,
    *,
    basename: str = "vdos",
    max_freq_thz: float = 10.0,
    show: bool = False,
) -> plt.Figure:
    """Plot VDOS vs frequency [THz] with cm^-1 top axis."""
    _apply_figure_style()
    fig, ax = plt.subplots()

    for label, (freq_thz, g) in vdos_results.items():
        color = VACF_COLORS.get(label, None)
        mask = freq_thz <= max_freq_thz
        ax.plot(freq_thz[mask], g[mask], label=label, color=color)

    ax.set_xlabel(r"$\nu$ [THz]")
    ax.set_ylabel(r"$g(\nu)$ [a.u.]")
    ax.set_xlim(0, max_freq_thz)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(bottom=max(0, ylo))

    ax2 = ax.secondary_xaxis("top", functions=(lambda f: f * THZ_TO_CM1, lambda w: w / THZ_TO_CM1))
    ax2.set_xlabel(r"$\tilde{\nu}$ [cm$^{-1}$]")

    _save_figure(fig, out_dir, basename)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def write_spectral_moments(
    path: str,
    moments: Mapping[str, Mapping[str, float]],
) -> None:
    """Write spectral moments table (human-readable)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("# VDOS spectral moments (frequency in THz)\n")
        f.write("# columns: <omega^2> [THz^2], <omega^4> [THz^4], <omega^4>/<omega^2>^2\n")
        for label, m in moments.items():
            f.write(
                f"{label}\t{m['omega2']:.8e}\t{m['omega4']:.8e}\t"
                f"{m['omega4_over_omega2_sq']:.8e}\n"
            )


def save_results_npz(
    path: str,
    *,
    vacf_results: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    vdos_results: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    metadata: Optional[dict] = None,
) -> None:
    """Save VACF/VDOS curves to a compressed .npz file."""
    payload = {}
    for label, (t, c) in vacf_results.items():
        payload[f"vacf_time_ps_{label}"] = t
        payload[f"vacf_{label}"] = c
    for label, (f, g) in vdos_results.items():
        payload[f"vdos_freq_thz_{label}"] = f
        payload[f"vdos_{label}"] = g
    if metadata:
        for k, v in metadata.items():
            payload[f"meta_{k}"] = np.array(v)
    np.savez_compressed(path, **payload)
