# Reference data of benchmark S, Part C

LAMMPS results for one fixed, distorted 500-atom fcc Ni-Al-Co cell (5 x 5 x 5
cubic cells, lattice constant 3.54 A, random occupation with Ni 396, Al 59,
Co 45, Gaussian displacements of 0.05 A, seed 7), with the potential
`Mishin_updated-Ni-Al-Co-2013.eam.alloy`.

Files:

- `NiAlCo.data` — the cell (LAMMPS atomic data file, positions in A).
- `in_NiAlCo.lmp` — the LAMMPS input that produced the reference values. To
  rerun it, copy the potential file from `potentials/` next to it.
- `NiAlCo_forces.dump` — the force on every atom (eV/A).
- `lammps_values.json` — energy per atom, virial (`compute pressure NULL
  virial`, bar), and the 21-component Born matrix (`compute born/matrix numdiff
  1e-6`, energy units), as printed by LAMMPS.

The values were computed with LAMMPS (22 Jul 2025 - Update 4) on one MPI task
by R. Shikler, in the check that verified the pair-block fix of `FIXES.md`,
entry 2. That run used the NIST ipr2 copy of the potential, which holds the
same numbers as the file in `potentials/`. Running benchmark S does not need
LAMMPS.
