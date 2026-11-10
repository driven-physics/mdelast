# Purja Pun-Yamakov-Mishin Ni-Al-Co potential notice

`potentials/Mishin_updated-Ni-Al-Co-2013.eam.alloy` contains the EAM potential for
the Ni-Al-Co system developed by Ganga P. Purja Pun, Vesselin Yamakov, and Yuri
Mishin:

G. P. Purja Pun, V. Yamakov, and Y. Mishin, "Interatomic potential for the
ternary Ni-Al-Co system and application to atomistic modeling of the B2-L1_0
martensitic transformation," *Modelling and Simulation in Materials Science and
Engineering* 23 (2015) 065006, <https://doi.org/10.1088/0965-0393/23/6/065006>.

Original distributions:

- NIST Interatomic Potentials Repository, file
  `Mishin_updated-Ni-Al-Co-2013.eam.alloy` (ipr2), sent by G. Purja Pun on
  12 October 2015 and posted on 15 December 2015; it corrects the cutoff
  distance of the Co interactions of the 2013 file:
  <https://www.ctcms.nist.gov/potentials/entry/2015--Purja-Pun-G-P-Yamakov-V-Mishin-Y--Ni-Al-Co/>
- OpenKIM model
  `EAM_Dynamo_PunYamakovMishin_2013_NiAlCo__MO_826591359508_001`, based on the
  NIST ipr2 file: <https://doi.org/10.25950/da2eab38>

The OpenKIM distribution states in its LICENSE file: "This work is dedicated to
the public domain." mdelast's MIT license does not apply to this file.

The file in this package is byte-for-byte identical to OpenKIM's
`Mishin_updated-Ni-Al-Co-2013.eam.alloy` (SHA-1
`c0851f77d0d7f411ddfb7b6d838c9df1eede6bea`, as listed in the model's
kimprovenance.edn); it has CRLF line endings. Apart from the line endings and a
final newline it is identical to the NIST ipr2 file. Its SHA-256 checksum is
`3118f1b0423210597c20450e79372d12f4b8a303a1280c2fd1191a3f148363f6`.

It is used by Part C of benchmark S (`benchmarks/smoke/`). Please cite the original
paper when you use the potential.
