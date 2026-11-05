# Clement-Auger Cu-Zn potential notice

`potentials/ZnCu.eam.alloy` contains the EAM potential developed by Antoine
Clement and Thierry Auger:

A. Clement and T. Auger, "An EAM potential for alpha-brass copper-zinc alloys:
application to plasticity and fracture," *Modelling and Simulation in Materials
Science and Engineering* 31 (2023) 015004,
<https://doi.org/10.1088/1361-651X/aca4ec>.

Original distribution:

- NIST Interatomic Potentials Repository:
  <https://www.ctcms.nist.gov/potentials/entry/2022--Clement-A-Auger-T--Cu-Zn/>
- OpenKIM model `EAM_Dynamo_ClementAuger_2022_CuZn__MO_875501152989_000`:
  <https://doi.org/10.25950/79dabed3>

The OpenKIM distribution identifies this model as LGPL-2.1 licensed. A copy of
that license is included as `Clement-Auger-LGPL-2.1.txt`.

The copy in this package is numerically identical to the NIST/OpenKIM
`CuZn.eam.alloy` data from its `2 Zn Cu` SETFL header onward. The deposited
mdelast version was modified before 15 July 2026 by converting line endings,
renaming the file to `ZnCu.eam.alloy`, removing leading blank lines, and adding
two provenance lines. Its SHA-256 checksum is
`bdb0497544fa28eff8ed81fc799da9bf8d6cf8b4e907658f40b18c13652ee2ce`.

The mdelast MIT license does not relicense this third-party potential.
