---
title: "Available MD Software"
teaching: 5
exercises: 0
questions:
- "What software is available on CC systems"
objectives:
- "List MD packages"
---

# Molecular Dynamics Software Available on Compute Canada Systems
## [AMBER](http://ambermd.org)

## [DL_POLY](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx)

## [GROMACS](http://gromacs.org)
### Force fields implemented in GROMACS:

- AMBER: 94, 96, 99, 99sb, 99sb-ildn, 03, GS (amberGS is amber94 with both backbone torsion potentials set to 0).
- CHARMM: 27 (optimized for proteins and nucleic acids).
- [CHARMM for GROMACS](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (CHARMM36, CGenFF)
- GROMOS: 43a1, 43a2, 45a3, 53a5, 53a6, 54a7.
- OPLS-AA (OPLS-AA implemented in GROMACS is actually OPLS-AA/L. OPLS-AA/L uses OPLS-AA atom types with the torsions and impropers refitted to QM calculations at the HF/6-31G** level followed by single-point LMP2/cc-pVTZ(-f))

## [NAMD](https://www.ks.uiuc.edu/Research/namd/)

### Force fields implemented in NAMD:
- [AMBER](http://ambermd.org/AmberModels.php) (amber format topology prepared with AMBERTOOLS). Currently AMBER recommends the following force fields: ff14SB for proteins, OL15 for DNA, OL3 for RNA, GLYCAM_06j for carbohydrates, lipid17 for lipids, and a general force field gaff2.
- [CHARMM](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (charmm format, topology can be prepared with the included **psfgen** program)
- [CHARMM-Drude](http://mackerell.umaryland.edu/charmm_drude_ff.shtml) (the polarizable force field based on the classical Drude oscillator model, charmm format)
- GROMOS (from GROMACS distribution, topology prepared with pdb2gmx)
- [OPLS-AA/M](http://zarbi.chem.yale.edu/oplsaam.html) New peptide dihedral parameters, significantly outperform the previous versions for proteins (charmm format)

## [LAMPPS](https://lammps.sandia.gov)
