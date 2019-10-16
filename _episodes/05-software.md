---
title: "Introducing molecular dynamics simulation"
teaching: 30
exercises: 0
questions:
- "What is Molecular Dynamics and how can I benefit from using it"
- "What is a force field"
objectives:
- "Explain strengths and weaknesses of MD simulations"
- "Explain the interactions between particles in a molecular dynamics simulation"
---

<!-- MDTOC maxdepth:6 firsth1:1 numbering:1 flatten:0 bullets:1 updateOnSave:1 -->

   - 0.1. [Molecular Dynamics Software Available on Compute Canada Systems](#molecular-dynamics-software-available-on-compute-canada-systems)   
      - 0.1.1. [AMBER](#amber)   
      - 0.1.2. [DL_POLY](#dl_poly)   
      - 0.1.3. [GROMACS](#gromacs)   
         - 0.1.3.1. [Force fields implemented in GROMACS:](#force-fields-implemented-in-gromacs)   
      - 0.1.4. [NAMD](#namd)   
         - 0.1.4.1. [Force fields implemented in NAMD:](#force-fields-implemented-in-namd)   
      - 0.1.5. [LAMMPS](#lammps)   
   - 0.2. [Water Models](#water-models)   

<!-- /MDTOC -->

[Link to the GitLab Repository](https://git.computecanada.ca/svassili/bst-md-theory-lesson-novice/blob/gh-pages/_episodes/01-introduction.md)

[Link to the GitHub Repository](https://github.com/ssvassiliev/bst-md-theory-lesson-novice/blob/gh-pages/_episodes/01-introduction.md#the-lennard-jones-potential)




## Molecular Dynamics Software Available on Compute Canada Systems
### AMBER
> [Web page](http://ambermd.org)
### DL_POLY
>[Web page](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx)
### GROMACS
>[Web page](http://gromacs.org)
#### Force fields implemented in GROMACS:
- AMBER: 94, 96, 99, 99sb, 99sb-ildn, 03, GS (amberGS is amber94 with both backbone torsion potentials set to 0).
- CHARMM: 27 (optimized for proteins and nucleic acids).
- [CHARMM for GROMACS](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (CHARMM36, CGenFF)
- GROMOS: 43a1, 43a2, 45a3, 53a5, 53a6, 54a7.
- OPLS-AA (OPLS-AA implemented in GROMACS is actually OPLS-AA/L. OPLS-AA/L uses OPLS-AA atom types with the torsions and impropers refitted to QM calculations at the HF/6-31G** level followed by single-point LMP2/cc-pVTZ(-f))
### NAMD
>[Web page](https://www.ks.uiuc.edu/Research/namd/)
#### Force fields implemented in NAMD:
- [AMBER](http://ambermd.org/AmberModels.php) (amber format topology prepared with AMBERTOOLS). Currently AMBER recommends the following force fields: ff14SB for proteins, OL15 for DNA, OL3 for RNA, GLYCAM_06j for carbohydrates, lipid17 for lipids, and a general force field gaff2.
- [CHARMM](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (charmm format, topology can be prepared with the included **psfgen** program)
- [CHARMM-Drude](http://mackerell.umaryland.edu/charmm_drude_ff.shtml) (the polarizable force field based on the classical Drude oscillator model, charmm format)
- GROMOS (from GROMACS distribution, topology prepared with pdb2gmx)
- [OPLS-AA/M](http://zarbi.chem.yale.edu/oplsaam.html) New peptide dihedral parameters, significantly outperform the previous versions for proteins  (charmm format)


----------------------------------------------------------
1. The organic solvents with OPLS force field generate slightly better properties than those with GAFF. [ C. Caleman, P.J. van Maaren, M. Hong, J. S. Hub, L. T. Costa and D. van der Spoel. Force Field Benchmark of Organic Liquids: Density, Enthalpy of Vaporization, Heat Capacities, Surface Tension, Isothermal Compressibility, Volumetric Expansion Coefficient, and Dielectric Constant, J. Chem. Theor. Comput., 8, 61-74 (2012) ].
### LAMMPS
>[Web page](https://lammps.sandia.gov)

## Water Models
OPC family water models: OPC, OPC3
The accuracy of OPC water model is dramatically better compared to the commonly used rigid models.

Good pdb files for the tutorial:
1bvi
1de3
1goa *
1h4g
1lni *
