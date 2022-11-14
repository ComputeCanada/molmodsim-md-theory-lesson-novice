---
title: "Supplemental: Overview of common Force Fields"
teaching: 0
exercises: 0
questions:
- "What are the categories of molecular dynamics force fields?"
objectives:
- "Be able to recognize the strengths and weaknesses of different types of force fields"
keypoints:
- "Using more accurate and diverse target data allows MD force fields to be improved."
---

#### Table of Contents

* Table of Contents
{:toc}

## History of Force Field Development

### Consistent Force Field (CFF)
CFF was the first modern force field [(Lifson, 1968)]({{ page.root }}/reference.html#lifson-1968) Introduced a concept of a 'consistent force field'. Introduced a methodology for deriving and validating force fields. The term 'consistent' emphasized importance of the ability to describe a wide range of compounds and physical observables (conformation, crystal structure, thermodynamic properties and vibrational spectra). After the initial derivation for hydrocarbons CFF was extended to proteins, but was very crude at that time.

### Allinger force fields: MM1, MM2, MM3, MM4.
MM1, MM2 - Class 1
MM3 - Class 2
MM4 - Class 3
Target data included electron diffraction, vibrational spectra, heats of formation, and crystal structures. The calculation results were verified via comparison with high-level _ab-initio_ quantum chemistry computations, and the parameters were additionally adjusted.

Intended for calculations on small and medium size organic molecules.

### Empirical Conformational Energy Program for Peptides (ECEPP)
ECEPP was the first force field targeting polypeptides and proteins. Crystal data of small organic compounds and semi-empirical QM calculations were used extensively in derivation of this force field. As more experimental data became available, force field parameters have been refined and modified versions ECEPP-2, ECEPP-3 were published.

Over the years of evaluations the problems originating from the deficiency of these FF became apparent and various approaches has been undertaken to correct them.

### Force Fields with development focused on application to proteins:
* AMBER
* CHARMM
* GROMOS
* OPLS

GROMOS and OPLS are focused on fitting to thermodynamic properties such as the heats of vaporization, liquid densities, and the solvation properties of small molecules.

Most relevant processes require very long simulations. Large systems create computational restraints. The goal is to develop a "minimalistic" force field to expand simulation time window as much as possible.

### Force fields with development focused on improving representation of molecular interactions
* CVFF
* CFF93
* CFF95
* MM1 - MM4

Development Phase I

Refinement after the initial introduction.
-Converted back to AA, except GROMOS

CHARMM22
AMBER, ff99, GAFF
OPLS-AA, OPLS-AA/L

Used large datasets for training. Training datasets were different for different FF.

Large deviations in different observables. Inability to predict conformations of peptides and proteins.
Simple 12‑6‑1 quadratic diagonal FFs (as used
in standard biomolecular FFs) are not adequate
to achieve quantitative accuracy. A major problem with all widely used protein force fields is the functional form of the potential energy.

2 paths:

1. expand and improve the rigor of the representation of the underlying physics.
2. Develop empirical corrections to compensate for deficiency of physical representation

Unaccounted physics:
-atomic charges depend on the geometry (charge flux)

AMBER, CHARMM, OPLS focused their efforts on empirical correction of the simple potential function


### COMPASS
COMPASS and COMPASS II forcefields are developed for simulations of organic molecules, inorganic small molecules, and polymers. The VDW parameters are obtained by fitting enthalpies of vaporization and densities, to experimental data. The atomic partial charges are derived using QM and empirically adjusted to take hydrogen bonding effects into account.

### OPLS
OPLS family force fields are created for liquid simulations containing organic molecules and proteins. The VDW parameters are optimized using experimental liquid properties, mainly enthalpies of vaporization and densities. The atomic partial charges are derived using QM and experimental condensed-phase properties.

[OPLS-AA](https://pubs.acs.org/doi/10.1021/ja9621760) (1996) This is the first all-atom OPLS force field. Bond stretching and angle bending parameters are taken from the AMBER force field. The torsional parameters were fit to the RHF/6-31G* calculations of about 50 organic molecules and ions. The charges are empirical and have been obtained from fitting to reproduce properties of organic liquids.

Jorgensen W, Maxwell D, Tirado-Rives J. Development and testing of the OPLS all-atom force field on conformational energetics and properties of organic liquids. J Am Chem Soc. 1996;118: 11225–11236.

[OPLS-AA/L](https://pubs.acs.org/doi/10.1021/jp003919d) (2001) Large data set, more than 2000 data points of energies for the 20 amino acids based on geometry optimization at the HF/6-31G** level followed by single-point LMP2/cc-pVTZ(-f) calculations. This level of theory is accurate within 0.25 kcal/mol.

Kaminski GA, Friesner RA, Tirado-Rives J, Jorgensen WL. Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides. J Phys Chem B. 2001;105: 6474–6487. doi:10.1021/jp003919d

[OPLS_2005](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20292) (2005) Main goal - to broaden the coverage of OPLS_2001 and refine torsion parameters using a larger dataset. The new Data set included torsion profiles from 637 compounds.

Banks JL, Beard HS, Cao Y, Cho AE, Damm W, Farid R, et al. Integrated Modeling Program, Applied Chemical Theory (IMPACT). J Comput Chem. 2005;26: 1752–1780. doi:10.1002/jcc.20292

[OPLS2.0](http://dx.doi.org/doi:10.1021/ct300203w) (2012) OPLS2 was developed to improve the accuracy of drug-like molecules. Substantially expanded data set contained QM data on more than 11,000 molecules. CM1A-BCC charges were used. The BCC terms were parameterized against the OPLS-AA charges for a core set of 112 molecules and the electrostatic potential at the HF/6-31G* level. The BCC terms were empirically adjusted to minimize the errors with regard to the ASFE using a training set of 153 molecules.

Shivakumar D, Harder E, Damm W, Friesner RA, Sherman W. Improving the Prediction of Absolute Solvation Free Energies Using the Next Generation OPLS Force Field. J Chem Theory Comput. 2012;8: 2553–2558. doi:10.1021/ct300203w

[OPLS-AA/M](http://dx.doi.org/doi:10.1021/acs.jctc.5b00356) (2015) New torsion parameters using higher level ωB97X-D(23)/6-311++G(d,p) and  B2PLYP-D3BJ/aug-cc-pVTZ(26) QM calculations.

Robertson MJ, Tirado-Rives J, Jorgensen WL. Improved Peptide and Protein Torsional Energetics with the OPLS-AA Force Field. J Chem Theory Comput. 2015;11: 3499–3509. doi:10.1021/acs.jctc.5b00356

[OPLS3](http://dx.doi.org/doi:10.1021/acs.jctc.5b00864) (2016) Added off-atom charge sites to represent halogen bonding and aryl nitrogen lone pairs. Complete refit of peptide dihedral parameters using an order of magnitude more data. Claimed 30% improvement. Still the same original VDW parameters.

Harder E, Damm W, Maple J, Wu C, Reboul M, Xiang JY, et al. OPLS3: A Force Field Providing Broad Coverage of Drug-like Small Molecules and Proteins. J Chem Theory Comput. 2016;12: 281–296. doi:10.1021/acs.jctc.5b00864

### AMBER force fields
AMBER forcefields are developed for simulations of proteins and nucleic acids and they are focused on accurate description of structures and non-bonded energies. The VDW parameters are obtained from crystal structures and lattice energies. The atomic partial charges are fitted to QM electrostatic potential without any empirical adjustments.

[ff99](http://dx.doi.org/doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F) (1999) The main idea is that the use of RESP charges, should lead to the need for fewer torsional potentials than in models with an empirical charge derivation. Potential energy surface scans were performed using four different ab initio methods, HF/6‐31G$*$, MP2/6‐31G$*$, MP2/6‐311+G (2d,p), and B3LYP/6‐311+G (2d,p).

Wang J, Cieplak P, Kollman PA. How well does a restrained electrostatic potential (RESP) model perform in calculating conformational energies of organic and biological molecules? Journal of Computational Chemistry. 2000;21: 1049–1074. doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F

After publication of ff99 a number of studies devoted primarily to modifying the torsion potentials in order to correct the observed discrepancies have been published:

[ff03](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.10349) (2003) Completely new charge set developed using B3LYP/cc-pvTZ quantum calculations in a polarizable continuum (PCM) solvent intended to mimic the interior of a protein. Backbone torsion Fourier series were derived specifically for this new charge set at the MP2/cc-pvTZ level of theory, also in the context of PCM solvent. MM gas-phase energies computed with charges derived in PCM solvent have been shown to double-count polarization effects, and ff03 force field has not become as widely used or further refined as ff99

Duan, Yong, Chun Wu, Shibasish Chowdhury, Mathew C. Lee, Guoming Xiong, Wei Zhang, Rong Yang, et al. “A point-charge force field for molecular mechanics simulations of proteins based on condensed-phase quantum mechanical calculations.” Journal of Computational Chemistry 24, no. 16 (2003): 1999–2012. https://doi.org/10.1002/jcc.10349.

[ff99sb] (https://doi.org/10.1002/prot.21123) 2006 is optimized for the correct description of the helix-coil equilibrium

**ff99SB-φ'**  targeted the reproduction of the intrinsic conformational preferences of tripeptides

[ff99SBnmr](https://doi.org/10.1002/anie.201001898) 2010 and FF99SB_φψ  target data during included protein NMR chemical shifts and residual dipolar couplings.

[ff99SBildn](https://doi.org/10.1002/prot.22711) 2010 targeted optimization of four amino acid side chains.

Lindorff‐Larsen, Kresten, Stefano Piana, Kim Palmo, Paul Maragakis, John L. Klepeis, Ron O. Dror, and David E. Shaw. “Improved Side-Chain Torsion Potentials for the Amber Ff99SB Protein Force Field.” Proteins: Structure, Function, and Bioinformatics 78, no. 8 (2010): 1950–58. doi:10.1002/prot.22711.

[ff12SB](https://doi.org/10.1021/acs.jctc.5b00255) 2012 This is the preliminary version of ff14SB described in the same paper.

[ff12SB-cMAP](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00662) 2015 force field for implicit-solvent simulations.

[ff99IDPs](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00043) 2015 force field for intrinsically disordered proteins

[ff14ipq](https://doi.org/10.1021/ct500643c) 2014 This force field is derived using another new approach aiming to derive charges implicitly polarized by the fixed charge explicit water (IPolQ method). The charges are derived in  self-consistent manner in presence of explicit water molecules represented by TIP4P-Ew water model at MP2/cc-pV(T+d)Z level. The weak point of the ff14ipq force field is overstabilization of salt bridges.

Cerutti-2014
:   Cerutti, David S., William C. Swope, Julia E. Rice, and David A. Case. *Ff14ipq: A Self-Consistent Force Field for Condensed-Phase Simulations of Proteins.* Journal of Chemical Theory and Computation 10, no. 10 (October 14, 2014): 4515–34. doi:10.1021/ct500643c.

[ff14SB](https://doi.org/10.1021/acs.jctc.5b00255) 2015 Another attempt to corect for deficiencies of ff99SB by using new side chain dihedral parameters based on MP2 level calculations followed by adjustment to the backbone φ energy profile. Used the old ff99SB charge set.

Maier, James A., Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling. “Ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from Ff99SB.” Journal of Chemical Theory and Computation 11, no. 8 (August 11, 2015): 3696–3713. https://doi.org/10.1021/acs.jctc.5b00255.

[ff15ipq](https://doi.org/10.1021/acs.jctc.6b00567) 2016 The motivation for the development was to address the salt bridge overstabilization problem of ff14ipq. However, this forcefiled is not a small correction applied to ff14ipq. It is is a complete semi-automatic rederivation of all parameters with the different water model SPC/Eb. The salt bridge overstabilization issue was addressed by increasing radius of polar hydrogens bonded to nitrogen.  The modified FF performed as well or better than the other fixed charge force fields. Polarizable CHARMM Drude-2013 and AMOEBA performed better in this respect, as expected.

Problems related to protein stability persist. Even 4 μs simulations “were not sufficiently long to obtain converged estimates of secondary structure of polypeptides”. In simulation tests some proteins deviated significantly near the end of several microsecond simulations, and it is not clear whether this is a transient fluctuation or transition to a different state.

Debiec-2016
:   Debiec, Karl T., David S. Cerutti, Lewis R. Baker, Angela M. Gronenborn, David A. Case, and Lillian T. Chong. *Further along the Road Less Traveled: AMBER Ff15ipq, an Original Protein Force Field Built on a Self-Consistent Physical Model.* Journal of Chemical Theory and Computation 12, no. 8 (August 9, 2016): 3926–47. doi:10.1021/acs.jctc.6b00567

[ff15FB](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b02320) 2017 The goal was to reoptimize intramolecular bond, angle, and dihedral parameters to fit MP2 level QM data without modifying the nonbonded parameters. The parameter optimization was done using [ForceBalance](https://pubs.acs.org/doi/10.1021/jz500737m) package. Same peformance as the earlier versions for equilibrium properties,  improved performance for temperature dependence. For best agreement with experiment recommended to use with the TIP3P-FB water. The TIP3P-FB model was parametrized to reproduce the temperature and pressure dependence of a wide range of thermodynamic properties.

[ff19SB]https://pubs.acs.org/doi/10.1021/acs.jctc.9b00591. 2020 Backbone dihedral parameters optimized using as reference data the entire 2D quantum mechanics (QM) energy surface. Both QM and MM calculations were done in aqueous solution. AMBER ff19SB uses CMAP torsional potentials. The authors concluded that "ff19SB, when combined with a more accurate water model such as OPC, should have better predictive power for modeling sequence-specific behavior, protein mutations, and also rational protein design".

### CHARMM

CHARMM22

CHARMM22/CMAP (CHARMM27)

CHARMM36 refined backbone CMAP potentials and introduced new side-chain dihedral parameters. The updated CMAP corrected the C22/CMAP FF bias towards alpha-helical conformations.


### Polarizable Force fields
**Drude-2013**
**AMOEBA**

The origins of FF based calculations, theory and methodology of FF development have been recently reviewed in [(Dauber-Osguthorpe, 2019)]({{ page.root }}/reference.html#dauber-osguthorpe-2019), and the latest developments in improvement of FF accuracy and robustness are discussed in [(Hagler, 2019)]({{ page.root }}/reference.html#hagler-2019).


## How to create ligand topology

[Automated Topology Builder](https://atb.uq.edu.au/index.py)
