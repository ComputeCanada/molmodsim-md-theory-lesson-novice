---
title: "Supplemental: Overview of the Common Force Fields"
teaching: 0
exercises: 0
questions:
- "What are the categories of molecular dynamics force fields?"
- "What force fields are available?" 
- "What systems they work best with?"
objectives:
- "Be able to recognize the strengths and weaknesses of different types of force fields."
- "Find out which force fields are available and which systems they are most suitable for."
- "Identify the original papers that introduced force fields."
keypoints:
- "There are different types of force fields designed for different types of simulations."
- "Induction effects are not accounted for by fixed-charge force fields."
- "Using more accurate and diverse target data allows MD force fields to be improved."
---

* Table of Contents
{:toc}

## Early Force Fields

#### CFF (Consistent Force Field) (1968)
CFF was the first modern force field [(Lifson, 1968)]({{ page.root }}/reference.html#lifson-1968) Introduced a concept of a 'consistent force field'. Introduced a methodology for deriving and validating force fields. The term 'consistent' emphasized importance of the ability to describe a wide range of compounds and physical observables (conformation, crystal structure, thermodynamic properties and vibrational spectra). After the initial derivation for hydrocarbons CFF was extended to proteins, but was very crude at that time.

#### Allinger Force Fields (MM1 - MM4) (1976-1996)
- MM1 [(N.L. Allinger, 1976)](https://www.sciencedirect.com/science/article/abs/pii/S0065316008602129) - Class 1
- MM2 [(N.L. Allinger, 1977)](https://pubs.acs.org/doi/10.1021/ja00467a001) - Class 1
- MM3 [N.L.Allinger et al., 1989](https://pubs.acs.org/doi/10.1021/ja00205a001) - Class 2
- MM4 [(N.L.Allinger et al., 1996)](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<642::AID-JCC6>3.0.CO;2-U) - Class 3

Target data included electron diffraction, vibrational spectra, heats of formation, and crystal structures. The calculation results were verified via comparison with high-level _ab-initio_ quantum chemistry computations, and the parameters were additionally adjusted. Intended for calculations on small and medium size organic molecules.

#### ECEPP (Empirical Conformational Energy Program for Peptides) (1975)
ECEPP was the first force field targeting polypeptides and proteins. Crystal data of small organic compounds and semi-empirical QM calculations were used extensively in derivation of this force field. As more experimental data became available, force field parameters have been refined and modified versions were published.
- ECEPP [(F. A. Momany, et al., 1975)](https://pubs.acs.org/doi/10.1021/j100589a006#) 
- UNICEPP [(L.G. Dunfield, 1978)](https://pubs.acs.org/doi/10.1021/j100513a014) The united atoms version of ECEPP, developed for the conformational analysis of large molecules
- ECEPP-02 [(G. Nemethy et al, 1983)](https://pubs.acs.org/doi/10.1021/j100234a011)
- ECEPP-03 [(G. Nemethy et al, 1992)](https://pubs.acs.org/doi/10.1021/j100194a068)
- ECEPP-05 [(Y. A. Arnautova et al., 2006)](https://pubs.acs.org/doi/full/10.1021/jp054994x)


## Evolution of Force Fields

### Early force field advances
Early force field development focused on developing mathematical forms for MM energy function and methods of deriving parameters. Researchers investigated various forms of potential energy functions, and experimented with hydrogen bonding potential, combination rules, and out of plane angle potentials during this period.

#### United atoms force fields. 
United Atoms Model was developed to speed up large-scale simulations. It represents nonpolar carbons and their bonded hydrogens as a single particle. United Atoms force fields can significantly reduce the size of most simulations, since roughly half of the atoms in biological or other organic macromolecules are hydrogens. Additional advantage is the efficiency gain in conformational sampling. The first united atoms force field was UNICEPP
[(L.G. Dunfield, 1978)](https://pubs.acs.org/doi/10.1021/j100513a014).

According to early comparisons between all-atom and united-atom simulations, united-atom force fields adequately represent molecular vibrations and bulk properties of small molecules.  

After this initial success all major developers of protein force fields implemented united atoms models:
- [CHARMM19](https://aip.scitation.org/doi/pdf/10.1063/1.472061) (1985) 
- GROMOS
- [AMBER-UA](https://pubs.acs.org/doi/10.1021/ja00315a051) (1984)
- [OPLS-UA](https://pubs.acs.org/doi/10.1021/ja00214a001) (1988)

It became apparent, however, that there were some limitations:
- In the absence of explicit hydrogens, hydrogen bonds cannot be accurately treated;  
- π-stacking cannot be represented without explicitly including hydrogens in aromatic groups; 
- when hydrogens were combined with polar heavy atoms, dipole and quadrupole moments were found to be inaccurate.
 
New approaches were found to overcome the limitations of united-atom force fields. For example, only aliphatic hydrogens, which are not significantly charged and do not participate in hydrogen bonds, are represented as united atoms while other hydrogens are represented explicitly. In this way, the limitations of the united-atom force field are partially mitigated while preserving most of the benefits of the united-atom force field.


### Refinement after the initial introduction.

Statistical errors caused by relatively short simulation lengths and systematic errors caused by inaccurate force fields limit the predictive power of MD simulations. Deficiencies of the force fields remained undetected when  statistical errors caused  by insufficient sampling prevailed. Increase of the computing power over last two decades allowed for much longer simulations and resulted in a significant reduction of statistical errors. This led to the detection of force field deficiencies such as large deviations in different observables and inability to predict conformations of proteins and peptides. Various approaches were undertaken to improve force fields such as:

- **Use an all-atom model to increase accuracy.** Most force fields (CHARMM22, AMBER ff99 and GAFF, OPLS-AA, OPLS-AA/L) converted back to all-atom model.

- **Increase the size of the target data.** Because each force field was derived with a different training set of atomic configurations, it was biased in one way or another. By using large atomic reference sets and careful selection of the atomic configurations, bias problems were reduced and accuracy was improved.

- **Improve representation of static electrostatic potential.** The atom-centered point charge model has two shortcomings: 
    - It is unable to describe the anisotropy of the charge distribution;
    - It does not account for the charge penetration effect (deviation of electrostatic interaction from the Coulomb form due to electron shielding when atomic electron clouds overlap).

    In molecular complexes, these effects determine equilibrium geometry and energy. Examples of anisotropic charge distribution are σ-holes, lone pairs, and aromatic systems.

    In covalent bonding an atom's side opposite its bond usually has a lower electron density region known as σ-hole. Through a positive electrostatic potential associated with a sigma-hole, an atom can interact attractively with negative sites.

    A simple solution to the σ-hole model is to attach an off-centered positive charge to the halogen atom (
   [W. Jorgensen & P. Schyman, 2012](https://pubs.acs.org/doi/10.1021/ct300180w), [F.Lin, A. MacKerell, 2018](https://sci-hub.se/10.1021/acs.jctc.7b01086)). An atomic multipole method provides a more thorough way of describing anisotropic charge distribution ([J.Ponder et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2918242/)).

- **Include polarization effects.** Early force fields employed fixed atomic charges to model the electrostatic interactions. Fixed-charge electrostatics does not account for the many-body polarization that can vary significantly depending on chemical and physical environments. Consequently, non-polarizable force fields fail to capture the conformational dependence of electrostatic properties. Polarizable force fields where the charges can be calculated from the energy equilibration have been developed ot address this problem ([see a recent review](https://www.annualreviews.org/doi/10.1146/annurev-biophys-070317-033349)). A drawback of including polarization is that simulations take longer to run due to the high computational cost. 


### Recent developments in force fields and future prospects
Despite extensive attempts to improve force files, they have often failed to achieve quantitative accuracy. There was a realization that the functional form of potential energy was the major problem with all widely used protein force fields.

Two strategies can be used to address this problem:

1. Expand and improve the rigor of the representation of the underlying physics.
2. Develop empirical corrections to compensate for deficiency of physical representation

AMBER, CHARMM, and OPLS focused their efforts on empirical correction of the simple potential function.  AMOEBA and COMPASS worked on improving the functional form.

It is well understood that charge distributions are affected by both chemical environments and local geometry changes. As we discussed above, the former is explicitly treated in polarizable force field models. AMOEBA+ (2019) improved representation of electrostatic interactions by incorporating charge penetration and intermolecular charge transfer.

Dependence of atomic charges on the local molecular geometry is ignored by almost all classical force fields even though it is well known that it causes issues. 

One of the few force fields that consider CF contributions is SDFF:

https://eurekamag.com/research/011/176/011176666.php

Palmo, K.; Mannfors, B.; Mirkin, N. G.; Krimm, S. Inclusion of charge and polarizability fluxes provides needed physical accuracy in molecular mechanics force fields. Chem. Phys. Lett. 2006, 429 (4), 628.

 Recently Geometry-Dependent Charge Flux (GDCF) model considering charge flux contributions along bond and angle was introduced in AMOEBA+(CF) force field. 

[Implementation of Geometry-Dependent Charge Flux into the Polarizable AMOEBA+ Potential](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.9b03489)

## Force Fields Aimed at Improving Quality of Molecular Interactions
#### CVFF (Consistent Valence Force Field) [(Maple & Hagler, 1988)](https://www.pnas.org/doi/epdf/10.1073/pnas.85.15.5350)
#### CFF93 (An ab initio all-atom force field for polycarbonates) [(Sun et al., 1994)](https://pubs.acs.org/doi/10.1021/ja00086a030)
#### CFF, formerly CFF95 [(Jonsdottir & Rasmussen, 2000)](https://doi.org/10.1039/A909475J)
#### MM3 [(N.L.Allinger et al., 1989)](https://pubs.acs.org/doi/10.1021/ja00205a001), MM4 [(N.L.Allinger et al., 1996)](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<642::AID-JCC6>3.0.CO;2-U) 

#### COMPASS (Condensed-phase Optimized Molecular Potentials for Atomistic Simulation Studies) 
##### COMPASS [(H. Sun, 1998)](https://pubs.acs.org/doi/10.1021/jp980939v)
 Developed for simulations of organic molecules, inorganic small molecules, and polymers

##### COMPASS II [(H. Sun et al., 2016)](https://link.springer.com/article/10.1007/s00894-016-2909-0) 
Extended the coverage to polymer and drug-like molecules found in popular databases. The VDW parameters are obtained by fitting enthalpies of vaporization and densities, to experimental data. The atomic partial charges are derived using QM and empirically adjusted to take hydrogen bonding effects into account. The COMPASS energy function offers six types of cross-terms: bond-bond, bond-angle, angle-angle, bond-torsion, angle-torsion, and angle-torsion-angle.


## Biomolecular Force Fields for Large-Scale Simulations.

Long simulations of large systems are required to study biologically relevant processes. It takes a lot of computer power to run such simulations. Due to this, the main objective is to develop a minimal force field that allows simulation sizes and time periods to be extended as much as possible while still keeping chemical structures, interaction energies, and thermodynamic properties within an acceptable range.

The workhorses of modern biomolecular simulations are all-atom, fixed-charge force fields:

* **AMBER** (Assisted Model Building with Energy Refinement)
* **CHARMM** (Chemistry at HARvard Macromolecular Mechanics) 
* **GROMOS** (GROningen MOlecular Simulation)
* **OPLS** (Optimized Potentials for Liquid Simulations) 

CHARMM and AMBER forcefields are developed for simulations of proteins and nucleic acids, and they focus on accurate description of structures and non-bonded energies. GROMOS and OPLS are geared toward accurate description of thermodynamic properties such as heats of vaporization, liquid densities, and molecular solvation properties. 

#### AMBER 
AMBER forcefields are developed for simulations of proteins and nucleic acids and they are focused on accurate description of structures and non-bonded energies. The VDW parameters are obtained from crystal structures and lattice energies. The atomic partial charges are fitted to QM electrostatic potential without any empirical adjustments.

##### [ff99](http://dx.doi.org/doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F) (1999) 
The main idea is that the use of RESP charges, should lead to the need for fewer torsional potentials than in models with an empirical charge derivation. Potential energy surface scans were performed using four different ab initio methods, HF/6‐31G<sup>*</sup>, MP2/6‐31G<sup>*</sup>, MP2/6‐311+G (2d,p), and B3LYP/6‐311+G (2d,p).

Wang J, Cieplak P, Kollman PA. How well does a restrained electrostatic potential (RESP) model perform in calculating conformational energies of organic and biological molecules? Journal of Computational Chemistry. 2000;21: 1049–1074. doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F

After publication of ff99 a number of studies devoted primarily to modifying the torsion potentials in order to correct the observed discrepancies have been published:

##### [ff03](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.10349) (2003)
Completely new charge set developed using B3LYP/cc-pvTZ quantum calculations in a polarizable continuum (PCM) solvent intended to mimic the interior of a protein. Backbone torsion Fourier series were derived specifically for this new charge set at the MP2/cc-pvTZ level of theory, also in the context of PCM solvent. MM gas-phase energies computed with charges derived in PCM solvent have been shown to double-count polarization effects, and ff03 force field has not become as widely used or further refined as ff99

Duan, Yong, Chun Wu, Shibasish Chowdhury, Mathew C. Lee, Guoming Xiong, Wei Zhang, Rong Yang, et al. “A point-charge force field for molecular mechanics simulations of proteins based on condensed-phase quantum mechanical calculations.” Journal of Computational Chemistry 24, no. 16 (2003): 1999–2012. https://doi.org/10.1002/jcc.10349.

##### [ff99sb](https://doi.org/10.1002/prot.21123) (2006) 
Optimized for the correct description of the helix-coil equilibrium

##### ff99SB-φ' 
Targeted the reproduction of the intrinsic conformational preferences of tripeptides

##### [ff99SBnmr](https://doi.org/10.1002/anie.201001898) (2010) and ff99SB_φψ
Target data included protein NMR chemical shifts and residual dipolar couplings.

##### [ff99SBildn](https://doi.org/10.1002/prot.22711) (2010) 
Targeted optimization of four amino acid side chains.

Lindorff‐Larsen, Kresten, Stefano Piana, Kim Palmo, Paul Maragakis, John L. Klepeis, Ron O. Dror, and David E. Shaw. “Improved Side-Chain Torsion Potentials for the Amber Ff99SB Protein Force Field.” Proteins: Structure, Function, and Bioinformatics 78, no. 8 (2010): 1950–58. doi:10.1002/prot.22711.

##### [ff12SB](https://doi.org/10.1021/acs.jctc.5b00255) (2012) 
This is the preliminary version of ff14SB described in the same paper.

##### [ff12SB-cMAP](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00662) (2015)
Force field for implicit-solvent simulations.

##### [ff99IDPs](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00043) (2015)
Force field for intrinsically disordered proteins

##### [ff14ipq](https://doi.org/10.1021/ct500643c) (2014) 
This force field is derived using another new approach aiming to derive charges implicitly polarized by the fixed charge explicit water (IPolQ method). The charges are derived in  self-consistent manner in presence of explicit water molecules represented by TIP4P-Ew water model at MP2/cc-pV(T+d)Z level. The weak point of the ff14ipq force field is overstabilization of salt bridges.

Cerutti, David S., William C. Swope, Julia E. Rice, and David A. Case. *Ff14ipq: A Self-Consistent Force Field for Condensed-Phase Simulations of Proteins.* Journal of Chemical Theory and Computation 10, no. 10 (October 14, 2014): 4515–34. doi:10.1021/ct500643c.

##### [ff14SB](https://doi.org/10.1021/acs.jctc.5b00255) (2015) 
Another attempt to corect for deficiencies of ff99SB by using new side chain dihedral parameters based on MP2 level calculations followed by adjustment to the backbone φ energy profile. Used the old ff99SB charge set.

Maier, James A., Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling. “Ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from Ff99SB.” Journal of Chemical Theory and Computation 11, no. 8 (August 11, 2015): 3696–3713. https://doi.org/10.1021/acs.jctc.5b00255.

##### [ff15ipq](https://doi.org/10.1021/acs.jctc.6b00567) (2016) 
The motivation for the development was to address the salt bridge overstabilization problem of ff14ipq. However, this forcefiled is not a small correction applied to ff14ipq. It is is a complete semi-automatic rederivation of all parameters with the different water model SPC/Eb. The salt bridge overstabilization issue was addressed by increasing radius of polar hydrogens bonded to nitrogen.  The modified FF performed as well or better than the other fixed charge force fields. Polarizable CHARMM Drude-2013 and AMOEBA performed better in this respect, as expected.

Problems related to protein stability persist. Even 4 μs simulations “were not sufficiently long to obtain converged estimates of secondary structure of polypeptides”. In simulation tests some proteins deviated significantly near the end of several microsecond simulations, and it is not clear whether this is a transient fluctuation or transition to a different state.

Debiec, Karl T., David S. Cerutti, Lewis R. Baker, Angela M. Gronenborn, David A. Case, and Lillian T. Chong. *Further along the Road Less Traveled: AMBER Ff15ipq, an Original Protein Force Field Built on a Self-Consistent Physical Model.* Journal of Chemical Theory and Computation 12, no. 8 (August 9, 2016): 3926–47. doi:10.1021/acs.jctc.6b00567

##### [ff15FB](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b02320) (2017) 
The goal was to reoptimize intramolecular bond, angle, and dihedral parameters to fit MP2 level QM data without modifying the nonbonded parameters. The parameter optimization was done using [ForceBalance](https://pubs.acs.org/doi/10.1021/jz500737m) package. Same peformance as the earlier versions for equilibrium properties,  improved performance for temperature dependence. For best agreement with experiment recommended to use with the TIP3P-FB water. The TIP3P-FB model was parametrized to reproduce the temperature and pressure dependence of a wide range of thermodynamic properties.

##### [ff19SB](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00591) (2020)
Backbone dihedral parameters optimized using as reference data the entire 2D quantum mechanics (QM) energy surface. Both QM and MM calculations were done in aqueous solution. AMBER ff19SB uses CMAP torsional potentials. The authors concluded that "ff19SB, when combined with a more accurate water model such as OPC, should have better predictive power for modeling sequence-specific behavior, protein mutations, and also rational protein design".

#### CHARMM 

##### [CHARMM19](https://aip.scitation.org/doi/pdf/10.1063/1.472061) (1996)
United-atom model, originally released in 1985. 

##### [CHARMM22](https://onlinelibrary.wiley.com/doi/10.1002/jcc.20077) (1998)
All-atom model.

##### <a href="https://doi.org/10.1002/(SICI)1096-987X(20000130)21:2%3C105::AID-JCC3%3E3.0.CO;2-P">CHARMM27</a> (2000)
Update to the CHARMM22 featuring the optimization of nucleic acid and lipid parameters, as well as the introduction of a number of new ions.

##### [CHARMM22/CMAP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1367299/#bib1) (2004)
This forcefield introduces a tabulated correction for the φ-, ψ-angular dependence of the potential energy. As a result of the application of CMAP correction, the dynamical and structural properties of proteins were significantly improved in molecular dynamics simulations.

##### [CHARMM27r](https://pubs.acs.org/doi/10.1021/jp0468096) (2005)
Improved lipid parameter set.

##### CHARMM35 (2008)
Carbohydrate parameter set.

##### [CGenFF](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21367) (2009)

##### CHARMM36 (2012)
refined backbone CMAP potentials and introduced new side-chain dihedral parameters. The updated CMAP corrected the C22/CMAP FF bias towards alpha-helical conformations.

##### [Review](https://doi.org/10.1016/j.bbagen.2014.08.004)

##### [CHARMM Force Field Files](http://mackerell.umaryland.edu/charmm_ff.shtml)

#### GROMOS

##### 53A5, 53A6 [(Oostenbrink et al., 2004)](https://onlinelibrary.wiley.com/doi/10.1002/jcc.20090)

##### 54A7, 54B7 [(N.Schmid et al., 2011)](https://link.springer.com/article/10.1007/s00249-011-0700-9)

##### 54A8 [(M.M.Reif et al., 2012)](https://pubs.acs.org/doi/10.1021/ct300156h)

Parameter sets for GROMOS force fields are specified in accordance with the following format: 

- The number of atom types
- A letter code that describes for what conditions the parameters are optimized (e.g. A - solution, B - vacuum).
- The version number 


#### OPLS
Force fields of the OPLS family are designed for simulating liquids that contain organic molecules and proteins. The VDW parameters are optimized using experimental liquid properties, mainly enthalpies of vaporization and densities. The atomic partial charges are derived using QM and experimental condensed-phase properties. An important part of the OPLS philosophy is balancing solvent-solvent and solute-solvent interactions. 

##### [OPLS-AA](https://pubs.acs.org/doi/10.1021/ja9621760) (1996) 
This is the first all-atom OPLS force field. Bond stretching and angle bending parameters are taken from the AMBER force field. The torsional parameters were fit to the RHF/6-31G* calculations of about 50 organic molecules and ions. The charges are empirical and have been obtained from fitting to reproduce properties of organic liquids.

Jorgensen W, Maxwell D, Tirado-Rives J. Development and testing of the OPLS all-atom force field on conformational energetics and properties of organic liquids. J Am Chem Soc. 1996;118: 11225–11236.

##### [OPLS-AA/L](https://pubs.acs.org/doi/10.1021/jp003919d) (2001) 
Large data set, more than 2000 data points of energies for the 20 amino acids based on geometry optimization at the HF/6-31G** level followed by single-point LMP2/cc-pVTZ(-f) calculations. This level of theory is accurate within 0.25 kcal/mol.

Kaminski GA, Friesner RA, Tirado-Rives J, Jorgensen WL. Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides. J Phys Chem B. 2001;105: 6474–6487. doi:10.1021/jp003919d

##### [OPLS_2005](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20292) (2005) 
Main goal - to broaden the coverage of OPLS_2001 and refine torsion parameters using a larger dataset. The new Data set included torsion profiles from 637 compounds.

Banks JL, Beard HS, Cao Y, Cho AE, Damm W, Farid R, et al. Integrated Modeling Program, Applied Chemical Theory (IMPACT). J Comput Chem. 2005;26: 1752–1780. doi:10.1002/jcc.20292

##### [OPLS-AAx & OPLS/CM1Ax](https://pubs.acs.org/doi/10.1021/ct300180w) (2012) 
The representation of chlorine, bromine, and iodine in aryl halides has been modified in the OPLS-AA and OPLS/CM1A force fields in order to incorporate halogen bonding. It was accomplished by adding a partial positive charge in the region of the σ-hole.


##### [OPLS2.0](http://dx.doi.org/doi:10.1021/ct300203w) (2012) 
OPLS2 was developed to improve the accuracy of drug-like molecules. Substantially expanded data set contained QM data on more than 11,000 molecules. CM1A-BCC charges were used. The BCC terms were parameterized against the OPLS-AA charges for a core set of 112 molecules and the electrostatic potential at the HF/6-31G* level. The BCC terms were empirically adjusted to minimize the errors with regard to the ASFE using a training set of 153 molecules.

Shivakumar D, Harder E, Damm W, Friesner RA, Sherman W. Improving the Prediction of Absolute Solvation Free Energies Using the Next Generation OPLS Force Field. J Chem Theory Comput. 2012;8: 2553–2558. doi:10.1021/ct300203w

##### [OPLS-AA/M](http://dx.doi.org/doi:10.1021/acs.jctc.5b00356) (2015) 
New torsion parameters using higher level ωB97X-D(23)/6-311++G(d,p) and  B2PLYP-D3BJ/aug-cc-pVTZ(26) QM calculations.

Robertson MJ, Tirado-Rives J, Jorgensen WL. Improved Peptide and Protein Torsional Energetics with the OPLS-AA Force Field. J Chem Theory Comput. 2015;11: 3499–3509. doi:10.1021/acs.jctc.5b00356

##### [OPLS3](http://dx.doi.org/doi:10.1021/acs.jctc.5b00864) (2016) 
Added off-atom charge sites to represent halogen bonding and aryl nitrogen lone pairs. Complete refit of peptide dihedral parameters using an order of magnitude more data. Claimed 30% improvement. Still the same original VDW parameters.

Harder E, Damm W, Maple J, Wu C, Reboul M, Xiang JY, et al. OPLS3: A Force Field Providing Broad Coverage of Drug-like Small Molecules and Proteins. J Chem Theory Comput. 2016;12: 281–296. doi:10.1021/acs.jctc.5b00864

## Polarizable Force fields
#### AMBER
##### ff02pol [(P. Cieplak et al., 2001)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.1065)
##### ff12pol [(J. Wang et al., 2011)](https://pubs.acs.org/doi/10.1021/jp112133g)

#### AMOEBA (Atomic Multipole Optimized Energetics for Biomolecular Applications) 
##### AMOEBA-2002 [(Ren and Ponder, 2002)](https://onlinelibrary.wiley.com/doi/10.1002/jcc.10127)
##### AMOEBA-2013 [(Shi et al., 2013)](https://pubs.acs.org/doi/10.1021/ct4003702)

AMOEBA-2013 uses permanent electrostatic multipole (dipole and quadrupole) moments at each atom and explicitly treats polarization effects under various chemical and physical conditions.

##### AMOEBA+ [(C. Liu et al., 2019)](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00261)

##### AMOEBA+(CF) [(C. Liu et al., 2020)](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.9b03489)

##### Q-AMOEBA [(N.Mauger et al., 2022)](https://pubs.acs.org/doi/10.1021/acs.jpcb.2c04454)

#### OPLS
##### OPLS/PFF https://onlinelibrary.wiley.com/doi/10.1002/jcc.10125

Force field for water with polarization capabilities. Intended for simulations that explicitly take nuclear quantum effects into account.

#### CHARMM 
##### CHARMM fluctuating charge model [(S.Patel et al., 2004)](https://onlinelibrary.wiley.com/doi/10.1002/jcc.20077)
##### CHARMM Drude model [(P.E.M. Lopez, 2013)](https://pubs.acs.org/doi/10.1021/ct400781b)

## Force Fields for Small Molecules
##### GAFF
##### CGenFF 

1. [(F.-U.Lin &  AD. MacKerell, 2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6733265/) Force Fields for Small Molecules

## Additional Reading

1. [(Dauber-Osguthorpe, 2019)](https://link.springer.com/article/10.1007/s10822-018-0111-4) Biomolecular force fields: where have we been, where are we now, where do we need to go and how do we get there? - Review of the origins of FF based calculations, theory and methodology of FF development.

2. [(Hagler, 2019)](https://doi.org/doi:10.1007/s10822-018-0134-x) Force field development phase II: Relaxation of physics-based criteria… or inclusion of more rigorous physics into the representation of molecular energetics. - The latest developments in improvement of FF accuracy and robustness are discussed.

3. [Tinker-HP](https://tinker-hp.org) is a multi-CPUs and multi-GPUs/multi-precision, MPI massively parallel package dedicated to long molecular dynamics simulations with classical and polarizable force fields, neural networks and advanced QM/MM.


