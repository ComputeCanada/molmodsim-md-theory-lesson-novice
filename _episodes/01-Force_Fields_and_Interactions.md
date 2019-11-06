---
title: "Force Fields and Interactions"
teaching: 30
exercises: 0
questions:
- "What is Molecular Dynamics and how can I benefit from using it?"
- "What is a force field?"
- "What mathematical functions are used to approximate interactions between atoms and molecules in MD simulations?"
objectives:
- "Explain strengths and weaknesses of MD simulations"
- "Explain the interactions between particles in a molecular dynamics simulation"
keypoints:
- "Molecular dynamics simulates atomic positions in time using forces acting between atoms"
- "The forces in molecular dynamics are approximated by simple empirical functions"
---

Molecular dynamics (MD) simulations are widely used as tools to investigate structure and dynamics of proteins, nucleic acids, carbohydrates, lipids, nanoparticles and liquid/solid interfaces under a wide variety of conditions. MD is the simulation of atomic positions in time accomplished by solving classical Newton's equation of motion stating that the rate of change of momentum $$ \vec{p} $$ of an object equals the force $$ \vec{F} $$ acting on it:

$ \vec{F}=\frac{d\vec{p}}{dt} $

The potential energy function *U* is a cornerstone of the MD simulations because it allows calculating the forces. The force on an object is the negative of the derivative of the potential energy function:

$\vec{F}=-\nabla{U}(\vec{r})$

- Molecular dynamics helps to link physics, chemistry and biology
- With the help of MD it is possible to model phenomena that cannot be studied experimentally.
-  Understand atomistic details of conformational changes, protein unfolding, interactions between proteins and drugs
- Study thermodynamics properties (free energies, binding energies)

## Force Fields

A force field (FF) is a set of empirical energy functions and parameters used to calculate the potential energy *U* of a system of atoms and/or molecules as a function of the molecular coordinates. Classical molecular mechanics (MM) potential energy function used in MD simulations is an empirical function comprised of non-bonded and bonded interactions:

$U(\vec{r})=\sum{U_{bonded}}(\vec{r})+\sum{U_{nonbonded}}(\vec{r})$

Typically MD simulations are confined to evaluating only interactions between pairs of atoms. In this approximation force fields are based on two-body potentials, and the energy of the whole system is described by the 2-dimensional force matrix.

For convenience force fields can be divided into 3 general classes based on how complex they are. In the class 1 force field dynamics of bond stretching and angle bending are described by simple harmonic motion, i.e. the magnitude of restoring force is assumed to be proportional to the displacement from the equilibrium position. As the energy of a harmonic oscillator is proportional to the square of the displacement, this approximation is called quadratic. Anharmonic higher-order energy terms are required for a more accurate description of molecular motions. Correlations between bond stretching and angle bending are omitted in the class 1 force field hence force matrix is diagonal.

Class 2 force fields add anharmonic cubic and/or quartic terms to the potential energy for bonds and angles. Besides, they contain cross-terms describing the coupling between adjacent bonds, angles and dihedrals.
Higher-order terms and cross terms allow for a better description of interactions resulting in a more accurate reproduction of bond and angle vibrations. However much more target data is needed for the determination of these additional parameters. Until recently biomolecular force fields were focused on nonbonded interactions and accurate reproduction of critical torsion potentials.

Class 3 force fields explicitly add special effects of organic chemistry. For example stereoelectronic effects, electronegativity effect, Jahn–Teller effect, etc.


## Energy Terms of Biomolecular Force Fields
Most force fields for biomolecular simulations are minimalistic class 1 forcefields trading off rigor of the physical representation for the ability to simulate large systems for a long period of time.

### Non-bonded Terms
#### The Lennard-Jones potential
The Lennard-Jones (LJ) potential approximates the potential energy of non-elecrostatic interaction between a pair of non-bonding atoms or molecules with a simple mathematical function:

$V_{LJ}(r)=\frac{C12}{r^{12}}-\frac{C6}{r^{6}}$

The $$r^{-12}$$ term approximates the strong Pauli repulsion originating from overlap of electron orbitals, while the $$r^{-6}$$ term describes weaker attractive forces acting between local dynamically induced dipoles in the valence orbitals. While the attractive term is physically realistic (London dispersive forces have $$r^{-6}$$ distance dependence), the repulsive term is a crude approximation of exponentially decaying repulsive interaction. The too steep repulsive part often leads to an overestimation of the pressure in the system.

The LJ potential is commonly expressed in terms of the well depth $$\epsilon$$ (the measure of the strength of the interaction) and the van der Waals radius $$\sigma$$ (the distance at which the intermolecular potential between the two particles is zero).

$V_{LJ}(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$

The LJ coefficients *C* are related to the $$\sigma$$  and the $$\epsilon$$  with the equations:

 $C12=4\epsilon\sigma^{12},C6=4\epsilon\sigma^{6}$

#### The Buckingham potential
The Buckingham potential replaces the repulsive $$r^{-12}$$ term in Lennard-Jones potential by exponential function of distance:

$V_{BU}(r)=Aexp(-Br) -\frac{C6}{r^{6}}$

Exponential function describes electron density more realistically but it is computationally more expensive to calculate. While using Buckingham potential there is a risk of "Buckingham Catastrophe", the condition when at short-range electrostatic attraction artificially overcomes the repulsive barrier and collision between atoms occurs. This can be remedied by the addition of $$r^{-12}$$ term.

#### The Combining rules
To describe all *LJ* interactions in a simulations system the matrix of the pairwise interactions is constructed. The *LJ* interactions between different types of atoms are computed by combining the *LJ* parameters. Different force fields use different combining rules.

**Geometric mean:**

$$C12_{ij}=\sqrt{C12_{ii}\times{C12_{jj}}},C6_{ij}=\sqrt{C6_{ii}\times{C6_{jj}}}$$  (GROMOS)

$$\sigma_{ij}=\sqrt{\sigma_{ii}\times\sigma_{jj}},\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}$$ (OPLS)

**Lorentz–Berthelot:**

$$\sigma_{ij}=\frac{\sigma_{ii}+\sigma_{jj}}{2},\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}$$ (CHARM, AMBER). This combining rule is known to overestimate the well depth

**Waldman–Hagler:**

$$\sigma_{ij}=\left(\frac{\sigma_{ii}^{6}+\sigma_{jj}^{6}}{2}\right)^{\frac{1}{6}}$$ , $$ \epsilon_{ij}=\sqrt{\epsilon_{ij}\epsilon_{jj}}\times\frac{2\sigma_{ii}^3\sigma_{jj}^3}{\sigma_{ii}^6+\sigma_{jj}^6}$$

This combining rule was developed specifically for simulation of noble gases.

**Hybrid** (the Lorentz–Berthelot for H and the Waldman–Hagler for other elements). Implemented in the [AMBER-ii](https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5b07233) force field for perfluoroalkanes, noble gases, and their mixtures with alkanes.

> ## Specifying Combining Rules in GROMACS
>
> Combining rule is specified in the **[nonbond_params]** section of the parameter file **ffnonbonded.itp** (in the column 'func').
>
> Geometric mean is selected by using rules 1 and 3;
>  Lorentz–Berthelot rule is selected using rule 2.
>
> GROMOS force field requires rule 1; OPLS requires rule 3; CHARM and AMBER require rule 2
{: .callout}

> ## Specifying Combining Rules in NAMD
>By default, Lorentz–Berthelot rules are used. Geometric mean can be turned on in the run parameter file:
>> **vdwGeometricSigma**=yes
{: .callout}


#### The electrostatic potential
To describe the elecrostatic interactions in MD the point charges are assigned to the positions of atomic nuclei. The atomic charges are derived using QM methods with the goal to approximate the electrostatic potential around a molecule. The electrostatic potential is described with the Coulomb's law:

$V_{Elec}=\frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon_{r}r_{ij}}$

where *r<sub>ij</sub>* is the distance between the pair of atoms, *q<sub>i</sub>* and *q<sub>j</sub>* are the charges on the atoms *i* and *j*,$$\epsilon_{0}$$ is the permittivity of vacuum. and $$\epsilon_{r}$$ is the relative permittivity.

> ## Short-range and Long-range Interactions
> The interaction is termed short-range if the potential decreases faster than *r<sup>-d</sup>*, where *r* is the distance between 2 particles and *d* is dimensionality. Otherwise the interaction is long-ranged. Accordingly the Lennard-Jones interactions are short-ranged, the Coulomb interactions are long-ranged.
{: .callout}

### Bonded Terms

#### The bond potential
The bond potential is used to model the interaction of covalently bonded atoms in a molecule. Bond stretch is approximated by a simple harmonic function describing oscillation about an equilibrium bond length *r<sub>0</sub>* with bond constant *k<sub>b</sub>*:

$V_{Bond}=k_b(r_{ij}-r_0)^2$

This is a fairly poor approximation at extreme stretching, but bonds are so stiff that it works for well moderate temperatures. A Morse potential is more accurate, but more expensive to calculate.

#### The angle potential
The angle potential describes the bond bending energy. It is defined for every triplet of bonded atoms. It is also approximated by a harmonic function describing oscillation about an equilibrium angle  $$\theta_{0}$$  with force constant $$k_\theta$$ :

$V_{Angle}=k_\theta(\theta_{ijk}-\theta_0)^2$

The force constants for angle potential are about 5 times smaller that for bond stretching.

#### The torsion (dihedral) angle potential
The torsion energy is defined for every 4 sequentially bonded atoms. The torsion angle $$\phi$$ is the angle of rotation about the covalent bond between the middle two atoms and the potential is given by:

$V_{Dihed}=k_\phi(1+cos(n\phi-\delta))$

Where the non-negative integer constant *n* defines periodicity and  $$\delta$$ is the phase shift angle.

#### The improper torsion potential
The improper torsion potentialis defined for a group of 4 bonded atoms where the central atom i is connected to the 3 peripheral atoms j,k, and l. Such group can be seen as a pyramid and the improper torsion potential is related to the distance of the central atom from the base of the pyramid. This potential is used mainly to keep molecular structures planar. As there is only one energy minimum the improper torsion term can be given by a harmonic function:

$V_{Improper}=k_\phi(\phi-\phi_0)^2$

Where the dihedral angle $$\phi$$ is the angle between planes ijk and jkl.

### Coupling Terms
#### The Ureu-Bradley potential
It is known that as a bond angle is decreased, the adjacent bonds stretch to reduce the interaction between the outer atoms of the bonded triplet. This means that there is a coupling between bond length and bond angle. This coupling can be decribed by the Ureu-Bradley potential. The Urey-Bradley term is defined as a (noncovalent) spring between the outer *i* and *k* atoms of a bonded triplet *ijk*. It is approximated by a harmonic function describing oscillation about an equilibrium distance *r<sub>ub</sub>* with force constant *k<sub>ub</sub>*:

$V_{UB}=k_{ub}(r_{ik}-r_{ub})^2$

U-B terms are used to improve agreement with vibrational spectra when a harmonic bending term alone would not adequately fit. These phenomena are largely inconsequential for the overall conformational sampling in a typical biomolecular/organic simulation. The Ureu-Bradley term is implemented in the CHARMM force fields.

### The CMAP potential
CMAP is a grid based error correction map designed by CHARMM force field developers to correct for errors in nonbonded interactions, electrostatics, lack of coupling terms, inaccurate combination rules and other force field deficiencies. The grid of energy correction factors is constructed using QM data for every combination of $$\phi/\psi$$ dihedral angles of the peptide backbone and further optimized using empirical data. CMAP potential was initially applied to improve CHARMM22 force field. CMAP corrections were later implemented in AMBER force fields ff99IDPs (force field for intrinsically disordered proteins) and ff12SB-cMAP (force field for implicit-solvent simulations).


## History of Force Field Develoment

### Consistent Force Field (CFF)
CFF was the first modern force field [(Lifson, 1968)]({{ page.root }}/reference.html#Lifson-1968) Introduced a concept of a 'consistent force field'. Introduced a methodology for deriving and validating force fields. The term 'consistent' emphasized importance of the ability to describe a wide range of compounds and physical observables (conformation, crystal structure, thermodynamic properties and vibrational spectra). After the initial derivation for hydrocarbons CFF was extended to proteins, but was very crude at that time.

### Allinger force fields: MM1, MM2, MM3, MM4.
MM1, MM2 - Class 1
MM3 - Class 2
MM4 - Class 3
Target data included electron diffraction, vibrational spectra, heats of formation, and crystal structures. The calculation results were verified via comparison with high-level ab initio quantum chemistry computations, and the parameters were additionally adjusted.

For calculations on small and medium size organic molecules.

### Empirical Conformational Energy Program for Peptides (ECEPP)
ECEPP was the first force field targeting polypeptides and proteins. Crystal data of small organic compounds and semi-empirical QM calculations were used extensively in derivation of this force field. As more experimental data became available, force field parameters have been refined and modified versions ECEPP-2, ECEPP-3 were published.


Over the years of evaluations the problems originating from the deficiency of these FF became apparent and various approaches has been undertaken to correct them.

Force Fields with development focused on application to proteins:

AMBER
CHARM

GROMOS and OPLS are focused on fitting to thermodynamic properties such as the heats of vaporization, liquid densities, and the solvation properties of small molecules.



Most relevant processes require very long simulations. Large systems create computational restraints. The goal is to develop a "minimalistic" force field to expand simulation time window as much as possible.

Force fields with development focused on improving representation of molecular interactions

CVFF
CFF93
CFF95
MM1 - MM4

Development Phase I

Refinement after the initial introduction.
-Converted back to AA, except GROMOS

CHARMM22
AMBER ff99, GAFF
OPLS-AA, OPLS-AA/L

Used large datasets for training. Training datasets were different for different FF.

Large deviations in different observables. Inability to pedict conformations of peptides and proteins.
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
COMPASS and COMPASS II forcefields are developed for simulations of organic molecules, inorganic smallmolecules, and polymers. The VDW parameters are obtained by fitting enthalpy of vaporizations and densities, to experimental data. The atomic partial charges are derived using QM and empirically adjusted to take hydrogen bonding effects into account.

### OPLS
OPLS family force fields are created for liquid simulations containing organic molecules and proteins. The VDW parameters are optimized using experimental liquid properties, mainly enthalpy of vaporizations and densities. The atomic partial charges are derived using QM and experimental condensed-phase properties.

**OPLS-AA**
Jorgensen W, Maxwell D, Tirado-Rives J. Development and testing of the OPLS all-atom force field on conformational energetics and properties of organic liquids. J Am Chem Soc. 1996;118: 11225–11236.

**OPLS-AA/L**
Kaminski GA, Friesner RA, Tirado-Rives J, Jorgensen WL. Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides. J Phys Chem B. 2001;105: 6474–6487. doi:10.1021/jp003919d

**OPLS_2005**
Banks JL, Beard HS, Cao Y, Cho AE, Damm W, Farid R, et al. Integrated Modeling Program, Applied Chemical Theory (IMPACT). J Comput Chem. 2005;26: 1752–1780. doi:10.1002/jcc.20292

**OPLS-AA/M**
Robertson MJ, Tirado-Rives J, Jorgensen WL. Improved Peptide and Protein Torsional Energetics with the OPLS-AA Force Field. J Chem Theory Comput. 2015;11: 3499–3509. doi:10.1021/acs.jctc.5b00356

**OPLS3**
Harder E, Damm W, Maple J, Wu C, Reboul M, Xiang JY, et al. OPLS3: A Force Field Providing Broad Coverage of Drug-like Small Molecules and Proteins. J Chem Theory Comput. 2016;12: 281–296. doi:10.1021/acs.jctc.5b00864


### AMBER
AMBER forcefields are developed for simulations of protein and nucleic acides. The VDW parameters are obtained from crystal structures and lattice energies. The atomic partial charges are derived using QM without any agjustments.

**ff99**
Wang J, Cieplak P, Kollman PA. How well does a restrained electrostatic potential (RESP) model perform in calculating conformational energies of organic and biological molecules? Journal of Computational Chemistry. 2000;21: 1049–1074. doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F

After publication of ff99 a number of studies devoted primarily to modifying the torsion potentials in order to correct the observed discrepancies have been publihed:

**ff99sb*** is optimized for the correct description of the helix-coil equilibrium

**ff99sb-φ'**  targeted the reproduction of the intrinsic conformational preferences of tripeptides

**ff99sb-nmr** and FF99SB_φψ  target data during included protein NMR chemical shifts and residual dipolar couplings.

**ff99sb-ILDN** targeted optimization of four amino acid side chains.

2003 New charges developed B3LYP/cc-pvTZ quantum calculations in a polarizable continuum (PCM) solvent intended to mimic the interior of a protein

**ff14ipq**
- overstabilization of salt bridges.

**ff15ipq** addessed salt bridge overstabilization by increasing radius of polar hydrogens bonded to nitrogen.  The modified FF performed as well or better than the other fixed charge FF. Polarizable CHARMM Drude-2013 and AMOEBA performed better in this respect, as expected.
Problems related to protein stability persist. Even 4 μs simulations “were not sufficiently long to obtain converged estimates of secondary structure of polypeptides”. In simulation tests some proteins deviated significantly near the end of several microsecond simulations, and it is not clear whether this is a transient fluctuation or transition to a different state.

### CHARMM

CHARMM22

CHARMM22/CMAP (CHARMM27)

CHARMM36 refined backbone CMAP potentials and introduced new side-chain dihedral parameters. The updated CMAP corrected the C22/CMAP FF bias towards alpha-helical conformations.


### Polarizable Force fields
**Drude-2013**
**AMOEBA**

The origins of FF based calculations, theory and methodology of FF development have been recently reviewed in [(Dauber-Osguthorpe, 2019)]({{ page.root }}/reference.html#Dauber-Osguthorpe-2019), and the latest developments in improvement of FF accuracy and robustness are discussed in [(Hagler, 2019)]({{ page.root }}/reference.html#Hagler-2019).


## How to create ligand topology

[Automated Topology Builder](https://atb.uq.edu.au/index.py)
