---
title: "Force Fields and Interactions"
teaching: 30
exercises: 10
questions:
- "What is Molecular Dynamics and how can I benefit from using it?"
- "What is a force field?"
- "What energy terms are used in biomolecular force fields mathematical?"
objectives:
- "Understand strengths and weaknesses of MD simulations"
- "Understand how the interactions between particles in a molecular dynamics simulation are modeled"
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

For convenience force fields can be divided into 3 general classes based on how complex they are.

Class 1. In the class 1 force field dynamics of bond stretching and angle bending are described by simple harmonic motion, i.e. the magnitude of restoring force is assumed to be proportional to the displacement from the equilibrium position. As the energy of a harmonic oscillator is proportional to the square of the displacement, this approximation is called quadratic. In general, bond stretching and angle bending are close to harmonic only near the equilibrium. Higher-order anharmonic energy terms are required for a more accurate description of molecular motions. In the class 1 force field force matrix is diagonal because correlations between bond stretching and angle bending are omitted.

Class 2 force fields add anharmonic cubic and/or quartic terms to the potential energy for bonds and angles. Besides, they contain cross-terms describing the coupling between adjacent bonds, angles and dihedrals. Higher-order terms and cross terms allow for a better description of interactions resulting in a more accurate reproduction of bond and angle vibrations. However much more target data is needed for the determination of these additional parameters.

Class 3 force fields explicitly add special effects of organic chemistry. For example stereoelectronic effects, electronegativity effect, Jahn–Teller effect, etc.

Until recently force fields for biomolecular simulations were focused on nonbonded interactions and accurate reproduction of critical torsion potentials.

## Energy Terms of Biomolecular Force Fields
Most force fields for biomolecular simulations are minimalistic class 1 force fields trading off rigor of the physical representation for the ability to simulate large systems for a long period of time.

### Non-Bonded Terms
#### The Lennard-Jones potential
The Lennard-Jones (LJ) potential approximates the potential energy of non-elecrostatic interaction between a pair of non-bonding atoms or molecules with a simple mathematical function:

$V_{LJ}(r)=\frac{C12}{r^{12}}-\frac{C6}{r^{6}}$

The $$r^{-12}$$ term approximates the strong Pauli repulsion originating from overlap of electron orbitals, while the $$r^{-6}$$ term describes weaker attractive forces acting between local dynamically induced dipoles in the valence orbitals. While the attractive term is physically realistic (London dispersive forces have $$r^{-6}$$ distance dependence), the repulsive term is a crude approximation of exponentially decaying repulsive interaction. The too steep repulsive part often leads to an overestimation of the pressure in the system.

The LJ potential is commonly expressed in terms of the well depth $$\epsilon$$ (the measure of the strength of the interaction) and the van der Waals radius $$\sigma$$ (the distance at which the intermolecular potential between the two particles is zero).

$V_{LJ}(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$

The LJ coefficients *C* are related to the $$\sigma$$  and the $$\epsilon$$  with the equations:

 $C12=4\epsilon\sigma^{12},C6=4\epsilon\sigma^{6}$


#### The Lennard-Jones Combining Rules
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


#### The Buckingham potential
The Buckingham potential replaces the repulsive $$r^{-12}$$ term in Lennard-Jones potential by exponential function of distance:

$V_{B}(r)=Aexp(-Br) -\frac{C}{r^{6}}$

Exponential function describes electron density more realistically but it is computationally more expensive to calculate. While using Buckingham potential there is a risk of "Buckingham Catastrophe", the condition when at short-range electrostatic attraction artificially overcomes the repulsive barrier and collision between atoms occurs. This can be remedied by the addition of $$r^{-12}$$ term.

There is only one combining rule for Buckingham potential in GROMACS:

$A_{ij}=\sqrt{(A_{ii}A_{jj})}$

$B_{ij}=2/(\frac{1}{B_{ii}}+\frac{1}{B_{jj}})$

$C_{ij}=\sqrt{(C_{ii}C_{jj})}$


> ## Specifying Combining Rules
> **GROMACS**
>
> Combining rule is specified in the **[defaults]** section of the **forcefield.itp** file (in the column 'comb-rule').
> ~~~
>[ defaults ]
>; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
>1               2               yes             0.5     0.8333
> ~~~
> {: .source}
> Geometric mean is selected by using rules 1 and 3;
>  Lorentz–Berthelot rule is selected using rule 2.
>
> GROMOS force field requires rule 1; OPLS requires rule 3; CHARM and AMBER require rule 2
>
> The type of potential function is specified in the 'nbfunc' column: 1 selects Lennard-Jones potential, 2 selects Buckingham potential.
{: .callout}
>
> **NAMD**
>
>By default, Lorentz–Berthelot rules are used. Geometric mean can be turned on in the run parameter file:
>~~~
>vdwGeometricSigma yes
>~~~
> {: .source}
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

### CHARMM CMAP potential
CMAP is a grid based error correction map used in CHARMM force field to  correct for errors in nonbonded interactions, electrostatics, lack of coupling terms, inaccurate combination rules and other force field deficiencies. The grid of energy correction factors is constructed using QM data for every combination of $$\phi/\psi$$ dihedral angles of the peptide backbone and further optimized using empirical data. CMAP potential was initially applied to improve CHARMM22 force field. CMAP corrections were later implemented in AMBER force fields [ff99IDPs](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00043) (force field for intrinsically disordered proteins) and [ff12SB-cMAP](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00662) (force field for implicit-solvent simulations).

### Exclusions from Non-Bonded Interactions
Pairs of atoms connected by chemical bonds are normally excluded from computation of non-bonded interactions because bonded energy terms replace non-bonded interactions. In biomolecular force fields all pairs of connected atoms separated by up to 2 bonds (1-2 and 1-3 pairs) are excluded from non-bonded interactions. Computation of the non-bonded interaction between 1-4 pairs depends on the specific force field. Some force fields exclude VDW interactions and scale down electrostatic (AMBER) while others may modify both or use electrostatic as is.

> ## Specifying Exclusions
> **GROMACS**
>
> The exclusions are generated by **grompp** as specified in the **[moleculetype]** section of the molecular topology **.top** file:
> ~~~
>[ moleculetype ]
>; name  nrexcl
>Urea         3
> ...
>[ exclusions ]
>;  ai    aj
>    1     2
>    1     3
>    1     4
> ~~~
>{: .source}
> In the example above non-bonded interactions between atoms that are no farther than 3 bonds are excluded (nrexcl=3). Extra exclusions may be added explicitly in the **[exclusions]** section.
>
> The scaling factors for 1-4 pairs, **fudgeLJ** and **fudgeQQ**, are specified in the **[defaults]** section of the **forcefield.itp** file. While **fudgeLJ** is used only when **gen-pairs** is set to 'yes', **fudgeQQ** is always used.
>
> ~~~
>[ defaults ]
>; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
>1               2               yes             0.5     0.8333
> ~~~
> {: .source}
{: .callout}
> **NAMD**
>
> Which pairs of bonded atoms should be excluded is specified by the **exclude** parameter.<br/> Acceptable values: **none, 1-2, 1-3, 1-4,** or **scaled1-4**
> ~~~
>exclude scaled1-4
>1-4scaling 0.83
>scnb 2.0
> ~~~
> {: .source}
> If **scaled1-4**  is chosen, the electrostatic interactions for 1-4 pairs are multiplied by a constant factor specified by the **1-4scaling** parameter. The LJ interacions for 1-4 pairs are divided by **scnb**.
{: .callout}

> ## Counting Non-Bonded Interactions
>
> 1. How many non-bonded interactions are in the system comprised of 10 Argon atoms?
>
> 2. How many VDW interactions are in the system comprised of 2 propane molecules?
>
> > ## Solution
> >
> > 1. Argon atoms are neutral, so there is no Coulomb interaction. Atoms don't interact with themselves and the interaction ij is the same as the interation ji.  Thus the total number of pairwise non-bonded interactions is (10x10 - 10)/2 = 45.
> >
> > 2. Propane molecule has 11 atoms. Each atom in one molecule interacts with each atom in another molecule, so the number of intermolecular VDW interactions is 11x11=121. In propane molecule only terminal hydrogens are separated by more than 3 bonds. So there are 3x3=9 intramolecular VDW interactions in each of the molecules. The total number of VDW interactions is 121+18=139
> >
> {: .solution}
{: .challenge}
