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

The potential energy function *U* acts as a cornerstone of the MD simulations because it allows calculating the forces. The force on an object is the negative of the derivative of the potential energy function:

$\vec{F}=-\nabla{U}(\vec{r})$

- Molecular dynamics helps to link physics, chemistry and biology
- With the help of MD it is possible to model phenomena that cannot be studied experimentally.
-  Understand atomistic details of conformational changes, protein unfolding, interactions between proteins and drugs
- Study thermodynamics properties (free energies, binding energies)

## Force Fields
A force field (FF) is a set of empirical energy functions and parameters used to calculate the potential energy *U* of a system of atoms and/or molecules in molecular dynamics simulations. The potential energy includes non-bonded and bonded interactions:

$U(\vec{r})=\sum{U_{bonded}}(\vec{r})+\sum{U_{nonbonded}}(\vec{r})$

The origins of FF based calculations, and the theory and methodology of FF development is reviewed in [(Dauber-Osguthorpe, 2019)]({{ page.root }}/reference.html#Dauber-Osguthorpe-2019), and the latest developments in improvement of FF rigor and robustness are discussed in [(Hagler, 2019)]({{ page.root }}/reference.html#Hagler-2019).

## Non-bonded interactions
### The Lennard-Jones potential
The Lennard-Jones (LJ) potential approximates the potential energy of non-elecrostatic interaction between a pair of non-bonding atoms or molecules with a simple mathematical function:

$V_{LJ}(r)=\frac{C12}{r^{12}}+\frac{C6}{r^{6}}$

The $$r^{-12}$$ term approximates the strong Pauli repulsion originating from overlap of electron orbitals, while the $$r^{-6}$$ term describes weaker attractive van der Waals forces acting between local dynamically induced dipoles in the valence orbitals.

The LJ potential is commonly expressed in terms of the well depth $$\epsilon$$ (the measure of the strength of the interaction) and the van der Waals radius $$\sigma$$ (the distance at which the intermolecular potential between the two particles is zero).

$V_{LJ}(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$

The LJ coefficients *C* are related to the $$\sigma$$  and the $$\epsilon$$  with the equations:

 $C12=4\epsilon\sigma^{12},C6=4\epsilon\sigma^{6}$

To describe all *LJ* interactions in a simulations system the matrix of the pairwise interactions is constructed. The *LJ* interactions between different types of atoms are computed by combining the *LJ* parameters. Different force fields use different combining rules.

### The Combining rules

**Geometric mean:**

$$C12_{ij}=\sqrt{C12_{ii}\times{C12_{jj}}},C6_{ij}=\sqrt{C6_{ii}\times{C6_{jj}}}$$  (GROMOS)

$$\sigma_{ij}=\sqrt{\sigma_{ii}\times\sigma_{jj}},\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}$$ (OPLS)

**Lorentz–Berthelot:**

$$\sigma_{ij}=\frac{\sigma_{ii}+\sigma_{jj}}{2},\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}$$ (CHARM, AMBER)

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


### The electrostatic potential
To describe the elecrostatic interactions in MD the point charges are assigned to the positions of atomic nuclei. The atomic charges are derived using QM methods with the goal to approximate the electrostatic potential around a molecule. The electrostatic potential is described with the Coulomb's law:

$V_{Elec}=\frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon_{r}r_{ij}}$

where *r<sub>ij</sub>* is the distance between the pair of atoms, *q<sub>i</sub>* and *q<sub>j</sub>* are the charges on the atoms *i* and *j*,$$\epsilon_{0}$$ is the permittivity of vacuum. and $$\epsilon_{r}$$ is the relative permittivity.

> ## Short-range and Long-range Interactions
> The interaction is termed short-range if the potential decreases faster than *r<sup>-d</sup>*, where *r* is the distance between 2 particles and *d* is dimensionality. Otherwise the interaction is long-ranged. Accordingly the Lennard-Jones interactions are short-ranged, the Coulomb interactions are long-ranged.
{: .callout}
## Bonded Interactions

### The bond potential
The bond potential is used to model the interaction of covalently bonded atoms in a molecule. Bond stretch is approximated by a simple harmonic function describing oscillation about an equilibrium bond length *r<sub>0</sub>* with bond constant *k<sub>b</sub>*:

$V_{Bond}=k_b(r_{ij}-r_0)^2$

This is a fairly poor approximation at extreme stretching, but bonds are so stiff that it works for well moderate temperatures. A Morse potential is more accurate, but more expensive to calculate.

### The angle potential
The angle potential describes the bond bending energy. It is defined for every triplet of bonded atoms. It is also approximated by a harmonic function describing oscillation about an equilibrium angle  $$\theta_{0}$$  with force constant $$k_\theta$$ :

$V_{Angle}=k_\theta(\theta_{ijk}-\theta_0)^2$

The force constants for angle potential are about 5 times smaller that for bond stretching.

### The torsion (dihedral) angle potential
The torsion energy is defined for every 4 bonded atoms.
The torsion angle $$\phi$$ is the angle between 2 planes defined by the first and the last 3 atoms of the 4 atoms involved in the torsion interaction.

$V_{Dihed}=k_\phi(1+cos(n\phi-\delta))$

 The non-negative integer constant *n* defines periodicity and  $$\delta$$ is the phase shift angle.

> ## The Ureu-Bradley potential
> The Ureu-Bradley potential describes cross-terms (correlation between bond length and bond angle). The presence of cross-terms in a force field reflects couplings between the internal coordinates. As a bond angle is decreased, it is found that the adjacent bonds stretch to reduce the interaction between the 1,3 atoms.
>
>The Urey-Bradley term is defined as a (noncovalent) spring between the outer *i* and *k* atoms of a bonded triplet *ijk*. It is approximated by a harmonic function describing oscillation about an equilibrium distance *r<sub>ub</sub>* with force constant *k<sub>ub</sub>*:
>
>$V_{UB}=k_{ub}(r_{ik}-r_{ub})^2$
>
> U-B terms have been used to improve agreement with vibrational spectra when a harmonic bending term alone would not adequately fit. These phenomena are largely inconsequential for the overall conformational sampling in a typical biomolecular/organic simulation.
{: .callout}
