---
title: "Introducing molecular dynamics simulation"
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

# Basic principles of molecular dynamics
Molecular dynamics (MD) simulations are widely used as tools to investigate structure and dynamics of proteins, nucleic acids, carbohydrates, lipids, nanoparticles and liquid/solid interfaces under a wide variety of conditions. MD is the simulation of atomic positions in time accomplished by solving classical Newton's equation of motion stating that the rate of change of momentum <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{p}" title="\vec{p}" /></a> of an object equals the force <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{F}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{F}" title="\vec{F}" /></a> acting on it:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=\frac{d\vec{p}}{dt}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=\frac{d\vec{p}}{dt}" title="\vec{F}=\frac{d\vec{p}}{dt}" /></a>

The potential energy function *U* acts as a cornerstone of the MD simulations because it allows calculating the forces. The force on an object is the negative of the derivative of the potential energy function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=-\nabla&space;U(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=-\nabla&space;U(\vec{r})" title="\vec{F}=-\nabla U(\vec{r})" /></a>

- Molecular dynamics helps to link physics, chemistry and biology
- With the help of MD it is possible to model phenomena that cannot be studied experimentally.
-  Understand atomistic details of conformational changes, protein unfolding, interactions between proteins and drugs
- Study thermodynamics properties (free energies, binding energies)


## Force Fields
A force field is a set of empirical energy functions and parameters used to calculate the potential energy *U* of a system of atoms and/or molecules in molecular dynamics simulations. The potential energy includes non-bonded and bonded interactions.

<a href="https://www.codecogs.com/eqnedit.php?latex=U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" title="U(\vec{r})=\sum U_{bonded}(\vec{r})+\sum U_{nonbonded}(\vec{r})" /></a>


### Non-bonded interactions
#### The Lennard-Jones potential
The Lennard-Jones (LJ) potential approximates the potential energy of non-elecrostatic interaction between a pair of non-bonding atoms or molecules with a simple mathematical function:

<img src="https://latex.codecogs.com/gif.latex?V_{LJ}(r)=\frac{C12}{r^{12}}&plus;\frac{C6}{r^{6}}" />

The <img src="https://latex.codecogs.com/gif.latex?r^{-12}"/> term approximates the strong Pauli repulsion originating from overlap of electron orbitals, while the <img src="https://latex.codecogs.com/gif.latex?r^{-6}"/> term describes weaker attractive van der Waals forces acting between local dynamically induced dipoles in the valence orbitals.

The LJ potential is commonly expressed in terms of the well depth <img src="https://latex.codecogs.com/gif.latex?\epsilon" /> (the measure of the strength of the interaction) and the van der Waals radius <img src="https://latex.codecogs.com/gif.latex?\sigma" /> (the distance at which the intermolecular potential between the two particles is zero).

<img src="https://latex.codecogs.com/gif.latex?V_{LJ}(r)=4\epsilon\left&space;[&space;\left&space;(&space;\frac{\sigma}{r}\right&space;)^{12}-&space;\left&space;(&space;\frac{\sigma}{r}\right&space;)^{6}&space;\right&space;]"/>

The LJ coefficients *C* are related to the <img src="https://latex.codecogs.com/gif.latex?\sigma"/>  and the <img src="https://latex.codecogs.com/gif.latex?\epsilon"/>  with the equations:

 <img src="https://latex.codecogs.com/gif.latex?&C12=4\epsilon\sigma^{12},C6=4\epsilon\sigma^{6}"/>

To describe all *LJ* interactions in a simulations system the matrix of the pairwise interactions is constructed. The *LJ* interactions between different types of atoms are computed by combining the *LJ* parameters. Different force fields use different combining rules.

#### The Combining rules

- **Geometric mean:**<br>
<img src="https://latex.codecogs.com/gif.latex?C12_{ij}=\sqrt{C12_{ii}\times{C12_{jj}}},&space;C6_{ij}=\sqrt{C6_{ii}\times{C6_{jj}}}"/>  , (GROMOS)<br><img src="https://latex.codecogs.com/gif.latex?\sigma_{ij}=\sqrt{\sigma_{ii}\times\sigma_{jj}},&space;\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}" />, (OPLS)<br><br>

- **Lorentz–Berthelot:**<br>
<img src="https://latex.codecogs.com/gif.latex?\sigma_{ij}=\frac{\sigma_{ii}&plus;\sigma_{jj}}{2},&space;\epsilon_{ij}=\sqrt{\epsilon_{ii}\times\epsilon_{jj}}"  />, (CHARM, AMBER)

- **Waldman–Hagler:**<br>
<img src="https://latex.codecogs.com/gif.latex?\sigma_{ij}=\left&space;(&space;\frac{&space;\sigma_{ii}^{6}&space;&plus;&space;\sigma_{jj}^{6}}{2}&space;\right&space;)^{\frac{1}{6}}"/>  <img src="https://latex.codecogs.com/gif.latex?\epsilon_{ij}=\sqrt{\epsilon_{ij}&space;\epsilon_{jj}}\times\frac{2\sigma_{ii}^3&space;\sigma_{jj}^3}{\sigma_{ii}^6&space;&plus;\sigma_{jj}^6&space;}"/><br>This combining rule was developed specifically for simulation of noble gases<br><br>

- **Hybrid** (the Lorentz–Berthelot for H and the Waldman–Hagler for other elements)
  Implemented in the [AMBER-ii](https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5b07233) force field for perfluoroalkanes, noble gases, and their mixtures with alkanes.


> ## Specifying Combining Rules
>
>  **GROMACS**
>
> In GROMACS combining rules are specified in the **[nonbond_params]** section of the parameter file **ffnonbonded.itp** (in the column 'func').
>
> Geometric mean is selected by using rules 1 and 3,  Lorentz–Berthelot rule is selected using rule 2.
>
> GROMOS requires rule 1, OPLS requires rule 3, CHARM and AMBER require rule 2
>
>  **NAMD**
>
> By default, NAMD uses Lorentz–Berthelot rules. Geometric mean can be turned on in the run parameter file:
> 
> **vdwGeometricSigma**=yes
{: .callout}


#### The electrostatic potential
To describe the elecrostatic interactions in MD the point charges are assigned to the positions of atomic nuclei. The atomic charges are derived using QM methods with the goal to approximate the electrostatic potential around a molecule. The electrostatic potential is described with the Coulomb's law:

<img src="https://latex.codecogs.com/gif.latex?V_{Elec}=\frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon_{r}&space;r_{ij}}" /><br>
where *r<sub>ij</sub>* is the distance between the pair of atoms, *q<sub>i</sub>* and *q<sub>j</sub>* are the charges on the atoms *i* and *j*, <img src="https://latex.codecogs.com/gif.latex?\epsilon_{0}"/> is the permittivity of vacuum. and  <img src="https://latex.codecogs.com/gif.latex?\epsilon_{r}"/> is the relative permittivity.
### Bonded Interactions

#### The bond potential
The bond potential is used to model the interaction of covalently bonded atoms in a molecule. Bond stretch is approximated by a simple harmonic function describing oscillation about an equilibrium bond length *r<sub>0</sub>* with bond constant *k<sub>b</sub>*:

<img src="https://latex.codecogs.com/gif.latex?V_{Bond}=k_b(r_{ij}-r_0)^2" /><br>

This is a fairly poor approximation at extreme stretching, but bonds are so stiff that it works for well moderate temperatures. A Morse potential is more accurate, but more expensive to calculate.

#### The angle potential
The angle potential describes the bond bending energy. It is defined for every triplet of bonded atoms. It is also approximated by a harmonic function describing oscillation about an equilibrium angle  <img src="https://latex.codecogs.com/gif.latex?\theta_{0}"/>  with force constant  <img src="https://latex.codecogs.com/gif.latex?k_\theta"/> :

<img src="https://latex.codecogs.com/gif.latex?V_{Angle}=k_\theta(\theta_{ijk}-\theta_0)^2" /><br>

The force constants for angle potential are about 5 times smaller that for bond stretching.

#### The torsion (dihedral) angle potential
The torsion energy is defined for every 4 bonded atoms.
The torsion angle <img src="https://latex.codecogs.com/gif.latex?\phi"/> is the angle between 2 planes defined by the first and the last 3 atoms of the 4 atoms involved in the torsion interaction.

<img src="https://latex.codecogs.com/gif.latex?V_{Dihed}=k_\phi(1+cos(n\phi-\delta))" /><br>

 The non-negative integer constant *n* defines periodicity and  <img src="https://latex.codecogs.com/gif.latex?\delta"/> is the phase shift angle.

#### The Ureu-Bradley potential

The Ureu-Bradley potential describes cross-terms (correlation between bond length and bond angle). The presence of cross-terms in a force field reflects couplings between the internal coordinates. As a bond angle is decreased, it is found that the adjacent bonds stretch to reduce the interaction between the 1,3 atoms.

The Urey-Bradley term is defined as a (noncovalent) spring between the outer *i* and *k* atoms of a bonded triplet *ijk*. It is approximated by a harmonic function describing oscillation about an equilibrium distance *r<sub>ub</sub>* with force constant *k<sub>ub</sub>*:

<img src="https://latex.codecogs.com/gif.latex?V_{UB}=k_{ub}(r_{ik}-r_{ub})^2" /><br>

 U-B terms have been used to improve agreement with vibrational spectra when a harmonic bending term alone would not adequately fit. These phenomena are largely inconsequential for the overall conformational sampling in a typical biomolecular/organic simulation.
