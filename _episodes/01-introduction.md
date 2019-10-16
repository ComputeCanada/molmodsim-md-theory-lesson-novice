---
title: "Force Fields"
teaching: 30
exercises: 0
questions:
- "What is a force field"
objectives:
- "Understand the interactions between particles in a molecular dynamics simulation"
keypoints:
- "First key point."
---

<!-- MDTOC maxdepth:6 firsth1:1 numbering:1 flatten:0 bullets:1 updateOnSave:1 -->

- 1. [Basic principles of molecular dynamics](#basic-principles-of-molecular-dynamics)
   - 1.1. [Force Fields](#force-fields)
      - 1.1.1. [Non-bonded interactions](#non-bonded-interactions)
         - 1.1.1.1. [The Lennard-Jones potential](#the-lennard-jones-potential)
         - 1.1.1.2. [The Combining rules](#the-combining-rules)
         - 1.1.1.3. [Specifying Combining Rules](#specifying-combining-rules)
            - 1.1.1.3.1. [GROMACS](#gromacs)
            - 1.1.1.3.2. [NAMD](#namd)
         - 1.1.1.4. [The electrostatic potential](#the-electrostatic-potential)
      - 1.1.2. [Bonded Interactions](#bonded-interactions)
         - 1.1.2.1. [The bond potential](#the-bond-potential)
         - 1.1.2.2. [The angle potential](#the-angle-potential)
         - 1.1.2.3. [The torsion angle potential](#the-torsion-angle-potential)
         - 1.1.2.4. [The Ureu-Bradley potential](#the-ureu-bradley-potential)
   - 1.2. [Truncation of Interactions](#truncation-of-interactions)
      - 1.2.1. [Neighbour Searching Methods](#neighbour-searching-methods)
         - 1.2.1.1. [Cell Lists](#cell-lists)
         - 1.2.1.2. [Verlet Lists](#verlet-lists)
      - 1.2.2. [Truncation of Lennard-Jones Interactions](#truncation-of-lennard-jones-interactions)
         - 1.2.2.1. [Cutoffs Used for Development of Common Force Fields](#cutoffs-used-for-development-of-common-force-fields)
      - 1.2.3. [Specifying Truncation of LJ Potential](#specifying-truncation-of-lj-potential)
         - 1.2.3.1. [GROMACS](#gromacs)
         - 1.2.3.2. [NAMD](#namd)
      - 1.2.4. [Truncation of the Electrostatic Interactions](#truncation-of-the-electrostatic-interactions)
      - 1.2.5. [Specifying Cutoff and Neighbour Searching Methods](#specifying-cutoff-and-neighbour-searching-methods)
         - 1.2.5.1. [GROMACS](#gromacs)
         - 1.2.5.2. [NAMD](#namd)
   - 1.3. [Balancing of charges](#balancing-of-charges)
   - 1.4. [Periodic boundary conditions](#periodic-boundary-conditions)
      - 1.4.1. [Specifying periodic box](#specifying-periodic-box)
         - 1.4.1.1. [GROMACS](#gromacs)
         - 1.4.1.2. [NAMD](#namd)
   - 1.5. [Integrating the Equations of Motion.](#integrating-the-equations-of-motion)
      - 1.5.1. [Integration Algorithms](#integration-algorithms)
         - 1.5.1.1. [The Euler Algorithm](#the-euler-algorithm)
         - 1.5.1.2. [The Verlet Algorithm](#the-verlet-algorithm)
         - 1.5.1.3. [The Velocity Verlet Algorithm](#the-velocity-verlet-algorithm)
         - 1.5.1.4. [The Leap Frog Algorithm](#the-leap-frog-algorithm)
      - 1.5.2. [Choosing Time Step](#choosing-time-step)
      - 1.5.3. [Specifying Time Step Parameters](#specifying-time-step-parameters)
      - 1.5.4. [Specifying Integration Method](#specifying-integration-method)
         - 1.5.4.1. [GROMACS](#gromacs)
         - 1.5.4.2. [NAMD](#namd)
   - 1.6. [Molecular Dynamics Software Available on Compute Canada Systems](#molecular-dynamics-software-available-on-compute-canada-systems)
      - 1.6.1. [AMBER](#amber)
      - 1.6.2. [DL_POLY](#dl_poly)
      - 1.6.3. [GROMACS](#gromacs)
         - 1.6.3.1. [Force fields implemented in GROMACS:](#force-fields-implemented-in-gromacs)
      - 1.6.4. [NAMD](#namd)
         - 1.6.4.1. [Force fields implemented in NAMD:](#force-fields-implemented-in-namd)
      - 1.6.5. [LAMMPS](#lammps)
   - 1.7. [Water Models](#water-models)

<!-- /MDTOC -->

[Link to the GitLab Repository](https://git.computecanada.ca/svassili/bst-md-theory-lesson-novice/blob/gh-pages/_episodes/01-introduction.md)

[Link to the GitHub Repository](https://github.com/ssvassiliev/bst-md-theory-lesson-novice/blob/gh-pages/_episodes/01-introduction.md#the-lennard-jones-potential)

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

#### Specifying Combining Rules
##### GROMACS
>In GROMACS combining rules are specified in the **[nonbond_params]** section of the parameter file **ffnonbonded.itp** (in the column 'func').
>
>Geometric mean is selected by using rules 1 and 3,  Lorentz–Berthelot rule is selected using rule 2.
>
>GROMOS requires rule 1,
>OPLS requires rule 3,
>CHARM and AMBER require rule 2
##### NAMD
>By default, NAMD uses Lorentz–Berthelot rules. Geometric mean can be turned on in the run parameter file:
>**vdwGeometricSigma**=yes

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

The presence of cross-terms in a force field reflects couplings between the internal coordinates. As a bond angle is decreased, it is found that the adjacent bonds stretch to reduce the interaction between the 1,3 atoms.

The Urey-Bradley term is defined as a (noncovalent) spring between the outer *i* and *k* atoms of a bonded triplet *ijk*. It is approximated by a harmonic function describing oscillation about an equilibrium distance *r<sub>ub</sub>* with force constant *k<sub>ub</sub>*:

<img src="https://latex.codecogs.com/gif.latex?V_{UB}=k_{ub}(r_{ik}-r_{ub})^2" /><br>



 U-B terms have been used to improve agreement with vibrational spectra when a harmonic bending term alone would not adequately fit. These phenomena are largely inconsequential for the overall conformational sampling in a typical biomolecular/organic simulation.

## Truncation of Interactions
The most computationally demanding part of a molecular dynamics simulation is the calculation of the nonbonded terms of the potential energy function. As non-bonded energy terms between every pair of atoms should be evaluated, the number of calculations increases as the square of the number of atoms. To speed up the computation, only the interactions between two atoms separated by a distance less than a pre-defined cutoff distance are evaluated.

### Neighbour Searching Methods
 The search for pairs of particles that are needed for calculation of the short-range nonbonded interactions is usually accelerated by maintaining a list of all particles within a predefined cutoff distance of each other.  Particle neighbours are determined either by dividing the simulation system into grid cells (cell lists) or by constructing a neighbour list for each particle (Verlet lists).

#### Cell Lists
The cell lists method divides the simulation domain into *n* cells within edge length greater or equal to the cutoff radius of the interaction to be computed.  The interaction potential for each particle is then computed as the sum of the pairwise interactions between the particle and all other particles in the same cell and all other particles in the neighbouring cells (26 cells for 3-dimensional simulation).

#### Verlet Lists
A Verlet list stores all particles within the cutoff distance of every particle plus some extra buffer distance. Although all pairwise distances must be evaluated to construct the Verlet list, it can be used for several consecutive time steps until any particle has moved more than half of the buffer distance. At this point the list is invalidated and the new list must be constructed. Verlet offer more efficient computation of pairwise interactions at the expence of relatively large memory requirement which can be a limiting factor. In practice, almost all simulations are run in parallel and use a combination of spatial decomposition and Verlet lists.

### Truncation of Lennard-Jones Interactions

The LJ potential is always truncated at the cutoff distance. How to choose the appropriate cutoff distance? A common practice is to truncate at <img src="https://latex.codecogs.com/gif.latex?2.5\sigma"/> and this practice has become a minimum standard for truncation.  At this distance, the LJ potential is about 1/60 of the well depth <img src="https://latex.codecogs.com/gif.latex?\epsilon"/> and it is assumed that errors arising from this truncation are small enough. The dependence of the cutoff on <img src="https://latex.codecogs.com/gif.latex?\sigma"/> means that the choice of the cutoff distance depends on the force field and atom types used in the simulation. For example for the O, N, C, S, and P atoms in the AMBER99 force field the values of <img src="https://latex.codecogs.com/gif.latex?\sigma"/> are in the range 1.7-2.1,  while for the Cs ions  <img src="https://latex.codecogs.com/gif.latex?\sigma=3.4"/>. Thus the minimum acceptable cutoff, in this case, is 8.5.

In practice, increasing cutoff does not necessarily improve accuracy. Each force field has been developed using a certain cutoff value, and effects of the truncation were compensated by adjustment of some other parameters. If you use cutoff 14 for the force field developed with the cutoff 9, then you cannot say that you used this forcefield. Thus to ensure consistency and reproducibility of simulation you should choose the cutoff appropriate for the force field.

#### Cutoffs Used for Development of Common Force Fields
> AMBER: 9,
> CHARMM: 12,
> GROMOS: 14,
> OPLS: 11-15 (depending on a molecule size)

There are several different ways to truncate the non-bonded interaction. The main option to control how LJ potential is truncated is the switching parameter. If the switching is turned on, the smooth switching function is applied to truncate the Lennard-Jones potential smoothly at the cutoff distance. If the switching function is applied the switching distance parameter specifies the distance at which the switching function starts to modify the LJ potential to bring it to zero at the cutoff distance.

### Specifying Truncation of LJ Potential
#### GROMACS
**vdw-modifier**
Acceptable values:
> **potential-shift**: shifts the Van der Waals potential by a constant such that it is zero at the **rvdw**.
>
> **force-switch**: smoothly switches the forces to zero between **rvdw-switch** and **rvdw**.
>
> **potential-switch**: smoothly switches the potential to zero between **rvdw-switch** and **rvdw**.
>
>**none**

**rvdw-switch**
> Where to start switching

**rvdw**
> Cut-off distance


#### NAMD
**cutoff**
> Cut-off distance common to both electrostatic and van der Waals calculations

**switching**
> Turn switching on/off

**switchdist**
> Where to start switching

**vdwForceSwitching**
> use force switching for VDW
> Default Value: off


### Truncation of the Electrostatic Interactions
Electrostatic interactions occurring over long distances are known to be important for biological molecules. Electrostatic interactions decay slowly and simple increase of the cutoff distance to account for long-range interactions can dramatically raise the computational cost. In periodic simulation systems, the most commonly used method for calculation of long-range electrostatic interactions is particle-mesh Ewald.  In this method, the electrostatic interaction is divided into two parts: a short-range contribution, and a long-range contribution. The short-range contribution is calculated by exact summation of all pairwise interactions of atoms separated by a distance that is less than cutoff in real space. The forces beyond the cutoff radius are approximated in Fourier space commonly by the Particle-Mesh Ewald (PME) method.

### Specifying Cutoff and Neighbour Searching Methods

#### GROMACS ####

**cutoff-scheme**
> Since version 5.1 **group** list has been deprecated and only **Verlet** scheme is available

**rlist**
> Cutoff distance for the short-range neighbour list. Active when **verlet-buffer-tolerance** = -1, otherwise ignored

**verlet-buffer-tolerance**
> The maximum allowed error for pair interactions per particle caused by the Verlet buffer. To achieve the predefined tolerance the cutoff distance **rlist** is adjusted indirectly. To override this feature set the value to -1
>
> Default value: 0.005 [kJ mol<sup>-1</sup> ps<sup>-1</sup>]

**nstlist**
> Frequency to update the neighbour list. If set to 0 the neighbour list is constructed only once and never updated.
>
> Default Value: 10

**ns-type**
> Neighbour search method.
> Acceptable values:
>> **grid**: make a grid in the box and only check atoms in neighboring grid cells.<br>
>> **simple**: loop over every atom in the box.

#### NAMD
 When run in parallel NAMD uses a combination of spatial decomposition into grid cells, "patches" and Verlet lists with extended cutoff distance

**cutoff**
> Local interaction distance. Same for both electrostatic and VDW interactions

**pairlistdist**
> Distance between pairs for inclusion in pair lists. Should be bigger or equal than the **cutoff**.
>
> Default value: **cutoff**
>
**stepspercycle**
> Number of timesteps in each cycle. Each cycle represents the number of timesteps between atom reassignments.
> Default Value: 20

**pairlistsPerCycle**
> Specifies how many times per cycle to regenerate pairlists.
>
> Default value: 2


## Balancing of charges
Neutralizing a system is a practice carried out for obtaining correct electrostatic energy during the simulation. This is done because under periodic boundary and using grid-based electrostatic the system has to be neutral. Otherwise, the electrostatic energy will essentially add to infinity from the interaction of the box with the infinite number of the periodic images. Simulation systems are most commonly neutralized by adding sodium or chloride ions.

## Periodic boundary conditions
Periodic boundary conditions (PBC) are used to approximate a large system by using a small part called a unit cell. The boundary to contain molecules in simulation is needed to preserve thermodynamic properties like temperature, pressure and density.

To implement PBC the unit cell is surrounded by translated copies in all directions to approximate an infinitely large system. When one molecule diffuses across the boundary of the simulation box it reappears on the opposite side. So each molecule always interacts with its neighbours even though they may be on opposite sides of the simulation box. This approach replaces the surface artifacts caused by the interaction of the isolated system with a vacuum with the PBC artifacts which are in general much less severe.

In simulations with PBC the non-bonded interaction cut-off radius should be smaller than half the shortest periodic box vector to prevent interaction of an atom with its image.


### Specifying periodic box ###

#### GROMACS ####
The box specification is integrated into structure file. The [editconf](http://manual.gromacs.org/archive/5.0/programs/gmx-editconf.html) utility is used to set the box parameters:

**-bt**
>box type (triclinic, cubic, dodecahedron, octahedron)<br>

**-box**
> box vectors lengths (a,b,c)<br>
>
**-angles**
> box vectors angles   (bc,ac,ab)<br>


#### NAMD ####
Periodic box is defined by three unit cell vectors:<br>

**cellBasisVector1**<br>
>Default value: 0 0 0

**cellBasisVector2**<br>
>Default value: 0 0 0

**cellBasisVector3**<br>
>Default value: 0 0 0

**extendedSystem**<br>
>NAMD generates a .xsc (eXtended System Configuration) file which contains the periodic cell parameters. If this keyword is used periodic box parameters will be read from .xsc file ignoring cellBasisVectors.
>Value: filename


## Integrating the Equations of Motion.

The integration algorithm advances simulation system by a small step <img src="https://latex.codecogs.com/gif.latex?\delta{t}"/> during which the forces are considered constant. If the time step is small enough the trajectory will be reasonably accurate.

A good integration algorithm for MD should be time-reversible and energy conserving.
### Integration Algorithms
#### The Euler Algorithm

The Euler algorithm uses the second order Taylor expansion to estimate position and velocity at the next time step:

<img src="https://latex.codecogs.com/gif.latex?\vec{r}(t&plus;\delta{t})=\vec{r}(t)&plus;\vec{v}(t)\delta{t}&plus;\frac{1}{2}a(t)\delta{t}^2"/><br>

<img src="https://latex.codecogs.com/gif.latex?\vec{v}(t&plus;\delta{t})=\vec{v}(t)&plus;\frac{1}{2}a(t)\delta{t}"/>

The Euler algorithm is neither time-reversible nor energy conserving and hence rather unfavourable. Nevertheless, the Euler scheme can be used to integrate other equations of motion. For example, GROMACS offers a Euler integrator for Brownian or position Langevin dynamics.

#### The Verlet Algorithm

Using the current positions and forces and the previous positions calculate the positions at the next time step:

<img src="https://latex.codecogs.com/gif.latex?\vec{r}(t&plus;\delta{t})=2\vec{r}(t)-\vec{r}(t-\delta{t})&plus;a(t)\delta{t}^2"/><br>

The Verlet algorithm requires positions at two time steps. It is inconvenient when starting a simulation. While velocities are not needed to compute trajectories, they are useful for calculating observables e.g. the kinetic energy. The velocities can only be computed once the next positions are calculated:

<img src="https://latex.codecogs.com/gif.latex?\vec{v}(t&plus;\delta{t})=\frac{r{(t&plus;\delta{t})-&space;r(t-\delta{t})&space;}}{2\delta{t}}"  />

The Verlet algorithm is time-reversible and energy conserving.

#### The Velocity Verlet Algorithm

The velocities, positions and forces are calculated at the same time according to:

<img src="https://latex.codecogs.com/gif.latex?\vec{r}(t&plus;\delta{t})=\vec{r}(t)&plus;\vec{v}(t)\delta{t}&plus;\frac{1}{2}a(t)\delta{t}^2"/>

<img src="https://latex.codecogs.com/gif.latex?\vec{v}(t&plus;\delta{t})=\vec{v}(t)&plus;\frac{1}{2}[a(t)&plus;a(t&plus;\delta{t})]\delta{t}"/>

The Velocity Verlet algorithm is mathematically equivalent to the original Verlet algorithm. It explicitly incorporates velocity, solving the problem of the first time step in the basic Verlet algorithm. Due to its simplicity and stability is has become the most widely used algorithm in the MD simulations.


#### The Leap Frog Algorithm

Using accelerations of the current time step, compute the velocities at half-time step:

<img src="https://latex.codecogs.com/gif.latex?\vec{v}(t&plus;\frac{1}{2}\delta&space;t)=\vec{v}(t-\frac{1}{2}\delta&space;t)\cdot&space;\delta&space;t&plus;\vec{a}(t)\cdot\delta{t}"  />

Then determine positions at the next time step:

<img src="https://latex.codecogs.com/gif.latex?\vec{r}(t&plus;\delta&space;t)=\vec{r}(t)&plus;\vec{v}(t&plus;\frac{1}{2}\delta&space;t))\cdot&space;\delta&space;t"/>

The Leap Frog algorithm is essentially the same as the Velocity Verlet. The Leap Frog and the Velocity Verlet integrators give equivalent trajectories. The only difference is that the velocities are not calculated at the same time as positions. Leapfrog integration is equivalent to updating positions and velocities at interleaved time points, staggered in such a way that they "leapfrog" over each other.


### Choosing Time Step
Mathematically Vertet family integrators are stable for time steps

<img src="https://latex.codecogs.com/gif.latex?\delta{t}\leq\frac{2}{w}"/>

 where <img src="https://latex.codecogs.com/gif.latex?\omega"/> is angular frequency.<br>
In molecular dynamics stretching of the bonds with the lightest atom H is usually the fastest motion. The period of oscillation of a C-H bond is ~10 fs. Hence Verlet integration will be stable for time steps < 3.2 fs. In practice, the time step of 1 fs is recommended to describe this motion reliably. If the dynamics of hydrogen atoms is not essential for a simulation, bonds with hydrogens can be constrained, and time step increased to 2 fs.

### Specifying Time Step Parameters
parameter             | GROMACS    |  NAMD
----------------------| -----------|----------
time step             | **dt**, ps   |  **timestep**, fs
number of steps       | **nstep**    |  **numsteps**
time of the first step| **tinit**    |  **firsttimestep**

### Specifying Integration Method
#### GROMACS
GROMACS offers several types of integration algorithms that can be selected using the **integrator** keyword.

**md**
> a leap frog algorithm
>
**md-vv**
> a velocity Verlet algorithm
>
**md-vv-avek**
> a velocity Verlet algorithm same as **md-vv** except the kinetic energy is calculated as the average of the two half step kinetic energies. It is more accurate than the md-vv
>
**sd**
> an accurate leap frog stochastic dynamics integrator.
>
**bd**
> a Euler integrator for Brownian or position Langevin dynamics.

#### NAMD

The only available integration method is Verlet. To further reduce the cost of computing short-range nonbonded interactions and full electrostatics, NAMD uses a multiple time-stepping integration scheme controlled by the following keywords:

**nonbondedFreq**
> number of timesteps between nonbonded evaluation<br>

**fullElectFrequency**
>number of timesteps between full electrostatic evaluations<br>




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
