<!-- MDTOC maxdepth:6 firsth1:1 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Basic principles of molecular dynamics](#basic-principles-of-molecular-dynamics)   
   - [Force Fields](#force-fields)   
      - [Non-bonded interactions](#non-bonded-interactions)   
         - [The Lennard-Jones potential](#the-lennard-jones-potential)   
         - [The electrostatic potential](#the-electrostatic-potential)   
      - [Bonded interactions](#bonded-interactions)   
         - [The bond potential](#the-bond-potential)   
         - [The angle potential](#the-angle-potential)   
         - [The torsion angle potential](#the-torsion-angle-potential)   
      - [Boundary conditions](#boundary-conditions)   
      - [Truncation of interactions](#truncation-of-interactions)   
      - [Balancing of Charges](#balancing-of-charges)   
      - [Integrating the equations of motion](#integrating-the-equations-of-motion)   

<!-- /MDTOC -->




# Basic principles of molecular dynamics
Molecular dynamics (MD) simulations are widely used as tools to investigate structure and dynamics of proteins, nucleic acids, carbohydrates, lipids, nanoparticles and liquid/solid interfaces under a wide variety of conditions. MD is the simulation of atomic positions in time accomplished by solving classical Newton's equation of motion stating that the rate of change of momentum <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{p}" title="\vec{p}" /></a> of an object equals the force <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{F}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{F}" title="\vec{F}" /></a> acting on it:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=\frac{d\vec{p}}{dt}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=\frac{d\vec{p}}{dt}" title="\vec{F}=\frac{d\vec{p}}{dt}" /></a>

 As the force on an object is the negative of the derivative of the potential energy function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=-\nabla&space;U(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=-\nabla&space;U(\vec{r})" title="\vec{F}=-\nabla U(\vec{r})" /></a>

  the potential energy function acts as a cornerstone of the MD simulations.

## Force Fields
A force field is a set of energy functions and parameters used to calculate the potential energy *U* of a system of atoms and/or molecules in molecular dynamics simulations. The potential energy includes intermolecular (non-bonded) and intramolecular (bonded) interactions.

<a href="https://www.codecogs.com/eqnedit.php?latex=U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" title="U(\vec{r})=\sum U_{bonded}(\vec{r})+\sum U_{nonbonded}(\vec{r})" /></a>

### Non-bonded interactions
#### The Lennard-Jones potential
#### The electrostatic potential
### Bonded interactions
#### The bond potential
#### The angle potential
#### The torsion angle potential
### Boundary conditions
### Truncation of interactions
### Balancing of Charges
### Integrating the equations of motion
