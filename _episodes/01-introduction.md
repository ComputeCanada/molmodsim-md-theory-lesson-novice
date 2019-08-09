<!-- MDTOC maxdepth:6 firsth1:1 numbering:1 flatten:0 bullets:1 updateOnSave:1 -->

- 1. [Basic principles of molecular dynamics](#basic-principles-of-molecular-dynamics)   
   - 1.1. [Force Fields](#force-fields)   
      - 1.1.1. [Non-bonded interactions](#non-bonded-interactions)   
         - 1.1.1.1. [The Lennard-Jones potential](#the-lennard-jones-potential)   
         - 1.1.1.2. [The electrostatic potential](#the-electrostatic-potential)   
      - 1.1.2. [Bonded interactions](#bonded-interactions)   
         - 1.1.2.1. [The bond potential](#the-bond-potential)   
         - 1.1.2.2. [The angle potential](#the-angle-potential)   
         - 1.1.2.3. [The torsion angle potential](#the-torsion-angle-potential)   
   - 1.2. [Boundary conditions](#boundary-conditions)   
   - 1.3. [Truncation of interactions](#truncation-of-interactions)   
   - 1.4. [Balancing of charges](#balancing-of-charges)   
   - 1.5. [Integrating the equations of motion](#integrating-the-equations-of-motion)   
   - 1.6. [MD software available on CC clusters](#md-software-available-on-cc-clusters)   
      - 1.6.1. [AMBER](#amber)   
      - 1.6.2. [GROMACS](#gromacs)   
         - 1.6.2.1. [Force fields implemented in GROMACS:](#force-fields-implemented-in-gromacs)   
      - 1.6.3. [NAMD](#namd)   
         - 1.6.3.1. [Force fields implemented in NAMD:](#force-fields-implemented-in-namd)   
      - 1.6.4. [LAMMPS](#lammps)   
      - 1.6.5. [DL_POLY](#dl_poly)   

<!-- /MDTOC -->

# Basic principles of molecular dynamics
Molecular dynamics (MD) simulations are widely used as tools to investigate structure and dynamics of proteins, nucleic acids, carbohydrates, lipids, nanoparticles and liquid/solid interfaces under a wide variety of conditions. MD is the simulation of atomic positions in time accomplished by solving classical Newton's equation of motion stating that the rate of change of momentum <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{p}" title="\vec{p}" /></a> of an object equals the force <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\vec{F}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\vec{F}" title="\vec{F}" /></a> acting on it:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=\frac{d\vec{p}}{dt}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=\frac{d\vec{p}}{dt}" title="\vec{F}=\frac{d\vec{p}}{dt}" /></a>

The potential energy function *U* acts as a cornerstone of the MD simulations because it allows to calculate the forces. The force on an object is the negative of the derivative of the potential energy function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}=-\nabla&space;U(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}=-\nabla&space;U(\vec{r})" title="\vec{F}=-\nabla U(\vec{r})" /></a>


## Force Fields
A force field is a set of empirical energy functions and parameters used to calculate the potential energy *U* of a system of atoms and/or molecules in molecular dynamics simulations. The potential energy includes non-bonded and bonded interactions.

<a href="https://www.codecogs.com/eqnedit.php?latex=U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U(\vec{r})=\sum&space;U_{bonded}(\vec{r})&plus;\sum&space;U_{nonbonded}(\vec{r})" title="U(\vec{r})=\sum U_{bonded}(\vec{r})+\sum U_{nonbonded}(\vec{r})" /></a>


### Non-bonded interactions
#### The Lennard-Jones potential
The Lennard-Jones potential approximates the potential energy of non-elecrostatic interaction between a pair of non-bonding atoms or molecules with a simple mathematical function:

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{LJ}(r)=\frac{C_{12}}{r^{12}}&plus;\frac{C_{6}}{r^{6}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V_{LJ}(r)=\frac{C_{12}}{r^{12}}&plus;\frac{C_{6}}{r^{6}}" /></a>
It is commonly expressed in terms of the well depth <img src="https://latex.codecogs.com/gif.latex?\epsilon" /> (the measure of the strength of the interaction) and the *van der Waals* radius <img src="https://latex.codecogs.com/gif.latex?\sigma" /> (the distance at which the intermolecular potential between the two particles is zero).
<img src="https://latex.codecogs.com/gif.latex?V_{LJ}(r)=4\epsilon\left&space;[&space;\left&space;(&space;\frac{\sigma}{r}\right&space;)^{12}-&space;\left&space;(&space;\frac{\sigma}{r}\right&space;)^{6}&space;\right&space;]" title="V_{LJ}(r)=4\epsilon\left [ \left ( \frac{\sigma}{r}\right )^{12}- \left ( \frac{\sigma}{r}\right )^{6} \right ]" />

The *LJ* coefficients *C* are related to the <img src="https://latex.codecogs.com/gif.latex?\sigma"/>  and the <img src="https://latex.codecogs.com/gif.latex?\epsilon"/>  with the equations:

 <img src="https://latex.codecogs.com/gif.latex?\inline&space;C_{12}=4\epsilon\sigma^{12}, C_{6}=4\epsilon\sigma^{6}"/>

To describe all *LJ* interactions in a simulations system the matrix of the pairwise interactions is constucted. The *LJ* interaction between different types of atoms are computed by combining the *LJ* parameters.

Combining rules:

1. geometric mean (OPLS)
2. Lorentz–Berthelot (CHARM, AMBER)
3. hybrid: Lorentz–Berthelot for H and the Waldman–Hagler for other elements  (AMBER-ii, FF for perfluoroalkanes)

#### The electrostatic potential
### Bonded interactions
#### The bond potential
#### The angle potential
#### The torsion angle potential
## Boundary conditions
## Truncation of interactions
The most computationally demanding part of a molecular dynamics simulation is the calculation of the nonbonded terms of the potential energy function. As non-bonded energy terms between every pair of atoms should be evaluated, the number of calculations increases as the square of the number of atoms. To speed up the computation, only the interactions between two atoms separated by a distance less than a pre-defined cutoff distance are evaluated. There are several different ways to truncate the non-bonded interaction.
## Balancing of charges
## Integrating the equations of motion

## MD software available on CC clusters
### AMBER
### GROMACS
#### Force fields implemented in GROMACS:
- AMBER: 94, 96, 99, 99sb, 99sb-ildn, 03, GS (amberGS is amber94 with both backbone torsion potentals set to 0).
- CHARMM: 27 (optimized for proteins and nucleic acids).
- [CHARMM for GROMACS](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (CHARMM36, CGenFF)
- GROMOS: 43a1, 43a2, 45a3, 53a5, 53a6, 54a7.
- OPLS-AA (OPLS-AA implemented in GROMACS is actually OPLS-AA/L. OPLS-AA/L uses OPLS-AA atom types with the torsions and impropers refitted to QM calculations at the HF/6-31G** level followed by single-point LMP2/cc-pVTZ(-f))
### NAMD
#### Force fields implemented in NAMD:
- [AMBER](http://ambermd.org/AmberModels.php) (amber format tolopogy prepared with AMBERTOOLS)
- [CHARMM](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm) (charmm format, topology prepared with psfgen)
- [CHARMM-Drude](http://mackerell.umaryland.edu/charmm_drude_ff.shtml) (the polarizable force field based on the classical Drude oscillator model, charmm format)
- GROMOS (from GROMACS distribution, topology prepared with pdb2gmx)
- [OPLS-AA/M](http://zarbi.chem.yale.edu/oplsaam.html) (charmm format)

1. The organic solvents with OPLS force field generate slightly better properties than those with GAFF. [ C. Caleman, P.J. van Maaren, M. Hong, J. S. Hub, L. T. Costa and D. van der Spoel. Force Field Benchmark of Organic Liquids: Density, Enthalpy of Vaporization, Heat Capacities, Surface Tension, Isothermal Compressibility, Volumetric Expansion Coefficient, and Dielectric Constant, J. Chem. Theor. Comput., 8, 61-74 (2012) ].
### LAMMPS
### DL_POLY
