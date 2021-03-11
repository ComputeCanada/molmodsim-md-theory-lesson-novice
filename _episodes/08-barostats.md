---
title: "Controlling Pressure"
teaching: 25
exercises: 0
questions:
- "Why do we want to control the temperature of our MD-simulations?"
- "What Pressure control algorithms are commonly used?"
- "What are the strengths and weaknesses of these common barostats?"
objectives:
- ""
keypoints:
- ""
---

The role of pressure control algorithms is to keep pressure in the simulation system constant or to apply an external stress to the simulated system. 

### How is pressure kept constant?
Pressure is kept on its target value by adjusting the volume of a periodic simulation system. In practice this is done by multiplying coordinates of each atom in the system by a small scaling factor proportional to the difference between the target and the real pressure. 

What defines pressure is a simulation system? Pressure is a force exerted by collision of particles with the walls of a closed container. The virial equation is commonly used to obtain the pressure from a molecular dynamics simulation. According to this equetion pressure in a simulation has two components:

$P=\frac{NK_{B}T}{V}+\frac{1}{3V}\langle\sum\{r_{ij}F_{ij}}\rangle$

The first term in this equation describes pressure of an ideal gas (meaning no interaction between molecules). The secong contribution comes from internal forces acting on each atom. The virial equation is well suited for MD since forces are evaluated at each simulation step anyway, and they are readily available.

As with temperature control, there are different algorithms of pressure control for MD simulation. 


# Overview of common Pressure Control Algorithms

 Barostats regulate pressure by adjusting the volume of the simulatted system. In pactice this done by scaling coordinates of each atom is a system by a small factor, so that the size of the system changes. The methods of maintaining pressure fall into categories similar to categories of temperature regulation.

## Weak coupling methods
### Berendsen pressure bath coupling. 
Berendsen barostat is conceptually similar to Berendsen thermostat. This thermostat is availablein all simulation packages. To regulate pressure this algorithm changes the volume by an increment proportional to the difference between the internal pressure and pressure in a weakly coupled bath. Berendsen thermostat is very efficient in equilibrating the system. However, it does not sample the exact NPT statistical ensemble. Often pressure fluctuations with this baristat are smaller then they should be. It is also known that this algorithm induces artifacts into simulations of inhomogeneous systems such as aqueous biopolymers or liquid/liquid interfaces.

The time constant for pressure bath coupling is the main parameter of the Berendsen thermostat. The pressure of the system is corrected such that the deviation exponentially decays with a lifetime defined by this constant. 

Reference: [Molecular dynamics with coupling to an external bath](https://aip.scitation.org/doi/10.1063/1.448118)

## Extended system methods

In the classical work Andersen proposed a method with an extended variable. It is based on including an additional degree of freedom corresponding to the volume of a simulation cell which adjusts itself to equalize the internal and external pressure was originally developed by [Andersen](https://aip.scitation.org/doi/abs/10.1063/1.439486). In this method an additional degree of freedom serves as a piston, and is given a fictious "mass". The choice of piston "mass" determines the decay time of the volume fluctuations. The Parrinello-Rahman barostat, the Nose-Hoover barostat, and the Martyna-Tuckerman-Tobias-Klein (MTTK) are all based on the Andersen barostat.

### Parinello-Rahman barostat
 The Andersen method was extended by [Parrinello and Rahman](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.45.1196) to allow changes in the shape of the simulation cell, and later it was further extended by the same authors to include external stresses [Parrinello and Rahman](https://aip.scitation.org/doi/10.1063/1.328693). This metod allows additionally for a dynamic shape change and is useful to study structural transformations in solids under external stress. Equations of motion are similar to Nose-Hoover barostat, and in most cases you would combine the Parrinello-Rahman barostat with the Nosé-Hoover thermostat.

The drawback is that decay of the volume fluctuations may oscillate with the frequency proportional to the piston mass. 

### Nose-Hoover barostat
[Nose and Klein](https://www.tandfonline.com/doi/abs/10.1080/00268978300102851) were the first to apply method analogous to Andersen's barostat for molecular simulation.  This variation was further improved by [Hoover](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.34.2499). It is widely used now under the name  "Nosé-Hoover barostat" [Constant pressure molecular dynamics algorithms](https://aip.scitation.org/doi/abs/10.1063/1.467468). 

### MTTK [Martyna-Tuckerman-Tobias-Klein](https://www.tandfonline.com/doi/abs/10.1080/00268979600100761) barostat.
However, Nose-Hoover method was not free of issues. It was observed that the Nose-Hoover equations of motion are only correct in the limit of large systems. To solve this problem [Martyna-Tuckerman-Tobias-Klein](https://www.tandfonline.com/doi/abs/10.1080/00268979600100761) introduced their own equations of motion. The MTTK method is a natural extension of the Nose−́Hoover and Nose−́Hoover chain thermostat,  and in most cases you would combine the MTTK barostat with the Nosé-Hoover thermostat. The nice feature of this thermostat is that it is time-reversible. This means that it can be used to integrate backwards  e.g., as needed for transition path sampling.
   
### Langevin piston pressure control.
Langevin piston barostat is based on Langevin thermostat. The equations of motion resemble MTTK equations, but an additional damping force and stochastic force are introduced. A suitable choice of collision frequency then eliminates the unphysical oscillation of the volume associated with the piston mass. In this way it is similar to the weak coupling Berendsen algorithm, but in contrast it yields the correct ensemble. 

Reference: [Constant pressure molecular dynamics simulation: The Langevin piston method](https://aip.scitation.org/doi/abs/10.1063/1.470648)


MTTK and Langevin produce identical ensembles. Langevin barostat oscillates less then MTTK and converges faster due to stochastic collisions and damping.

[A Comparison of Barostats for the Mechanical Characterization of Metal−Organic Frameworks](https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b00748)


### Monte-Carlo pressure control. 
Recently several efficient modified hybrid Monte Carlo methods have been introduced.

Monte Carlo pressure control samples volume fluctuations every at a predefined number of steps at a given constant external pressure. It involves generation of a random volume change from a uniform distribution followed by evaluation of the potential energy of the trial system. The volume move is then accepted with the standard Monte-Carlo probability.  

Reference: [Molecular dynamics simulations of water and biomolecules with a Monte Carlo constant pressure algorithm](https://www.sciencedirect.com/science/article/abs/pii/S0009261403021687)

[Constant pressure hybrid Molecular Dynamics–Monte Carlo simulations](https://aip.scitation.org/doi/10.1063/1.1420460)

### Selecting barostats in molecular dynamics packages

| Thermostat/MD package | GROMACS                      |  NAMD                    | AMBER         |
|-----------------------|------------------------------|--------------------------|---------------|
| Berendsen             | pcoupl = Berendsen           |  BerendsenPressure on    | barostat = 1  |
| Langevin              |                              |  LangevinPiston on       |               |   
| Monte-Carlo           |                              |                          | barostat = 2  |   
| Parrinello-Rahman     | pcoupl = Parrinello-Rahman   |                          |               |   
| MTTK                  | pcoupl = MTTK                |                          |               | 



### Pitfalls

To ensure stabitily of a simulation volume must be adjusted very slowly with a small increments at each simulations step. Rapid change of the system size may lead to simulation crash. This can occur, for example when pressure coupling is turned on when you begin simulation from a cold start and turn pressure coupling too early in the heating process. In this case, the difference between the target and the real pressure will be large, the program will try to adjust the density too quickly, and bad things (such as SHAKE failures) are likely to happen.

#### Conclusion
Each barostat or thermostat technique has its own limitations and it is your responsibility to choose the most appropriate method or their combination for the problem of interest.
