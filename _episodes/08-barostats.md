---
title: "Controlling Pressure"
teaching: 10
exercises: 0
questions:
- "Why do we want to control the pressure of our MD-simulations?"
- "What pressure control algorithms are commonly used?"
- "What are the strengths and weaknesses of these common barostats?"
objectives:
- ""
keypoints:
- ""
---
## Introduction 
The role of pressure control algorithms is to keep pressure in the simulation system constant or to apply an external stress to the simulated system. 

- Pressure is kept on its target value by adjusting the volume of a periodic simulation system.
- Pressure is a force exerted by collision of particles with the walls of a closed container. 
- The virial equation is used to obtain the pressure:

$ \qquad {P}=\frac{NK_{B}T}{V}+\frac{1}{3V}\langle\sum\{r_{ij}F_{ij}}\rangle$

- The first term in this equation describes pressure of an ideal gas (no interaction between molecules). 
- The second contribution comes from internal forces acting on each atom. 
- Well suited for MD because forces are evaluated at each simulation step.

## Pressure Control Algorithms
- Regulate pressure by adjusting the volume 
- In practice barostats do that by scaling coordinates of each atom by a small factor. 
- The methods of maintaining pressure fall into four categories:

1. Weak coupling methods
2. Extended system methods
3. Stochastic methods
4. Monte-Carlo methods

### 1. Weak coupling methods
#### Berendsen pressure bath coupling. 
- Conceptually similar to Berendsen thermostat. 
- Available in all simulation packages. 
- Change the volume by an increment proportional to the difference between the internal pressure and pressure in a weakly coupled bath. 
- Very efficient in equilibrating the system. 

##### Downsides:
- Does not sample the exact NPT statistical ensemble.
- Induces artifacts into simulations of inhomogeneous systems such as aqueous biopolymers or liquid/liquid interfaces.
- Should be avoided for production MD simulations.

The time constant for pressure bath coupling is the main parameter of the Berendsen thermostat. The pressure of the system is corrected such that the deviation exponentially decays with a lifetime defined by this constant. 

Reference: [Molecular dynamics with coupling to an external bath](https://aip.scitation.org/doi/10.1063/1.448118)

### 2. Extended system methods
- Extended system methods originate from the classical theoretical work of [Andersen](https://aip.scitation.org/doi/abs/10.1063/1.439486).
- He included an additional degree of freedom, the volume of a simulation cell.
- Volume adjusts itself to equalize the internal and external pressure. 
- Volume serves as a piston, and is given a fictitious "mass" controlling the decay time of pressure fluctuations.
- Extended system methods are time-reversible. They can be used to integrate backwards, for example, for transition path sampling.
{: .instructor_notes :}

#### Parrinello-Rahman barostat
- Extension of the Andersen method allowing changes in the shape of the simulation cell [[Parrinello and Rahman, 1980]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.45.1196).
- Further extended to include external stresses [[Parrinello and Rahman, 1981]](https://aip.scitation.org/doi/10.1063/1.328693).
- Useful to study structural transformations in solids under external stress.
- Equations of motion are similar to Nosé-Hoover barostat, and in most cases it is used with the Nosé-Hoover thermostat.

##### Downsides:
- Volume may oscillate with the frequency proportional to the piston mass. 

#### Nosé-Hoover barostat
- First application of the method analogous to Andersen's barostat for molecular simulation [[Nosé and Klein, 1983]](https://www.tandfonline.com/doi/abs/10.1080/00268978300102851). 
- The Nosé-Hoover equations of motion are only correct in the limit of large systems.

References: [[Hoover, 1986]](https://pubmed.ncbi.nlm.nih.gov/9897546/),
[[Martyna, 1994]](https://aip.scitation.org/doi/abs/10.1063/1.467468). 

#### MTTK (Martyna-Tuckerman-Tobias-Klein) barostat.
- Extension of the Nosé-Hoover and Nosé-Hoover chain thermostat, performs better for small systems [[Martyna et al., 1996]](https://www.tandfonline.com/doi/abs/10.1080/00268979600100761).

### 3. Stochastic methods
#### Langevin piston pressure control.
- Based on Langevin thermostat. 
- The equations of motion resemble MTTK equations.
- An additional damping (friction) force and stochastic force are introduced. 
- Random collisions eliminate oscillation of the volume associated with the piston mass.

Reference: [Constant pressure molecular dynamics simulation: The Langevin piston method](https://aip.scitation.org/doi/abs/10.1063/1.470648)

|:-:|:-:|
|MTTK and Langevin barostats produce identical ensembles | Langevin barostat oscillates less then MTTK and converges faster due to stochastic collisions and damping.|

![Comparison of Barostats]({{ page.root }}/fig/barostats_comp.png)

Reprinted with permission from [Rogge et al. 2015]({{ page.root }}/reference.html#Rogge-2015), *A Comparison of Barostats for the Mechanical Characterization of Metal−Organic Frameworks*, J Chem Theory Comput. 2015;11: 5583-97. [doi:10.1021/acs.jctc.5b00748](https://doi.org/10.1021/acs.jctc.5b00748). Copyright 2015 American Chemical Society.
{% comment %}
See "fig/barostats_comp -  Copyright Clearance Center.pdf"
{% endcomment %}

#### Stochastic Cell Rescaling
- Improved version of the Berendsen barostat.
- Adds stochastic term to rescaling matrix.
- Produces correct fluctuations of local pressure for NPT ensemble.
- Pressure converges fast without oscillations.
- Can be used for all stages of MD, including production.

Reference: [[Bernetti and Bussi (2020)]][Bernetti-2020]

### 4. Monte-Carlo pressure control. 
- Recently several efficient Monte Carlo methods have been introduced. 
- Sample volume fluctuations at a predefined number of steps at a given constant external pressure. 
- Generate a random volume change, evaluate the potential energy. The volume move is then accepted with the standard Monte-Carlo probability.  
- Do not compute virial, so pressure is not available at the runtime, and not printed in energy files. 


References: 
1. [Molecular dynamics simulations of water and biomolecules with a Monte Carlo constant pressure algorithm](https://www.sciencedirect.com/science/article/abs/pii/S0009261403021687)
2. [Constant pressure hybrid Molecular Dynamics–Monte Carlo simulations](https://aip.scitation.org/doi/10.1063/1.1420460)


## Pitfalls
To ensure stability of a simulation volume must be adjusted very slowly with a small increments at each simulations step. Rapid change of the system size may lead to simulation crash. This can occur, for example when pressure coupling is turned on when you begin simulation from a cold start and turn pressure coupling too early in the heating process. In this case, the difference between the target and the real pressure will be large, the program will try to adjust the density too quickly, and bad things (such as SHAKE failures) are likely to happen.
{: .self_study_text :}
- If the difference between the target and the real pressure is large, the program will try to adjust the density too quickly.
- Rapid change of the system size may lead to simulation crash.
- To ensure stability of a simulation volume must be adjusted very slowly with a small likely 
{: .instructor_notes :}


## Conclusion
Each barostat or thermostat technique has its own limitations and it is your responsibility to choose the most appropriate method or their combination for the problem of interest.


### Selecting barostats in molecular dynamics packages

| Thermostat\MD package | GROMACS                      |  NAMD                    | AMBER         |
|-----------------------|------------------------------|--------------------------|---------------|
| Berendsen             | pcoupl = Berendsen           |  BerendsenPressure on    | barostat = 1  |
| Stoch. cell rescaling | pcoupl = C-rescale           |                          |               |
| Langevin              |                              |  LangevinPiston on       |               |
| Monte-Carlo           |                              |                          | barostat = 2  |
| Parrinello-Rahman     | pcoupl = Parrinello-Rahman   |                          |               |
| MTTK                  | pcoupl = MTTK                |                          |               |


{% comment %}
### References
Below here we resolve reference-style links so that 
[Refname-YYYY] points to the anchor #Refname-YYYY on the {{ page.root }}/reference.html page.

Example:
[Refname-YYYY]: {{ page.root }}/reference.html#Refname-YYYY
{% endcomment %}

[Bernetti-2020]:   {{ page.root }}/reference.html#Bernetti-2020
[Martyna-1996]:    {{ page.root }}/reference.html#Martyna-1996
[Parrinello-1981]: {{ page.root }}/reference.html#Parrinello-1981
