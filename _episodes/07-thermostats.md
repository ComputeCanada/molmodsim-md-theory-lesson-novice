---
title: "Controlling Temperature"
teaching: 30
exercises: 0
questions:
- "What is temperature on the molecular level?"
- "Why do we want to control the temperature of our MD-simulations?"
- "What temperature control algorithms are commonly used?"
- "What are the strengths and weaknesses of these common thermostats?"
objectives:
- "Remind us how Maxwell-Boltzmann distributions relate to temperature."
- "Quickly review thermodynamic ensembles."
- "Learn about different thermostats and get an idea how they work."
- "Learn where thermostats that don't produce correct thermodynamic ensembles can still be very useful."
keypoints:
- "On the molecular level, temperature manifests itself as a number of particles having a certain average kinetic energy."
- "Thermodynamic ensembles are ... FIXME"
- "Some temperature control algorithms (e.g. the Berendsen thermostat) fail to produce kinetic energy distributions that represent a correct thermodynamic ensemble."
- "Other thermostats, like Nosé-Hoover, produce correct thermodynamic ensembles but can take long to converge."
- "Even though the the Berendsen thermostat fails to produce correct thermodynamic ensembles, it can be useful for system relaxation as it is robust and converges fast."
---

## Temperature at the molecular level (Maxwell-Boltzmann distributions)

Temperature is a macroscopic property and thermodynamics tells us that on the molecular level, the temperature of a system is  defined by the average kinetic energy of all the particles (atoms, molecules) that make up the system. Statistical thermodynamics tells us that a system that contains a certain amount of energy and is in equilibrium, will have all of its energy distributed in the most probable way. This distribution is called the Maxwell-Boltzmann distribution. This means that in a system of particles, these particles won’t all have the same velocity but rather the velocities will follow a distribution that depends on their mass and the temperature of the system.

The speed distribution of ideal gas obeys the relationship:

$ f_v(v) = \left(\frac{m}{2\pi k_B T}\right)^{3/2} \; 4\pi v^2 \; \exp({-mv^2/2k_B T}) $

![Plot of velocity distributions]({{ page.root }}/fig/MB_ideal_gas.svg)

The plot above shows how the Maxwell-Boltzmann distributions of ideal gas depend on temperature (left) and mass (right). The python code is in: code/thermostats/MB_dist_gas.py


![Plot of Maxwell-Boltzmann distributions]({{ page.root }}/fig/Maxwell_Boltzmann_distributions.svg)

The plot above shows the distributions for the velocities of oxygen atoms from water molecules that have been simulated at three different temperatures (280K, 320K and 360K).

## Thermodynamic ensembles

MD simulations are usually simulate one of the following thermodynamic ensembles:

1. The *microcanonical* or constant **NVE** ensemble, where the number of particles (*N*), the
   system's Volume (*V*) and the energy (*E*) are kept constant during the simulation.
2. The *canonical* or constant **NVT** ensemble, where the number of particles, the Volume and
   the temperature (*T*) are kept constant.
3. The *isothermal-isobaric* or **NPT** ensemble, where number of particles, the system's pressure
   (*P*) and temperature are kept constant.

Simulation of *NVE* ensemble is relatively easy to achieve, as long as the MD code manages to conserve the energy of the system. The *NVT* ensemble is more practically relevant, as in the real world it is much easier to manage the temperature of a system rather than it's energy.  To achieve this in our simulation, we need to use a temperature control algorithm, which is commonly called a *thermostat*. As many experiments in the lab are carried out at constant ambient pressure, rather than in a fixed/confined volume, the *NPT* ensemble is also widely used for simulations. In addition to using a thermostat to control the temperature, we also need to use a pressure control algorithm, often called a *barostat*.

The *grand canonical ensemble* however, where the chemical potential *&mu;* is kept constant, requires that the number of particles is allowed to change, which is not supported by most MD packages.


## Overview of common Temperature Control Algorithms

Temperature control algorithms or *thermostats* have been described as the "necessary evil".
[[Wong-ekkabut-2016][Wong-ekkabut-2016]].

Their role is to allow energy to enter and leave the simulated system to keep its temperature constant. In practice thermostats do that by adjusting the velocities of a subset of particles. The methods of maintaining temperature fall into several categories:

- Strong coupling methods
  - scale velocity to give exact desired temperature

- Weak coupling methods
  - scale velocity in direction of desired temperature

- Stochastic methods
  - change a system variable stochastically, so that the result matches the desired velocity distribution function 

- Extended system dynamics
  - extend degrees of freedom to include temperature

## Strong coupling methods
### Velocity rescaling 
The idea is that we rescale the velocities at each step (or after a preset number of steps) so that we get the desired target temperature. 
### Velocity reassignment
Using this method, all of the velocities in the system are periodically reassigned (by assigning new randomized velocities) so that the entire system is set to the desired temperature. 

Both methods do not generate the correct canonical ensemble. If the velocities are rescaled at every time step the kinetic energy will not fluctuate in time. This is inconsistent with statistical mechanics of the canonical ensemble. These methods are not recommended for equilibrium dynamics, but they are useful for system heating or cooling. Both have advantages and downsides. If you want to heat up a system, then rescaling will make hot spots even hotter, and this is not desired. Temperature reassignment avoids this problem, but the kinetic energy of particles is no longer consistent with their potential energy, and thus needs to be redistributed. 

## Weak coupling methods
### Berendsen thermostat
The Berendsen thermostat at every simulation step rescales the velocities of all particles to remove a predefined fraction of the difference from the predefined temperature. [[Berendsen-1984][Berendsen-1984]. Conceptually this thermostat works analogous to coupling the simulated system to a fictitious heat bath kept at some constant temperature. The rate of temperature equilibration with this thermostat is controlled by strength of the coupling. This makes the Berendsen thermostat a predictably converging and robust thermostat, which can be very useful when allowing the system to relax, for example when starting a molecular dynamics simulation after after an energy minimization.

However, a drawback of Berendsen thermostat is that it cannot be mapped onto a specific thermodynamic ensemble. It has been shown that it produces an energy distribution with a lower variance than of a true canonical ensemble because it disproportionally samples kinetic energies closer to $$T_0$$ than would be observed in the true Maxwell-Boltzmann distribution [[Basconi-2013][Basconi-2013] and [Shirts-2013][Shirts-2013]].  Therefore, the Berendsen thermostat should be avoided for production MD simulations in most cases.

Heat flows between the simulation system and the heat bath with the rate defined by a time constant $$\tau_T$$ where small values of $$\tau_T$$ mean tight coupling and fast temperature equilibration between the simulation and the bath. Large values of $$\tau_T$$ mean slow equilibration and weak coupling. The time constant for heat bath coupling for the system is measured in picoseconds, and the default valu is usually 1 ps. 

## Stochastic methods

### Andersen thermostat
The Andersen thermostat controls the temperature by assigning a subset of atoms new velocities that are randomly selected from the Maxwell-Boltzmann distribution for the target temperature. The probability for a given particle to have it's velocity reassigned at each step is small and can be expressed as $$\Delta t / \tau_T$$ where $$\Delta t$$ is the time step.  This means that effectively on average every atom experiences a stochastic collision with a virtual particle every $$\Delta t$$ [[Andersen-1980][Andersen-1980]]. A variant of the Andersen thermostat in which the velocities of all atoms are randomized every $$\Delta t$$ is termed the  "massive Andersen" thermostat [[Basconi-2013][Basconi-2013]].

The Andersen thermostat correctly samples the canonical ensemble, however momentum is not conserved with this thermostat.  Due to velocity randomization it can impair some correlated motions and thus slow down the kinetics of the system. This algorithm therefore is not recommended when studying kinetics or diffusion properties of the system. This applies to all stochastic methods.

The number of steps between randomization of velocities to a distribution is an important parameter. Too high a collision rate (short interval between randomization) will slow down the speed at which the molecules explore configuration space, whereas too low a rate (long interval between randomization)  means that the canonical distribution of energies will be sampled slowly. 

### The Lowe-Andersen thermostat 
A variant of the Andersen thermostat that conserves momentum. This method perturbs the system dynamics to a far less than the original Andersen method. This alleviates suppressed diffusion in the system. [[Koopman-2006][Koopman-2006]].

### Bussi stochastic velocity rescaling thermostat
Extension of the Berendsen method corrected for sampling the canonical distribution. The velocities of all the particles are rescaled by a properly chosen random factor.
[[Bussi-2007][Bussi-2007]]

### Langevin thermostat
Langevin dynamics mimics the viscous aspect of a solvent and interaction with the environment by adding two force terms to the equation of motion: a frictional force and a random force. The frictional force and the random force combine to give the correct canonical ensemble.

The amount of friction is controlled by the damping coefficient.  If its value is high, atoms will experience too much 
unnatural friction, however if the coefficient is too low, your system will fluctuate too far away from the desired temperature. The default value is usually 1/ps.

## Extended system thermostats
### Nosé-Hoover thermostat
The extended system method originally introduced by Nose and subsequently developed by Hoover. The idea is to consider the heat bath as an integral part of the system by addition of an artificial variable associated with a fictional "heat bath mass" to the equations of motion. An important feature of this method is that the temperature can be controlled without involving random numbers. Thus correlated motions are not impaired and this method describes kinetics and diffusion properties better. Because the time-evolution of the added variable is described by a second-order equation, heat may flow in and out of the system in an oscillatory fashion, leading to nearly periodic temperature fluctuations with the frequency proportional to the "heat bath mass". [Nose-1984][Nose-1984], [Hoover-1985][Hoover-1985]. The drawback of this thermostat is that it was shown to impart the canonical distribution as well as ergodicity  (space-filling) of the thermostatted system. 

The time constant parameter in this thermostat controls the period of temperature fluctuations at equilibrium. 

### Nosé-Hoover-chains
A modification of the Nosé-Hoover thermostat which includes not a single thermostat variable but a chain of variables.  [Martyna-1992][Martyna-1992]. Nose-Hoover thermostat with one variable does not guarantee ergodicity, especially for small or stiff systems. Chaining variables behaves better for small or stiff cases, however an infinite chain is required to completely correct these issues. 

### Selecting thermostats in molecular dynamics packages

| Thermostat/MD package | GROMACS                      |  NAMD                    | AMBER         |
|-----------------------|------------------------------|--------------------------|---------------|
| velocity rescaling    |                              |  reascaleFreq (steps)    |               |
| velocity reassignment |                              |  reassignFreq (steps)    |               |
| Andersen              | tcoupl = andersen            |                          |               |
| massive-Andersen      | tcoupl = andersen-massive    |                          | ntt = 2       |
| Lowe-Andersen         |                              |  loweAndersen on         |               |
| Berendsen             | tcoupl = berendsen           |  tCouple on              | ntt = 1       |
| Langevin              |                              |  langevin on             | ntt = 3       |
| Bussi                 | tcoupl = V-rescale           |  stochRescale  on        |               |
| Nose-Hoover           | tcoupl = nose-hoover         |                          |               |
| Nose-Hoover-chains    | nh-chain-length (default 10) |                          |               |


#### Global and local thermostats
Global thermostats control temperature of all atom in a system uniformly. This may lead to cold solute and hot solvent due to a slow heat transfer.

With local thermostats it is possible to control temperature in  selected groups of atoms independently. This works well for large solutes, but if solute is small this approach may result in large fluctuations and hence unphysical dynamics.

With NAMD it is possible to set coupling coefficients for each atom.

langevinFile 
langevinCol

tCoupleFile
tCoupleFCol

With GROMACS temperature of a selected groups of atoms can be controlled independently tc-grps groups to couple separately to temperature bath

## Conclusions

[Andersen-1980]: {{ page.root }}/reference.html#andersen-1980
[Basconi-2013]: {{ page.root }}/reference.html#basconi-2013
[Berendsen-1984]: {{ page.root }}/reference.html#berendsen-1984
[Bussi-2007]: {{ page.root }}/reference.html#bussi-2007
[Hoover-1985]: {{ page.root }}/reference.html#hoover-1985
[Koopman-2006]: {{ page.root }}/reference.html#koopman-2006
[Martyna-1992]: {{ page.root }}/reference.html#martyna-1992
[Nose-1984]: {{ page.root }}/reference.html#nose-1984
[Shirts-2013]: {{ page.root }}/reference.html#shirts-2013
[Wong-ekkabut-2016]: {{ page.root }}/reference.html#wong-ekkabut-2016
