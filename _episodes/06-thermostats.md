---
title: "Controlling Temperature"
teaching: 25
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

Temperature is a macroscopic property and thermodynamics tells us that on the molecular level, 
the temperature of a system is basically defined by the average kinetic energy of all the 
particles (atoms, molecules) that make up the system.
Statistical thermodynamics tells us that a system that contains a certain amount of energy 
and is in equilibrium, will have all of its energy distributed in the most probable way. 
This distribution is called the Maxwell-Boltzmann distribution. 
This means that in a system of particles, these particles won’t all have the same velocity 
but rather the velocities will follow a distribution that depends on their mass and the 
temperature of the system.


![Plot of Maxwell-Boltzmann distributions]({{ page.root }}/fig/Maxwell_Boltzmann_distributions.svg)
The plot above shows the Maxwell-Boltzmann distributions for the velocities of oxygen atoms
from water molecules that have been simulated at three different temperatures (280K, 320K and 360K).

## Thermodynamic ensembles

MD simulations are usually simulate one of the following thermodynamic ensembles:

1. The *microcanonical* or constant **NVE** ensemble, where the number of particles (*N*), the
   system's Volume (*V*) and the energy (*E*) are kept constant during the simulation.
2. The *canonical* or constant **NVT** ensemble, where the number of particles, the Volume and
   the temperature (*T*) are kept constant.
3. The *isothermal-isobaric* or **NPT** ensemble, where number of particles, the system's pressure
   (*P*) and temperature are kept constant.

Simulations using the *NVE* ensemble is relatively easy to achieve, as long as the MD code manages 
to conserve the energy of the system.
The *NVT* ensemble is more practically relevant, as in the real world it is much easier to manage 
the temperature of a system rather than it's energy.  To achieve this in our simulation, we need
to use a temperature control algorithm, which is commonly called a *thermostat*.
As many experiments in the lab are carried out at constant ambient pressure, rather than in a
fixed/confined volume, the *NPT* ensemble is also widely used for simulations.  
In addition to using a thermostat to control the temperature, we also need to use a pressure
control algorithm, often called a *barostat*.

The *grand canonical ensemble* however, where the chemical potential *&mu;* is kept constant,
requires that the number of particles is allowed to change, which is not supported by most 
MD packages.


## Overview of common Temperature Control Algorithms

Temperature control algorithms or *thermostats* have been described as the "necessary evil"
[[Wong-ekkabut-2016]({{ page.root }}/reference.html#Wong-ekkabut-2016)].

Their role is to allow energy to enter and leave the simulated system to keep its temperature
constant.  Thermostats do that by adjusting the velocities of one or more particles.



### Andersen thermostat
[Andersen-1980]({{ page.root }}/reference.html#Andersen-1980)

### Berendsen thermostat
[Berendsen-1984]({{ page.root }}/reference.html#Berendsen-1984)

### Bussi's stochastic velocity rescaling thermostat
[Bussi-2007]({{ page.root }}/reference.html#Bussi-2007)

### Nosé-Hoover thermostat
[Nose-1984]({{ page.root }}/reference.html#Nose-1984), [Hoover-1985]({{ page.root }}/reference.html#Hoover-1985)

### Nosé-Hoover-chains
[Martyna-1992]({{ page.root }}/reference.html#Martyna-1992)

## Important parameters


## Conclusions
