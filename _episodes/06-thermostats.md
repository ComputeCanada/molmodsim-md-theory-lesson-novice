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

![Plot of Maxwell Boltzmann distributions at different temperatures]({{ page.root }}/fig/Maxwell_Boltzmann_distributions.svg)


## Thermodynamic ensembles


## Overview of common Temperature Control Algorithms

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
