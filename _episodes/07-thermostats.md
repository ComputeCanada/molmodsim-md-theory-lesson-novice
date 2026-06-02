---
title: "Controlling Temperature"
teaching: 20
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
- "Some temperature control algorithms (e.g. the Berendsen thermostat) fail to produce kinetic energy distributions that represent a correct thermodynamic ensemble."
- "Other thermostats, like Nosé-Hoover, produce correct thermodynamic ensembles but can take long to converge."
- "Even though the the Berendsen thermostat fails to produce correct thermodynamic ensembles, it can be useful for system relaxation as it is robust and converges fast."
---

## Temperature at the molecular level.
- Temperature is defined by the average kinetic energy of all the particles.
- A system in equilibrium, will have all of its energy distributed in the most probable way. 
- Particles in a system at equilibrium don't all have the same velocity. 
- Velocities follow a distribution that depends on their mass and the temperature of the system:

$$f_v(v)=\left(\frac{m}{2\pi k_B T}\right)^{3/2}\cdot4\pi v^2\cdot\exp({-mv^2/2k_B T})$$
{: .math-center }

#### <center>Maxwell-Boltzmann distributions</center>

|:--:|:--:|:--:|
|Krypton at different temperatures||Noble gases with different mass at 298K|

![Plot of velocity distributions]({{ page.root }}/fig/MB_ideal_gas.svg)

Velocity distributions obtained from MD simulations of water at different temperature.

![Plot of Maxwell-Boltzmann distributions]({{ page.root }}/fig/Maxwell_Boltzmann_distributions.svg){: width="320"} 


### Thermodynamic Ensembles
MD simulations typically model one of the following thermodynamic ensembles:

1. **Microcanonical ensemble (NVE)** — constant number of particles (N), volume (V), and energy (E)
2. **Canonical ensemble (NVT)** — constant number of particles (N), volume (V), and temperature (T)
3. **Isothermal–isobaric ensemble (NPT)** — constant number of particles (N), pressure (P), and temperature (T)

### Which ensemble should I use?
- NVE simulations are straightforward if total energy is conserved.
- NVT simulations use a **thermostat** to maintain constant temperature.
- NPT simulations use both a **thermostat** and a **barostat** to maintain constant temperature and pressure.

### Why do we need thermostats?
- Allow energy to enter and leave the simulated system to keep its temperature constant. 
- In practice thermostats do that by adjusting the velocities of a subset of particles. 
- The methods of maintaining temperature fall into four categories:
### Categories of thermostats
1. Strong coupling methods
2. Weak coupling methods
3. Stochastic methods
4. Extended system dynamics


| Method  | Thermostat         | Good for equilibration? | Good for production? | Samples correct NVT? |
| ------------------ | ----------------------- | -------------------- | -------------------- |
| Strong  coupling | Velocity rescaling | ✓                       | ✗                    | ✗                    |
| Weak coupling | Berendsen          | ✓                       | ✗                    | ✗                    |
| Stochastic | Andersen           | ✓                       | Sometimes            | ✓                    |
| Stochastic | Langevin           | ✓                       | ✓                    | ✓                    |
| Stochastic | Bussi              | ✓                       | ✓                    | ✓                    |
| Extended   | Nosé-Hoover        | ✓                       | ✓                    | ✓                    

### 1. Strong coupling methods
#### Velocity rescaling 
- Rescale the velocities at each step (or after a preset number of steps) to get the desired target temperature. 

#### Velocity reassignment
- Periodically assign new randomized velocities so that the entire system is set to the desired temperature. 

- Both methods do not generate the correct canonical ensemble. 
- Not recommended for equilibrium dynamics
- Useful for heating or cooling. 

##### Downsides:
- Rescaling will make hot spots even hotter.
- Temperature reassignment avoids this problem, but the kinetic energy of particles is no longer consistent with their potential energy, and thus needs to be redistributed. 

### 2. Weak coupling methods
#### Berendsen thermostat
- Rescale the velocities of all particles to remove a fraction of the difference from the predefined temperature.
- The rate of temperature equilibration is controlled by strength of the coupling. 
- The Berendsen thermostat a predictably converging and robust thermostat.
- Very useful when allowing the system to relax.

##### Downsides:
- Cannot be mapped onto a specific thermodynamic ensemble. 
- Produces an energy distribution with a lower variance than of a true canonical ensemble [[Basconi-2013][Basconi-2013] and [Shirts-2013][Shirts-2013]]. 
- Should be avoided for production MD simulations.

Heat flows between the simulation system and the heat bath with the rate defined by a time constant $$\tau_T$$ 

### 3. Stochastic methods
Randomly assign a subset of atoms new velocities based on Maxwell-Boltzmann distributions for the target temperature. Randomization interferes with correlated motion and thus slows down the system's kinetics.

#### Andersen thermostat

- Assign a subset of atoms new velocities that are randomly selected from the Maxwell-Boltzmann distribution for the target temperature [[Andersen-1980][Andersen-1980]]. 
- "massive Andersen" thermostat randomizes the velocities of all atoms [[Basconi-2013][Basconi-2013]]. 

The Andersen thermostat: 
- Correctly samples the canonical ensemble
- Does not conserve momentum.  
- Can impair correlated motions and thus slow down the kinetics of the system. 
- Not recommended when studying kinetics or diffusion properties of the system. 

#### The Lowe-Andersen thermostat 
- A variant of the Andersen thermostat that conserves momentum [[Koopman-2006][Koopman-2006]]. 
- Perturbs the system dynamics to a far less than the original Andersen method. 
- Improves suppressed diffusion in the system relative to the original Andersen.

#### Bussi stochastic velocity rescaling thermostat
- Extension of the Berendsen method corrected for sampling the canonical distribution. 
- The velocities of all the particles are rescaled by a properly chosen random factor [[Bussi-2007][Bussi-2007]].

#### Langevin thermostat
- Mimics the viscous aspect of a solvent and interaction with the environment. 
- Adds a frictional force and a random force. 
- The frictional force and the random force combine to give the correct canonical ensemble.
- The amount of friction is controlled by the damping coefficient.

### 4. Extended system thermostats
#### Nosé-Hoover thermostat
- The heat bath is integrated with the system by addition of an artificial variable associated with a fictional "heat bath mass" to the equations of motion. 
- The temperature can be controlled without involving random numbers. 
- Correlated motions are not impaired
- Better description of kinetics and diffusion properties [Nose-1984][Nose-1984], [Hoover-1985][Hoover-1985].

##### Drawbacks:
 - Periodic temperature fluctuations with the frequency proportional to the "heat bath mass". 
 - Imparts the canonical distribution as well as ergodicity (space-filling). 

The time constant parameter in this thermostat controls the period of temperature fluctuations at equilibrium. 

#### Nosé-Hoover-chains
- A modification of the Nosé-Hoover thermostat which includes not a single thermostat variable but a chain of variables with different "masses" [Martyna-1992][Martyna-1992]. 
- Chaining variables with different masses helps to suppress oscillations. 

### Global and local thermostats
- Global thermostats control temperature of all atom in a system uniformly. 
- This may lead to cold solute and hot solvent due to a slow heat transfer.


- Local thermostats allow to control temperature in selected groups of atoms independently. 
- Local thermostats work well for large solutes.
- Temperature of small solutes this approach may significantly fluctuate leading to unrealistic dynamics.

<br>

> ## Challenge: Thermodynamic ensembles
>
> What thermodynamic ensemble describes an isolated system?
> 1. Canonical
> 2. Grand canonical
> 3. Isothermal-isobaric
> 4. Microcanonical 
>
> > ## Solution
> >
> > Microcanonical
> >
> {: .solution}
{: .challenge}

<br>

> ## Challenge: Thermostats
>
> Which of the following statements is correct? 
> 1. Conformational transitions are not affected by stochastic temperature control methods 
> 2. The Berendsen thermostat is very useful for heating simulation systems
> 3. Extended system thermostats control temperature by random velocity rescaling
> 4. Local thermostats work well for small groups of atoms
>
> > ## Solution
> >
> > The Berendsen thermostat is very useful for heating simulation systems
> >
> {: .solution}
{: .challenge}

>## Specifying local thermostats
>With NAMD it is possible to set coupling coefficients for each atom in occupancy or beta column of a pdb file:
>  
>langevinFile  
>langevinCol  
>
>tCoupleFile  
>tCoupleFCol
>
>In GROMACS temperature of a selected groups of atoms can be controlled independently using *tc-grps*.  Temperature coupling groups are coupled separately to temperature bath.
{: .callout}

>## Selecting thermostats in molecular dynamics packages
>
>| Thermostat/MD package | GROMACS                      |  NAMD                    | AMBER         |
>|-----------------------|------------------------------|--------------------------|---------------|
>| velocity rescaling    |                              |  reascaleFreq (steps)    |               |
>| velocity reassignment |                              |  reassignFreq (steps)    |               |
>| Andersen              | tcoupl = andersen            |                          |               |
>| massive-Andersen      | tcoupl = andersen-massive    |                          | ntt = 2       |
>| Lowe-Andersen         |                              |  loweAndersen on         |               |
>| Berendsen             | tcoupl = berendsen           |  tCouple on              | ntt = 1       |
>| Langevin              |                              |  langevin on             | ntt = 3       |
>| Bussi                 | tcoupl = V-rescale           |  stochRescale  on        |               |
>| Nose-Hoover           | tcoupl = nose-hoover         |                          |               |
>| Nose-Hoover-chains    | nh-chain-length (default 10) |                          |               |
>
{: .callout}

{% comment %}
### References
Below here we resolve reference-style links so that 
[Refname-YYYY] points to the anchor #Refname-YYYY on the {{ page.root }}/reference.html page.

Example:
[Refname-YYYY]: {{ page.root }}/reference.html#Refname-YYYY
{% endcomment %}

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
