---
title: "Advancing Simulation in Time"
teaching: 10
exercises: 0
questions:
- "How is simulation time advanced?"
- "How to choose an appropriate simulation time step?"
objectives:
- "Explain how to choose simulation timestep and integration method"
keypoints:
- "Simulation time step must be short enough to describe the fastest motion"
- "The most widely used integration method is the velocity Verlet"
---
To simulate evolution of the system in time the integration algorithm advances positions of all atomistic by a small step $$\delta{t}$$ during which the forces are considered constant. If the time step is small enough the trajectory will be reasonably accurate. A good integration algorithm for MD should be time-reversible and energy conserving.

## Integration Algorithms
### The Euler Algorithm
The Euler algorithm uses the second order Taylor expansion to estimate position and velocity at the next time step:

$\vec{r}(t+\delta{t})=\vec{r}(t)+\vec{v}(t)\delta{t}+\frac{1}{2}a(t)\delta{t}^2$

$\vec{v}(t+\delta{t})=\vec{v}(t)+\frac{1}{2}a(t)\delta{t}$

The Euler algorithm is neither time-reversible nor energy conserving and hence rather unfavourable. Nevertheless, the Euler scheme can be used to integrate some other than classical MD equations of motion. For example, GROMACS offers a Euler integrator for Brownian or position Langevin dynamics.

### The Verlet Algorithm
Using the current positions and forces and the previous positions calculate the positions at the next time step:

$\vec{r}(t+\delta{t})=2\vec{r}(t)-\vec{r}(t-\delta{t})+a(t)\delta{t}^2$

The Verlet algorithm  [(Verlet, 1967)]({{ page.root }}/reference.html#Verlet-1967) requires positions at two time steps. It is inconvenient when starting a simulation. While velocities are not needed to compute trajectories, they are useful for calculating observables e.g. the kinetic energy. The velocities can only be computed once the next positions are calculated:

$\vec{v}(t+\delta{t})=\frac{r{(t+\delta{t})-r(t-\delta{t})}}{2\delta{t}}$

The Verlet algorithm is time-reversible and energy conserving.

### The Velocity Verlet Algorithm
The velocities, positions and forces are calculated at the same time according to:

$\vec{r}(t+\delta{t})=\vec{r}(t)+\vec{v}(t)\delta{t}+\frac{1}{2}a(t)\delta{t}$

$\vec{v}(t+\delta{t})=\vec{v}(t)+\frac{1}{2}[a(t)+a(t+\delta{t})]\delta{t}$

The Velocity Verlet algorithm is mathematically equivalent to the original Verlet algorithm. It explicitly incorporates velocity, solving the problem of the first time step in the basic Verlet algorithm. Due to its simplicity and stability is has become the most widely used algorithm in the MD simulations.


### The Leap Frog Algorithm
 The leap frog algorithm is a modified version of the Verlet algorithm. Using accelerations of the current time step, compute the velocities at half-time step:

$\vec{v}(t+\frac{1}{2}\delta+t)=\vec{v}(t-\frac{1}{2}\delta{t})\cdot\delta{t}+\vec{a}(t)\cdot\delta{t}$

Then determine positions at the next time step:

$\vec{r}(t+\delta t)=\vec{r}(t)+\vec{v}(t+\frac{1}{2}\delta{t}))\cdot\delta{t}$

The Leap Frog algorithm is essentially the same as the Velocity Verlet. The Leap Frog and the Velocity Verlet integrators give equivalent trajectories. The only difference is that the velocities are not calculated at the same time as positions. Leapfrog integration is equivalent to updating positions and velocities at interleaved time points, staggered in such a way that they "leapfrog" over each other.

## Choosing Time Step
Mathematically Vertet family integrators are stable for time steps

$$\delta{t}\leq\frac{2}{w}$$ where $$\omega$$ is angular frequency.

In molecular dynamics stretching of the bonds with the lightest atom H is usually the fastest motion. The period of oscillation of a C-H bond is about 10 fs. Hence Verlet integration will be stable for time steps < 3.2 fs. In practice, the time step of 1 fs is recommended to describe this motion reliably. If the dynamics of hydrogen atoms is not essential for a simulation, bonds with hydrogens can be constrained, and time step increased to 2 fs.

On using a too large integration time step in molecular dynamics simulations of coarse-grained molecular models [(Winger, 2009)]({{ page.root }}/reference.html#Winger-2009).

> ## Specifying Time Parameters
> **GROMACS**
>
> Time parameters are specified in the run parameter file **.mdp**
>~~~
> dt = 0.001
>; Time step, ps
>
> nsteps = 10000
>; Number of steps to simulate
>
> tinit = 0
>; Time of the first step
> ~~~
> {: .source}
> **NAMD**
>
> Time parameters are specified in the run parameter file **mdin**
> ~~~
> TimeStep = 1
># Time step, fs
>
> NumSteps = 10000
># Number of steps to simulate
>
>FirstTimeStep = 0
># Time of the first step
> ~~~
> {: .source}
{: .callout}

> ## Specifying Integration Method in GROMACS
> **integrator** Selects integration algorithm. Acceptable values:
>
>> **md** A leap frog algorithm
>>
>> **md-vv** A velocity Verlet algorithm
>>
>> **md-vv-avek** A velocity Verlet algorithm same as **md-vv** except the kinetic energy is calculated as the average of the two half step kinetic energies. It is more accurate than the md-vv
>>
>> **sd** an accurate leap frog stochastic dynamics integrator.
>>
>> **bd** a Euler integrator for Brownian or position Langevin dynamics.
{: .callout}

> ## Specifying Integration Method in  NAMD
>The only available integration method is Verlet. To further reduce the cost of computing short-range nonbonded interactions and full electrostatics, NAMD uses a multiple time-stepping integration scheme controlled by the following keywords:
>
> **nonbondedFreq** Number of timesteps between nonbonded evaluation<br>
>
> **fullElectFrequency** Number of timesteps between full electrostatic evaluations
{: .callout}


## Constraints
The bonds involving hydrogen atoms are often constrained in MD simulations. By replacing bond vibrations with holonomic (not changing in time) constraints the simulation step can be doubled since the next fastest motions (bond vibrations involving only heavy atoms and angles involving hydrogen atoms) have a period of about 20 fs. Further increase of the simulation step requires constraining bonds between all atoms and angles involving hydrogen atoms. Then the next fastest bond vibration will have 45 fs period allowing for another doubling of the simulation step.

To constrain bond length in a simulation the equations of motion must be modified. This is often accomplished by the application of constraint forces acting along a bond in opposite directions. The total energy of the simulation system is not affected in this case because the total work done by constraint forces is zero. In constrained simulation first the unconstrained step is done, then corrections are applied to satisfy constraints.

Because bonds in molecules are coupled satisfying all constraints is a non-linear problem. Is it fairly easy to solve it for a small molecule like water but as the number of coupled bonds increases, the problem becomes more difficult. Different algorithms have been developed for use specifically with small or large molecules.

**SETTLE** is very fast analytical solution for small molecules. It is widely used to constrain bonds in water molecules.

**SHAKE** is an iterative algorithm that resets all bonds to the constrained values sequentially until the desired tolerance is achieved. SHAKE is simple and stable, it can be applied for large molecules and it works with both bond and angle constraints. However it is substantially slower than SETTLE and hard to parallelize. SHAKE may fail to find the constrained positions when displacements are large. The original SHAKE algorithm was developed for use with a leap-frog integrator.  Later on, the extension of SHAKE for use with a velocity Verlet integrator called RATTLE has been developed. Several other extensions of the original SHAKE algorithm exist (QSHAKE, WIGGLE, MSHAKE, P-SHAKE).

**LINCS**  algorithm (linear constraint solver),  employs a power series expansion to determine how to move the atoms such that all constraints are satisfied.  It is 3-4 times faster than SHAKE and easy to parallelize. The parallel LINCS (P-LINKS) allows to constrain all bonds in large molecules. The drawback is that it is not suitable for constraining both bonds and angles.
