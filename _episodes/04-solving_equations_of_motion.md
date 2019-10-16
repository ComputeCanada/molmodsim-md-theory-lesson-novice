---
title: "Integrating the Equations of Motion"
teaching: 30
exercises: 0
questions:
- "How to advance simulation time?"
objectives:
- "Explain how to choose simulation timestep"
---

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
