---
title: "Advancing Simulation in Time"
teaching: 20
exercises: 0
questions:
- "How is simulation time advanced in a molecular dynamics simulation?"
- "What factors are limiting a simulation time step?"
- "How to accelerate a simulation?"
objectives:
- "Understand how simulation is advanced in time."
- "Learn how to choose time parameters and constraints for an efficient and accurate simulation."
- "Learn how to specify time parameters and constraints in GROMACS and NAMD."
keypoints:
- "A good integration algorithm for MD should be time-reversible and energy conserving."
- "The most widely used integration method is the velocity Verlet."
- "Simulation time step must be short enough to describe the fastest motion."
- "Time step can be increased if bonds involving hydrogens are constrained."
- "Additional time step increase can be achieved by constraining all bonds and angles involving hydrogens."
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

The Velocity Verlet algorithm is mathematically equivalent to the original Verlet algorithm. It explicitly incorporates velocity, solving the problem of the first time step in the basic Verlet algorithm.
- *Due to its simplicity and stability the Velocity Verlet has become the most widely used algorithm in the MD simulations.*


### The Leap Frog Algorithm
 The leap frog algorithm is a modified version of the Verlet algorithm. Using accelerations of the current time step, compute the velocities at half-time step:

$\vec{v}(t+\frac{1}{2}\delta+t)=\vec{v}(t-\frac{1}{2}\delta{t})\cdot\delta{t}+\vec{a}(t)\cdot\delta{t}$

Then determine positions at the next time step:

$\vec{r}(t+\delta t)=\vec{r}(t)+\vec{v}(t+\frac{1}{2}\delta{t}))\cdot\delta{t}$

The Leap Frog algorithm is essentially the same as the Velocity Verlet. The Leap Frog and the Velocity Verlet integrators give equivalent trajectories. The only difference is that the velocities are not calculated at the same time as positions. Leapfrog integration is equivalent to updating positions and velocities at interleaved time points, staggered in such a way that they "leapfrog" over each other.


> ## Selecting the Intergator
> **GROMACS**
>
>Several integration algorithms available in GROMACS are specified in the run parameter **mdp** file.
> ~~~
> integrator = md
>; A leap frog algorithm
>
>integrator = md-vv
>;  A velocity Verlet algorithm
>
>integrator = md-vv-avek
>; A velocity Verlet algorithm same as md-vv except the kinetic energy is calculated as the average of the two half step kinetic energies. More accurate than the md-vv.
>
>integrator = sd
>;  An accurate leap frog stochastic dynamics integrator.
>
>integrator = bd
>; A Euler integrator for Brownian or position Langevin dynamics.
> ~~~
> {: .file-content}
> **NAMD**
>
>The only available integration method is Verlet.
{: .callout}


## How to Choose Simulation Time Step?
Mathematically Verlet family integrators are stable for time steps

$$\delta{t}\leq\frac{2}{w}$$ where $$\omega$$ is angular frequency.

In molecular dynamics stretching of the bonds with the lightest atom H is usually the fastest motion. The period of oscillation of a C-H bond is about 10 fs. Hence Verlet integration will be stable for time steps < 3.2 fs. In practice, the time step of 1 fs is recommended to describe this motion reliably.

If the dynamics of hydrogen atoms is not essential for a simulation, bonds with hydrogens can be constrained. By replacing bond vibrations with holonomic (not changing in time) constraints the simulation step can be doubled since the next fastest motions (bond vibrations involving only heavy atoms and angles involving hydrogen atoms) have a period of about 20 fs. Further increase of the simulation step requires constraining bonds between all atoms and angles involving hydrogen atoms. Then the next fastest bond vibration will have 45 fs period allowing for another doubling of the simulation step.

To accelerate a simulation the electrostatic interactions outside of a specified cutoff distance can be computed less often than the short range bonded and non-bonded interactions. It is also possible to employ an intermediate timestep for the short-range non-bonded interactions, performing only bonded interactions every timestep.

> ## Specifying Time Parameters
> **GROMACS**
>
> Time parameters are specified in the **mdp** run parameter file.
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
> {: .file-content}
> **NAMD**
>
> Time parameters are specified in the **mdin** run parameter file.
> ~~~
> TimeStep = 1
># Time step, fs
>
> NumSteps = 10000
># Number of steps to simulate
>
>FirstTimeStep = 0
># Time of the first step
>
>#  Multiple time stepping parameters
>nonbondedFreq 2
>#  Number of timesteps between short-range non-bonded evaluation.
>
>fullElectFrequency 4
># Number of timesteps between full electrostatic evaluations
>~~~
> {: .file-content}
{: .callout}

### Constraint Algorithms
To constrain bond length in a simulation the equations of motion must be modified. This is often accomplished by the application of constraint forces acting along a bond in opposite directions. The total energy of the simulation system is not affected in this case because the total work done by constraint forces is zero. In constrained simulation first the unconstrained step is done, then corrections are applied to satisfy constraints.

Because bonds in molecules are coupled satisfying all constraints is a non-linear problem. Is it fairly easy to solve it for a small molecule like water but as the number of coupled bonds increases, the problem becomes more difficult. Several algorithms have been developed for use specifically with small or large molecules.

**SETTLE** is very fast analytical solution for small molecules. It is widely used to constrain bonds in water molecules.

**SHAKE** is an iterative algorithm that resets all bonds to the constrained values sequentially until the desired tolerance is achieved. SHAKE is simple and stable, it can be applied for large molecules and it works with both bond and angle constraints. However it is substantially slower than SETTLE and hard to parallelize. SHAKE may fail to find the constrained positions when displacements are large. The original SHAKE algorithm was developed for use with a leap-frog integrator.  Later on, the extension of SHAKE for use with a velocity Verlet integrator called RATTLE has been developed. Several other extensions of the original SHAKE algorithm exist (QSHAKE, WIGGLE, MSHAKE, P-SHAKE).

**LINCS**  algorithm (linear constraint solver),  employs a power series expansion to determine how to move the atoms such that all constraints are satisfied.  It is 3-4 times faster than SHAKE and easy to parallelize. The parallel LINCS (P-LINKS) allows to constrain all bonds in large molecules. The drawback is that it is not suitable for constraining both bonds and angles.

> ## Specifying Constraints
> **GROMACS**
>
>SHAKE, LINKS and SETTLE constraint algorithms are implemented. They are selected via keywords in mdp input files
> ~~~
>constraints = h-bonds
>; Constrain bonds with hydrogen atoms
>
>constraints = all-bonds
>; Constrain all bonds
>
>constraints = h-angles
>; Constrain all bonds and additionally the angles that involve hydrogen atoms
>
>constraints = all-angles
>; Constrain all bonds and angles
>
>constraint-algorithm = LINKS
>; Use LINKS
>
>constraint-algorithm = SHAKE
>; Use SHAKE
>
>shake-tol = 0.0001
>;  Relative tolerance for SHAKE, default value is 0.0001.
> ~~~
> {: .file-content}
>SETTLE can be selected in the topology file:
>~~~
>[ settles ]
>; OW    funct   doh     dhh
>1       1       0.1     0.16333
>~~~
> {: .file-content}
> **NAMD**
>
>SHAKE and SETTLE constraint algorithms are implemented. They are selected via keywords in simulation input file.
> ~~~
>rigidBonds water
># Use SHAKE to constrain bonds with hydrogens in water molecules.
>
> rigidBonds all
># Use SHAKE to constrain bonds with hydrogens in all molecules.
>
>rigidBonds none
># Do not constrain bonds. This is the default.
>
>rigidTolerance 1.0e-8
># Stop iterations when all constrained bonds differ from the nominal bond length by less than this amount. Default value is 1.0e-8.
>
>rigidIterations 100
># The maximum number of iterations. If the bond lengths do not converge, a warning message is emitted. Default value is 100.
>
>rigidDieOnError on
># Exit and report an error if rigidTolerance is not achieved after rigidIterations. The deault value is on.
>
>useSettle on
># If rigidBonds are enabled then use the SETTLE algorithm to constrain waters. The default value is on.
> ~~~
> {: .file-content}
{: .callout}
