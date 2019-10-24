---
title: "Fast Methods to Evaluate Forces"
teaching: 20
exercises: 0
questions:
- "Why the computation of nonbonded interactions is the speed limiting factor?"
- "How are nonbonded interactions computed efficiently?"
- "What is a cutoff distance?"
- "How to choose the appropriate cutoff distance?"
objectives:
- "Explain how nonbonded interactions are truncated."
- "Determine the appropriate cutoff distance and truncation method."
- "Set cutoff distance and truncation method in MD input files."
keypoints:
- Nonbonded interactions are truncated to speed up simulations of large systems.
- To ensure consistency and reproducibility of simulation you should choose the cutoff appropriate for the force field and the size of the periodic box.
---

The most computationally demanding part of a molecular dynamics simulation is the calculation of the nonbonded terms of the potential energy function. As non-bonded energy terms between every pair of atoms should be evaluated, the number of calculations increases as the square of the number of atoms. To speed up the computation, only the interactions between two atoms separated by a distance less than a pre-defined cutoff distance are evaluated.

## Neighbour Searching Methods
 The search for pairs of particles that are needed for calculation of the short-range nonbonded interactions is usually accelerated by maintaining a list of all particles within a predefined cutoff distance of each other.  Particle neighbours are determined either by dividing the simulation system into grid cells (cell lists) or by constructing a neighbour list for each particle (Verlet lists).

### Cell Lists
The cell lists method divides the simulation domain into *n* cells within edge length greater or equal to the cutoff radius of the interaction to be computed.  The interaction potential for each particle is then computed as the sum of the pairwise interactions between the particle and all other particles in the same cell and all other particles in the neighbouring cells (26 cells for 3-dimensional simulation).

### Verlet Lists
A Verlet list stores all particles within the cutoff distance of every particle plus some extra buffer distance. Although all pairwise distances must be evaluated to construct the Verlet list, it can be used for several consecutive time steps until any particle has moved more than half of the buffer distance. At this point, the list is invalidated and the new list must be constructed. Verlet offers more efficient computation of pairwise interactions at the expense of relatively large memory requirement which can be a limiting factor. In practice, almost all simulations are run in parallel and use a combination of spatial decomposition and Verlet lists.

## Problems with Truncation of Lennard-Jones Interactions and How to Avoid Them
We have learned that the LJ potential is always truncated at the cutoff distance. A cutoff introduces a discontinuity in the potential energy at the cutoff value. As forces are computed by differentiating potential, a sharp difference in potential may result in nearly infinite forces at the cutoff distance (Figure 1A). There are several approaches to minimize the impact of the cutoff.

![Cutoff Methods](../fig/Cutoff_Methods.svg)

### Shifted Potential
The standard solution is to shift the whole potential uniformly by adding a constant at values below cutoff (shifted potential method, Figure 1B). This ensures continuity of the potential at the cutoff distance and avoids infinite forces. The addition of a constant term does not change forces at the distances below cutoff because it disappears when the potential is differentiated. However, it introduces a discontinuity in the force at the cutoff distance. Particles will experience sudden unphysical acceleration when other particles cross their respective cutoff distance. Another drawback is that when potential is shifted the total potential energy changes.
### Shifted Forces
One way to address discontinuity in forces is to shift the whole force so that it vanishes at the cutoff distance (Figure 1C).  As opposed to the potential shift method the shifted forces cutoff modifies equations of motion at all distances. Nevertheless, the shifted forces method has been found to yield better results at shorter cutoff values compared to the potential shift method (Toxvaerd, 2011)
### Switching Function
Another solution is to modify the shape of the potential function near the cutoff boundary to truncate the non-bonded interaction smoothly at the cutoff distance. This can be achieved by the application of a switching function, for example, polynomial function of the distance. If the switching function is applied the switching parameter specifies the distance at which the switching function starts to modify the LJ potential to bring it to zero at the cutoff distance. The advantage is that the forces are modified only near the cutoff boundary and they approach zero smoothly.

### How to Choose the Appropriate Cutoff Distance?
A common practice is to truncate at 2.5 $$\sigma$$ and this practice has become a minimum standard for truncation.  At this distance, the LJ potential is about 1/60 of the well depth $$\epsilon$$, and it is assumed that errors arising from this truncation are small enough. The dependence of the cutoff on $$\sigma$$ means that the choice of the cutoff distance depends on the force field and atom types used in the simulation. For example for the O, N, C, S, and P atoms in the AMBER99 force field the values of $$\sigma$$ are in the range 1.7-2.1,  while for the Cs ions  $$\sigma=3.4$$. Thus the minimum acceptable cutoff, in this case, is 8.5.

In practice, increasing cutoff does not necessarily improve accuracy. There are documented cases showing opposite tendency [(Yonetani, 2006)]({{ page.root }}/reference.html#Yonetani-2006).  Each force field has been developed using a certain cutoff value, and effects of the truncation were compensated by adjustment of some other parameters. If you use cutoff 14 for the force field developed with the cutoff 9, then you cannot say that you used this forcefield. Thus to ensure consistency and reproducibility of simulation you should choose the cutoff appropriate for the force field.

> ## Cutoffs Used for Development of Common Force Fields
> AMBER: 9 <span>&#8491;</span>
>
> CHARMM: 12 <span>&#8491;</span>
>
> GROMOS: 14 <span>&#8491;</span>
>
> OPLS: 11-15 <span>&#8491;</span> (depending on a molecule size)
{: .callout}

Different molecular properties are affected differently by various cutoff approximations. Examples of properties that are very sensitive to the choice of cutoff include the surface tension [(Ahmed, 2010)]({{ page.root }}/reference.html#Ahmed-2010), the solid–liquid coexistence line [(Grosfils, 2009)]({{ page.root }}/reference.html#Grosfils-2009), the lipid bilayer properties [(Huang, 2014)]({{ page.root }}/reference.html#Huang-2014), and the structural properties of proteins [(Piana, 2012)]({{ page.root }}/reference.html#Piana-2012). For such quantities even a cutoff at 2.5 $$ \sigma $$ gives inaccurate results, and in some cases the cutoff must be larger than 6 $$ \sigma $$ was required for reliable simulations [(Grosfils, 2009)]({{ page.root }}/reference.html#Grosfils-2009).

Cutoff problems are especially pronounced when energy conservation is required. The net effect could be an increase in the temperature of the system over time. The best practice is to carry out trial simulations of the system under study without temperature scaling to test it for energy leaks or sources before a production run.


> ## Specifying Truncation of LJ Potential in GROMACS
> **vdw-modifier** Controls how to truncate LJ potential. Acceptable values:
> >**potential-shift** Shifts the Van der Waals potential by a constant such that it is zero at the **rvdw**.
>>
>>**force-switch** Smoothly switches the forces to zero between **rvdw-switch** and **rvdw**.
>>
>> **potential-switch** Smoothly switches the potential to zero between **rvdw-switch** and **rvdw**.
>>
>>**none**
>
> **rvdw-switch** Where to start switching.
>
> **rvdw**
> Cut-off distance
{: .callout}
> ## Specifying Truncation of LJ Potential in NAMD
> **cutoff** Cut-off distance common to both electrostatic and van der Waals calculations
>
> **switching** Turn switching on/off
>
> **switchdist** Where to start switching
>
> **vdwForceSwitching** Use force switching for VDW
>> Default Value: off
{: .callout}

## Truncation of the Electrostatic Interactions
Electrostatic interactions occurring over long distances are known to be important for biological molecules. Electrostatic interactions decay slowly and simple increase of the cutoff distance to account for long-range interactions can dramatically raise the computational cost. In periodic simulation systems, the most commonly used method for calculation of long-range electrostatic interactions is particle-mesh Ewald.  In this method, the electrostatic interaction is divided into two parts: a short-range contribution, and a long-range contribution. The short-range contribution is calculated by exact summation of all pairwise interactions of atoms separated by a distance that is less than cutoff in real space. The forces beyond the cutoff radius are approximated in Fourier space commonly by the Particle-Mesh Ewald (PME) method.


> ## Specifying Cutoff and Neighbour Searching Methods in GROMACS
> **cutoff-scheme**
> Since version 5.1 **group** list has been deprecated and only **Verlet** scheme is available
>
> **rlist**
> Cutoff distance for the short-range neighbour list. Active when **verlet-buffer-tolerance** = -1, otherwise ignored
>
> **verlet-buffer-tolerance**
> The maximum allowed error for pair interactions per particle caused by the Verlet buffer. To achieve the predefined tolerance the cutoff distance **rlist** is adjusted indirectly. To override this feature set the value to -1
>
>> Default value: 0.005 [kJ mol<sup>-1</sup> ps<sup>-1</sup>]
>
>**nstlist**
> Frequency to update the neighbour list. If set to 0 the neighbour list is constructed only once and never updated.
>
>> Default Value: 10
>
>**ns-type**
> Neighbour search method. Acceptable values:
>> **grid**: make a grid in the box and only check atoms in neighboring grid cells.
>>
>> **simple**: loop over every atom in the box.
{: .callout}

> ## Specifying Cutoff and Neighbour Searching Methods in NAMD
> When run in parallel NAMD uses a combination of spatial decomposition into grid cells, "patches" and Verlet lists with extended cutoff distance
>
> **cutoff**
> Local interaction distance. Same for both electrostatic and VDW interactions
>
> **pairlistdist**
> Distance between pairs for inclusion in pair lists. Should be bigger or equal than the **cutoff**.
>
>> Default value: **cutoff**
>
> **stepspercycle**
> Number of timesteps in each cycle. Each cycle represents the number of timesteps between atom reassignments.
>> Default Value: 20
>
>**pairlistsPerCycle**
> Specifies how many times per cycle to regenerate pairlists.
>> Default value: 2
{: .callout}
