---
layout: page
title: "Lesson Design"
permalink: /design/
---

> ## Help Wanted
> {:.no_toc}
>
> **We are filling in the exercises [below](#stage-3-learning-plan)
> in order to make the lesson plan more concrete.
> Contributions (both in the form of pull requests with filled-in exercises,
> and comments on specific exercises, ordering, and timings) are greatly appreciated.**
{: .callout}

## Motivation for the BST-Workshop series

The motivation behind this 3-part workshop/tutorial series is to give new grad 
students an introduction into preparing, running and analyzing MD simulations.
The series is not meant to be a replacement for university courses in thermodynamics 
and statistical thermodynamics.  Instead this series takes an applied approach
and touches on theory only as much it is necessary to aid choosing reasonable
simulation parameters.

The **first lesson** is a hands-on example on how to use GROMACS to prepare and 
submit a MD simulation.  There are numerous good GROMACS tutorials available and 
maybe one of them could be adapted for this.  In addition to general instructions 
and commands, there should be specific instructions on how to do this on Compute 
Canada clusters.  As in Software Carpentry lessons, certain steps should be 
implemented as exercises for the learner - the correct answer can be hidden in 
an answer box.  
Theory should be skipped over as much as possible and should be covered in the 
second module.  
This lesson could be adapted for other MD-codes (e.g. NAMD) as well.


The **second lesson** is a small theory review to remind the learners of important 
concepts and how they influence the choice of simulation parameters.  The aim is
that people can avoid pit-falls and misunderstandings that are commonly made
by novice MD users. Topics to cover are Periodic Boundary Conditions (and why 
there is no "outside" with a periodic box), thermostats/barostats, cut-offs, etc.  
This should be no substitute for a formal course in statistical thermodynamics 
but should help users make more informed choices of simulation settings/parameters.


The **third lesson** is again a hands-on tutorial to use Python and the Python 
packages [MDAnalysis](https://www,mdanalysis.org), [MDtraj](http://mdtraj.org/),
and [NGLview](http://nglviewer.org/nglview/latest/) to write their own analysis 
tools.  While, for example GROMACS, comes with a very large number of analysis 
tools out-of-the-box, users often limit themselves to those tools and the options
and variations that they offer.  The above frameworks make it very easy to read 
MD-trajectories in different formats and get access to the coordinates and come 
up with fully customized analysis methods.


The mentioned lessons *one* and *three* are currently put on hold and are subject
to be created at a later time.  In the meantime we can recommend to work through
online tutorials that are already available such as:
* http://www.mdtutorials.com/ 
* https://www.mdanalysis.org/MDAnalysisTutorial/

This module here will implement the *second lesson* focusing on theory and giving
guidance on choosing good parameters for MD-Simulations.


## Process Used

* 

## Stage 1: Assumptions

*   Audience
    * Graduate students in with some background in chemistry (e.g. Biochemistry, 
      Chemical Engineering, Condensed Matter Physics or similar),
      who are beginning a research project with Molecular Dynamics.
    * Had a course in thermodynamics and possibly statistical thermodynamics.
    * Have either attended lesson 1 of this series or worked through a MD 
      tutorial or two.  
      At this point they are likely overwhelmed by the number of simulation 
      parameters that need to be set and need some practical guidance.

*   Constraints
    * A half day or maybe full day?
    * Difficult to make this interactive.  
        * One could have the students run simulations with bad settings, 
          which then either are unstable and fail, result on obviously bad 
          results or compare observables with simulations with better settings.
        * This however will take a lot longer than a traditional lecture.
        * One could prepare pairs of simulations with bad and good settings
          and compare the results.

*   Motivating Example
    * 

*   Data
    * 

## Stage 2: Desired Results

### Topics to cover 
(in no particular order; might need to be pruned)

* System preparation
    * complete input structure (check for missing atoms/residues)
    * balance charges (the world does not have net-charge)
    * non-integer charges point to broken topology
    * interfaces (?)
    * membranes (?)

* Periodic Boundary Conditions
    * there is no "outside of the box" with PBC
    * minimum image convention
    * cut-off restrictions

* Force fields and Cut-Offs
    * Genealogy of common Forcefields (?)
    * bonded interactions
    * non-bonded interactions
        * cut-offs
        * Electrostatics
        * long range interactions (PME)
    * Water models (?)

* The Global MD algorithm
    * choice of time steps
    * applying constraints
    * (but probably not going into differences between Leap-Frog and Velocity-Verlet)

* Phases of MD Workflow
    * energy minimization
    * relaxation / position restrained MD
    * equilibration
    * production MD

* Thermostats/Barostats
    * Maxwell-Boltzmann distribution
    * Why use different thermo-/barostats for equilibration and production.
        * Berendsen 's algorithms are suitable for equilibration (exponential 
          approach to target temperature/pressure) but result in a non-physical
          ensemble.
        * Parrinello-Rahman and Nose-Hoover are suitable for production MD
          (correct physical ensemble), but approach the target temperature/pressure
          in a dampened oscillation, which might be unstable and takes much longer.
    * Size of temperature coupling groups
    * Flying ice cube -> resetting COM-movement
    * Hot-Solvent/Cold-Solute -> tcoupl groups
    * Literature:
        * [Basconi2013]({{ page.root }}/reference.html#Basconi-2013)
        * [Wong-ekkabut2016]({{ page.root }}/reference.html#Wong-ekkabut-2016)

* Parallelization/Performance
    * Particle decomposition vs. Domain Decomposition
    * MPI / OpenMP / GPU
    * PME nodes /  shifting of real-space/PME cutoff
    * Load balancing
    * Literature:
        * [Larsson2011]({{ page.root }}/reference.html#Larsson-2011)

### Questions

How do I...

  * 

### Skills

I can...

  * 

### Concepts

I know...

  * 

## Stage 3: Learning Plan

### Summative Assessment

*   Midpoint: 
    * show input-files with mistakes and let users spot them
*   Final:

### Lesson 1:
