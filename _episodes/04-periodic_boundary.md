---
title: "Periodic boundary conditions"
teaching: 30
exercises: 0
questions:
- "How to simlate a large system?"
- "Why the simulation system should be neutralized"
objectives:
- "Explain how to setup periodic boundary"
---

## Periodic boundary conditions
Periodic boundary conditions (PBC) are used to approximate a large system by using a small part called a unit cell. The boundary to contain molecules in simulation is needed to preserve thermodynamic properties like temperature, pressure and density.

To implement PBC the unit cell is surrounded by translated copies in all directions to approximate an infinitely large system. When one molecule diffuses across the boundary of the simulation box it reappears on the opposite side. So each molecule always interacts with its neighbours even though they may be on opposite sides of the simulation box. This approach replaces the surface artifacts caused by the interaction of the isolated system with a vacuum with the PBC artifacts which are in general much less severe.

In simulations with PBC the non-bonded interaction cut-off radius should be smaller than half the shortest periodic box vector to prevent interaction of an atom with its image.

## Balancing of charges
Neutralizing a system is a practice carried out for obtaining correct electrostatic energy during the simulation. This is done because under periodic boundary and using grid-based electrostatic the system has to be neutral. Otherwise, the electrostatic energy will essentially add to infinity from the interaction of the box with the infinite number of the periodic images. Simulation systems are most commonly neutralized by adding sodium or chloride ions.

> ## Specifying periodic box in GROMACS
> The box specification is integrated into structure file. The [editconf](http://manual.gromacs.org/archive/5.0/programs/gmx-editconf.html) utility is used to set the box parameters:
>
> **-bt**  Box type (triclinic, cubic, dodecahedron, octahedron)<br>
>
> **-box** Box vectors lengths (a,b,c)<br>
>
> **-angles** Box vectors angles   (bc,ac,ab)<br>
{: .callout}

> ## Specifying periodic box in NAMD
> Periodic box is defined by three unit cell vectors:
>
> **cellBasisVector1**
>>
>> Default value: 0 0 0
>
> **cellBasisVector2**
>
>> Default value: 0 0 0
>
> **cellBasisVector3**
>
>> Default value: 0 0 0
>
> **extendedSystem**
>> NAMD generates a .xsc (eXtended System Configuration) file which contains the periodic cell parameters. If this keyword is used periodic box parameters will be read from .xsc file ignoring cellBasisVectors.
>>Value: filename
{: .callout}
