---
title: "Periodic Boundary Conditions"
teaching: 10
exercises: 0
questions:
- "How to simlate a large system?"
- "Why the simulation system should be neutralized"
objectives:
- "Explain why and when periodic boundary conditions are used"
- "Explain how to setup periodic box"
- "Explain why it is necessary to neutlalize the simulation system"
keypoints:
- "Periodic boundary conditions are used to approximate an infinitely large system"
-  "Simulation system should be neutralized by adding counterions"
---
Periodic boundary conditions (PBC) are used to approximate a large system by using a small part called a unit cell. The boundary to contain molecules in simulation is needed to preserve thermodynamic properties like temperature, pressure and density. Application of PBC to simulations allows to include the influence of bulk solvent or crystalline environments.

To implement PBC the unit cell is surrounded by translated copies in all directions to approximate an infinitely large system. When one molecule diffuses across the boundary of the simulation box it reappears on the opposite side. So each molecule always interacts with its neighbours even though they may be on opposite sides of the simulation box. This approach replaces the surface artifacts caused by the interaction of the isolated system with a vacuum with the PBC artifacts which are in general much less severe.

In simulations with PBC the non-bonded interaction cut-off radius should be smaller than half the shortest periodic box vector to prevent interaction of an atom with its image.

> ## Specifying periodic box
>  **GROMACS**
>
> The box specification is integrated into structure file. The box parameters can be set using the [editconf](http://manual.gromacs.org/archive/5.0/programs/gmx-editconf.html) or manually. The **editconf** program accepts the following options:
>
> **-bt**  Box type. Acceptable values: **triclinic, cubic, dodecahedron, octahedron**<br>
>
> **-box** Box vectors lengths **a,b,c** in nm.<br>
>
> **-angles** Box vectors angles **bc,ac,ab** in degrees.
>
> **-d** Distance between the solute and the box in nm.
>
>~~~
> gmx editconf -f system.gro -o system_wbox.gro -d 1.0 -bt cubic
>~~~
> {: .source}
> The **editconf** program appends box vectors to the structure (**.gro**) file. The 9 components of the three box vectors are saved in the last line of the structure file in the order: xx yy zz xy xz yx yz zx zy. Three of the values (xy, xz, and yz) are always zeros because they are duplicates of (yx, zx, and zy).  The values of the box vectors components are related to the unit cell vectors $$a,b,c,\alpha,\beta,\gamma$$ from the **CRYST1** record of a PDB file with the equations:
>
>$$xx=a, yy=b\cdot\sin(\gamma), zz=\frac{v}{(a*b*\sin(\gamma))}$$
>
>$$xy=0, xz=0, yx=b\cdot\cos(\gamma)$$
>
>$$yz=0, zx=c\cdot\cos(\beta), zy=\frac{c}{\sin(\gamma)}\cdot(cos(\alpha)-cos(\beta)\cdot\cos(\gamma))$$
>
>$$v=\sqrt{1-\cos^2(\alpha)-cos^2(\beta)-\cos^2(\gamma) +2.0\cdot\cos(\alpha)\cdot\cos(\beta)\cdot\cos(\gamma)}\cdot{a}\cdot{b}\cdot{c}$$
>
> **NAMD**
>
> Periodic box is specified in the run parameter file **mdin** by three unit cell vectors:
>
> **cellBasisVector1**, **cellBasisVector2**, **cellBasisVector3** in <span>&#8491;</span>.
>~~~
> # cubic box
> cellBasisVector1 100 0 0
> cellBasisVector2 0 100 0
> cellBasisVector3 0 0 100
>~~~
>{: .source}
> Alternatively periodic box parameters can be read from the **.xsc** (eXtended System Configuration) file by using the **extendedSystem** keyword.  If this keyword is used **cellBasisVectors** are ignored.  NAMD always generates  **.xsc** files at runtime.
>~~~
> extendedSystem restart.xsc
>~~~
>{: .source}
{: .callout}

## Balancing of charges
Neutralizing a system is a practice carried out for obtaining correct electrostatic energy during the simulation. This is done because under periodic boundary and using grid-based electrostatic the system has to be neutral. Otherwise, the electrostatic energy will essentially add to infinity from the interaction of the box with the infinite number of the periodic images. Simulation systems are most commonly neutralized by adding sodium or chloride ions.

gmx genion -s file.tpr -p cg-something.top -o output.gro -neutral.
