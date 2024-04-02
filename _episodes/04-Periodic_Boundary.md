---
title: "Periodic Boundary Conditions"
teaching: 20
exercises: 0
questions:
- "How to simulate a bulk (gas, liquid or solid) system by using only a small part?"
objectives:
- "Understand why and when periodic boundary conditions are used."
- "Understand how shape and size of a periodic box can affect simulation."
- "Learn how to set periodic box parameters in GROMACS and NAMD."
keypoints:
- "Periodic boundary conditions are used to approximate an infinitely large system."
- "Periodic box should not restrict molecular motions in any way."
- "The macromolecule shape, rotation and conformational changes should be taken into account in choosing the periodic box parameters."
---
### What is PBC and why is it important?

- In most cases we want to simulate a system in realistic environment, such as solution.
- Try simulating a droplet of water, it will simply evaporate.
- We need a boundary to contain water and control temperature, pressure, and density. 
- Periodic boundary conditions allow to approximate an infinite system by using a small part (unit cell). 

![Figure: Periodic Boundary Conditions]({{ page.root }}/fig/periodic_boundary-3.svg){: width="200" }

- Unit cell is surrounded by an infinite number of translated copies in all directions (images). 
- When a particle in unit cell moves across the boundary it reappears on the opposite side. 
- Each molecule always interacts with its neighbors even though they may be on opposite sides of the simulation box. 
- Artifacts caused by the interaction of the isolated system with a vacuum are replaced with the PBC artifacts which are in general much less severe.

## Choosing periodic box size and shape.
### Box shape
#### Cubic periodic box

- A cubic box is the most intuitive and common choice
- Cubic box is inefficient due to irrelevant water molecules in the corners. 

 ![Figure: Cubic box]({{ page.root }}/fig/cubic_box.svg){: width="200" }
 
- Ideal simulation system is a sphere of water surrounding the macromolecule, but spheres can't be packed to fill space.

#### Octahedral and dodecahedral periodic boxes

- The dodecahedron (12 faces) or the truncated octahedron (14 faces) are closer to sphere.

|                              | Space filling with truncated octahedrons |
|:----------------------------:|:-----------------------------------------|
| ![Figure: truncated Octahedron]({{ page.root }}/fig/trunc-octa.svg){: width="64" } | ![]({{ page.root }}/fig/truncated_octahedron_group.svg){: width="140" } |

- These shapes work reasonably well for globular macromolecules.

#### Triclinic periodic boxes
- Any repeating shape that fills all of space has an equivalent triclinic unit cell.
- A periodic box of any shape can be represented by a triclinic box with specific box vectors and angles.

![Figure: Triclinic Cell]({{ page.root }}/fig/triclinic_cell.gif){: width="200" }

- The optimal triclinic cell has about 71% the volume of the optimal rectangular cell.

###  Box size
- The minimum box size should extend at least 10 nm from the solute.

![10 nm margin]({{ page.root }}/fig/box_size-2.svg){: width="200" }

 - The shortest periodic box vector should be at least twice bigger than the cuf-off radius.  

![Figure: Periodic Boundary Conditions]({{ page.root }}/fig/periodic_boundary-4.svg){: width="300" }

- In simulations with macromolecules solvent molecules should not "feel" both sides of a solute.

![]({{ page.root }}/fig/box_size.svg){: width="400" }


#### Pitfalls
 - A simulation system with elongated solute in cubic or dodecahedral box  will have a large amount of water located far away from the solute.
 - Consider using a narrow rectangular box. 
 - Rotation of elongated macromolecules and/or conformational changes must be taken in consideration.
 
![Figure: Periodic Boundary Conditions]({{ page.root }}/fig/PBC.svg){: width="300" }

- Constrain the rotational motion. 
- The box shape itself may influence conformational dynamics by restricting motions in certain directions [[1]](https://aip.scitation.org/doi/10.1063/1.480557), [[2]](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20341).


> ## Specifying Periodic Box
>  **GROMACS**
>
> The box specification is integrated into structure files. The box parameters can be set using the [*editconf*](http://manual.gromacs.org/archive/5.0/programs/gmx-editconf.html) program or manually. The *editconf* program accepts the following options:
>
>-----------------|------------------|
> *-bt* &emsp;    | Box type            | *triclinic, cubic, dodecahedron, octahedron*
> *-box* &emsp;   | Box vectors lengths, *a, b, c* | nm
> *-angles* &emsp;| Box vectors angles, *bc, ac, ab* | degrees
> *-d* &emsp;     | Distance between the solute and the box | nm
>
>Example:
>~~~
>module load StdEnv/2020 gcc gromacs
>wget http://files.rcsb.org/view/1lyz.pdb
>gmx pdb2gmx -f 1lyz.pdb -ff amber99sb-ildn -water spce -ignh
>gmx editconf -f conf.gro -o conf_boxed.gro -d 1.0 -bt cubic
>~~~
> {: .language-bash}
> In the example above the *editconf* program will append box vectors to the structure file *'conf.gro'* and save it in the file *'conf_boxed.gro'*. The 9 components of the three box vectors are saved in the last line of the structure file in the order: xx yy zz xy xz yx yz zx zy. Three of the values (xy, xz, and yz) are always zeros because they are duplicates of (yx, zx, and zy).  The values of the box vectors components are related to the unit cell vectors $$a,b,c,\alpha,\beta,\gamma$$ from the *CRYST1* record of a PDB file with the equations:
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
> Periodic box is specified in the run parameter file by three unit cell vectors, the units are <span>&#8491;</span>.
>~~~
> # cubic box
> cellBasisVector1 100 0 0
> cellBasisVector2 0 100 0
> cellBasisVector3 0 0 100
>~~~
>{: .file-content}
> Alternatively periodic box parameters can be read from the *.xsc* (eXtended System Configuration) file by using the *extendedSystem* keyword.  If this keyword is used *cellBasisVectors* are ignored.  NAMD always generates  *.xsc* files at runtime.
>~~~
> extendedSystem restart.xsc
>~~~
>{: .file-content}
{: .callout .self_study_text }

>## Comparing periodic boxes
>Using the structure file *'conf.gro'* from the example above generate triclinic, cubic, dodecahedral and truncated octahedral boxes with the 15 <span>&#8491;</span> distance between the solute and the box edge.
>
>Which of the boxes will be the fastest to simulate?
{: .challenge }

References:  
1. [Molecular dynamics simulations with constrained roto-translational motions: Theoretical basis and statistical mechanical consistency](https://aip.scitation.org/doi/10.1063/1.480557) 
2. [The effect of box shape on the dynamic properties of proteins simulated under periodic boundary conditions](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20341)
3. [Periodic box types in Gromacs manual](https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html?highlight=periodic%20boundary%20conditions)
