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
- "The macromolecule shape, rotation and conformational changes should be taken into account in chosing the periodic box parameters."
---
Periodic boundary conditions (PBC) are used to approximate a large system by using a small part called a unit cell. The boundary to contain molecules in simulation is needed to preserve thermodynamic properties like temperature, pressure and density. Application of PBC to simulations allows to include the influence of bulk solvent or crystalline environments.

To implement PBC the unit cell is surrounded by translated copies in all directions to approximate an infinitely large system. When one molecule diffuses across the boundary of the simulation box it reappears on the opposite side. So each molecule always interacts with its neighbours even though they may be on opposite sides of the simulation box. This approach replaces the surface artifacts caused by the interaction of the isolated system with a vacuum with the PBC artifacts which are in general much less severe.

In simulations with PBC the non-bonded interaction cut-off radius should be smaller than half the shortest periodic box vector to prevent interaction of an atom with its image.

> ## Specifying Periodic Box
>  **GROMACS**
>
> The box specification is integrated into structure files. The box parameters can be set using the [*editconf*](http://manual.gromacs.org/archive/5.0/programs/gmx-editconf.html) or manually. The *editconf* program accepts the following options:
>
> *-bt* &emsp; Box type: *triclinic, cubic, dodecahedron, octahedron*<br>
>
> *-box* &emsp; Box vectors lengths *a, b, c* in nm.<br>
>
> *-angles* &emsp; Box vectors angles *bc, ac, ab* in degrees.
>
> *-d* &emsp; Distance between the solute and the box in nm.
>
>Example:
>~~~
>wget http://files.rcsb.org/view/1lyz.pdb
>gmx pdb2gmx -f 1lyz.pdb -ff amber99sb-ildn -water spce -ignh
>gmx editconf -f conf.gro -o conf_boxed.gro -d 1.0 -bt cubic
>~~~
> {: .bash}
> in the example above the *editconf* program will append box vectors to the structure file *'conf.gro'* and save it in the file *'conf_boxed.gro'*. The 9 components of the three box vectors are saved in the last line of the structure file in the order: xx yy zz xy xz yx yz zx zy. Three of the values (xy, xz, and yz) are always zeros because they are duplicates of (yx, zx, and zy).  The values of the box vectors components are related to the unit cell vectors $$a,b,c,\alpha,\beta,\gamma$$ from the *CRYST1* record of a PDB file with the equations:
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
> Periodic box is specified in the run parameter file by three unit cell vectors. the units are <span>&#8491;</span>.
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
{: .callout}

## What size/shape of a periodic box should I use?

Solvated macromolecules rotate during simulations. Furthermore macromolecules may undergo conformational changes. Often these changes are of major interest and shoud not be restricted in any way. If the molecule is not spherical and the box dimension is not large enough rotation will result in the interaction between copies. This artefactual interaction may influence the motions of the system and affect the outcome of the simulation. To avoid these problems the minimum box dimension should be larger than the largest dimension of the macromolecule plus at least 10 <span>&#8491;</span>.

 A cubic box is the most intuitive and common choice, but it is inefficient due to irrelevant water molecules in the corners. The extra water will make your simulation run slower. Ideally you need a sufficiently large sphere of water surrounding the macromolecule, but that's impossible because spheres can't be packed to fill space. A common alternatives that are closer to spherical are the dodecahedron or the truncated octahedron. These shapes work reasonably well for globular macromolecules, but if the solute is elongated there will be a large amount of the unnecessary water located far away from the solute. In this case you may consider constraining the rotational motion [[ref]](https://aip.scitation.org/doi/10.1063/1.480557) and using a smaller rectangular box. But be aware that the box shape itself may influence conformational dynamics by restricting motions in certain directions [[ref]](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20341). This effect may be significant when the amount of solvent is minimal.



 >## Comparing periodic boxes
 >Using the structure file *'conf.gro'* from the example above generate triclinic, cubic, dodecahedral and octahedral boxes with the 15<span>&#8491;</span> distance between the solute and the box edge.
>
>Compare the volumes of the boxes.
>Which of the boxes will be the fastest to simulate?
{: .challenge}

Truncated octahedron (implemented in AMBER)
V = x * y * z * sqrt(1.0 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 
        2.0 * cos(alpha) * cos(beta) * cos(gamma)) 

x, y, z are the box lengths: x = y = z = a in triclinic system. alpha, beta, 
gamma are the box angles which equal 109.47 (angles between each sides) 

Gromacs manual:
[Periodic box types](https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html?highlight=periodic%20boundary%20conditions)

Truncated octahedron periodic boundary conditions are isomorphic to 
parallelepiped boundary conditions with specific cell vectors.
If you view the simulation in VMD, it will look like a parallelepiped and not a 
truncated octahedron, but they have exactly the same properties, so it 
doesn't matter. 
