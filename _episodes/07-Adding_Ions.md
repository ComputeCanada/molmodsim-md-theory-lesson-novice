---
title: "Solvating a System and Adding Ions"
teaching: 30
exercises: 0
questions:
- "Why the simulation system should be neutralized"
- "How to add ions?"
objectives:
- "Explain why it is necessary to neutlalize the simulation system"
keypoints:
- "Simulation system must be neutralized by adding counterions"
---
There are two main reasons to add ions to a simulation system:

I. Under periodic boundary conditions and using grid-based methods to compute Coulomb energy the box interacts with the infinite number of its periodic images. If simulation system is charged the electrostatic energy will essentially add to infinity. We need to neutralize the system by adding counterions to obtain the correct electrostatic energy during the simulation.

II. The conformations, dynamics and function of biological macromolecules are sensitive to salt concentration and composition of the local environment.


#### How Many Ions Do I Need to Add?

Fist we should add enough counterions to neutralize the system. The neutralized system  would represent a salt-free solution.

To mimic the macroscopic salt concentration in the environment being studied we will need to add more ions to a simulation system. The number of ion pairs can be estimated using the formula:

$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$

The problem with this approach is that it does not take into account the charge of macromolecules. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high. For more accurate salt concentration you can calculate the number of ions corrected for screening effects in biomolecule using [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) web server.


- Add ions to neutralize proteins
- Solvate the system
- Add more ions to simulate salt randomly


Solvation.
Standard practice is to tile a pre-equilibrated solvent box across the system and eliminate solvent molecules which clash with the solute.


The addions leap commands adds ions to simulation cells near the minima of the solute's electrostatic potential field. We will use it prior to solvation so that the potential from unequilibrated water does not interfere with it.



~~~
module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.1.1 amber/18
tleap

source leaprc.protein.ff14SB
source leaprc.water.spce
s = loadpdb 1RGG_chain_A_prot.pdb
loadoff spcebox.off
addions s Na+ 0
solvatebox s SPCBOX 15
addionsrand s Na+ 10 Cl- 10
~~~






Let's apply the first method to neutralize the protein. The **genion** program can automatically determine the number of ions required to neutralize a system when the option **-neutral** us used:
~~~
$ echo "SOL" | gmx genion -s solvated.tpr -p topol.top -neutral -o neutralized.pdb
~~~
{: .bash}
~~~
...
Processing topology
Replacing 2 solute molecules in topology file (topol.top)  by 0 NA and 2 CL ions.

Back Off! I just backed up topol.top to ./#topol.top.4#
Using random seed 2093762991.
Replacing solvent molecule 3947 (atom 14303) with CL
Replacing solvent molecule 2365 (atom 9557) with CL
...
~~~
{: .output}
The **genion** program replaced 2 randomly chosen solvent molecules with 2 chloride ions:
~~~
$ tail -n 4 neutralized.gro
~~~
{: .bash}
~~~
13942SOL    HW243823   7.610   7.577   7.744
13943CL      CL43824   4.747   6.352   6.552
13944CL      CL43825   6.451   0.348   7.643
   7.67638   7.67638   7.67638
~~~
{: .output}
By default genion uses Na+ and Cl- ions. Other ions can be chosen by selecting options **-pname [positive ion]** and **-nname [negative ion]**


The topology file **topol.top** has been updated to include 2 chloride ions:
~~~
$ tail -n 2 topol.top
~~~
{: .bash}
~~~
SOL         13787
CL               2
~~~
{: .output}

#### Caveats and limitations of the random ion placement
As genion places ions randomly the generated system will be most likely in the completely dissociated, energetically unfavourable state. The random placement of charges is particularly problematic if the electric charge of the molecule is big (for example DNA) because ions form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time required for the system equilibration and may affect structural stability. A much better approach would be to place ions according to the electrostatic potential of the molecule. Such method is implemented in the **tleap** program from AMBERTOOLS package.
