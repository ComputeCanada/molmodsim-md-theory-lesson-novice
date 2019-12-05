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
So far ee have learned to identify common problems in PDB files, correct them and fix protonation states. Now we can start adding ions and solvent to complete the simulation system setup. There are two main reasons to add ions to a simulation system:

I. Under periodic boundary conditions and using grid-based methods to compute Coulomb energy the box interacts with the infinite number of its periodic images. If simulation system is charged the electrostatic energy will essentially add to infinity. We need to neutralize the system by adding counterions to obtain the correct electrostatic energy during the simulation.

II. The conformations, dynamics and function of biological macromolecules are sensitive to salt concentration and composition of the local environment.


## Neutralizing a system

Fist we should add enough counterions to neutralize the system. The neutralized system would represent a salt-free solution. Ions can be added using two approaches:
1. Solvate the system first, then replace random solvent molecules with ions.
2. Place ions according to the electrostatic potential of the molecule before solvation.


### Caveats and limitations of the random ion placement
Random placement of ions will generate a system in the completely dissociated, energetically unfavourable state. The random placement of charges is particularly problematic if the electric charge of a macromolecule is big (for example DNA) because ions tend to form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time required for the system equilibration and may affect structural stability of a macromolecule. A better approach is be to place ions according to the electrostatic potential of the molecule. Such method is implemented in the **leap** module of AmberTools package.The addions command adds ions to simulation cells near the minima of the solute's electrostatic potential field.

Let's neutralize 1RGG proteins using **leap**. We will do it prior to solvation so that the potential from unequilibrated water does not interfere with ion placement:

~~~
module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.1.1 amber/18
tleap
>source leaprc.water.spce
>source leaprc.protein.ff14SB
>s = loadpdb 1RGG_chain_A_prot.pdb
>addions s Na+ 0
~~~
{: .bash}


## Adding Ions to Mimic the Macroscopic Salt Concentration
To mimic the macroscopic salt concentration in the environment being studied we will need to add more ions to a simulation system. The number of ion pairs can be estimated using the formula:

$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$

The problem with this approach is that it does not take into account the charge of macromolecules. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high. For more accurate salt concentration you can calculate the number of ions corrected for screening effects in biomolecule using [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) web server.

As you can see from the equation above, to calculate the number of ions we need to know the nuber of water molecules in the simulation system. So we continue our leap session and solvate the simulation system:

~~~
> solvatebox s SPCBOX 15 iso
> savemol2 s 1RGG_chain_A_solvated.mol2 1
~~~
{: .bash}

~~~
  Solute vdw bounding box:              40.514 32.235 37.352
  Total bounding box for atom centers:  70.514 70.514 70.514
      (box expansion for 'iso' is  18.6%)
  Solvent unit box:                     18.774 18.774 18.774
  Volume: 399256.044 A^3
  Total mass 194369.824 amu,  Density 0.808 g/cc
  Added 10202 residues.
 ~~~
 {: .output}


Standard practice is to tile a pre-equilibrated solvent box across the system and eliminate solvent molecules which clash with the solute. **solvateBox** command has many options. Here we created a cuboid water box around the neutralized 1RGG protein. We set the minimum distance between any atom in the solute and the edge of the periodic box to 15<span>&#8491;</span>, and we requested an isometric box. The solvated system contains 10202 water molecules.

> ## Preparing Aqueous Salt Solution
> How many Na+ and Cl- ions do we need to add to the simulation box with 1RGG protein and 10202 water molecules to prepare 0.15 M salt solution?
> Calculate the number of ions using two methods: the formula above and [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) websever.
>
>>## Solution
>> 1. Using the formula above: N_ions = 0.0187 x 0.15 x 10202 = 29. We need to add 35 Na+ and 29 Cl- ions
>> 2. SLTCAP calculation using MW=11 KDa, 10202 water molecules, charge = 6 and 150 mM salt concentration yields 30 Na+ and 24 Cl- ions.
>> [Calculate MW of a protein in kilodaltons](https://www.bioinformatics.org/sms/prot_mw.html)
>>
> {: .solution}
{: .challenge}

We already have the neutralized and solvated simulation system, so we will replace 48 randomly selected water molecules with ions to prepare 150 mM salt solution:
~~~
> addionsrand s Na+ 24 Cl- 24
~~~
{: .bash}

Setup of our simulation is completed and we can now save the topology and the coordinate files. .. Wait we forgot something important. Disulfide bond:

~~~
>bond s.7.SG s.96.SG
>saveamberparm s prmtop inpcrd
>quit
~~~
{: .bash}


## To summarize the whole solvation process:
Paste the following commants into solvate_1RGG.leap
~~~
source leaprc.water.spce
source leaprc.protein.ff14SB
s = loadpdb 1RGG_chain_A_prot.pdb
addions s Na+ 0
solvatebox s SPCBOX 15 iso
addionsrand s Na+ 24 Cl- 24
bond s.7.SG s.96.SG
saveamberparm s prmtop inpcrd
quit
~~~

tleap -f solvate_1RGG.leap
