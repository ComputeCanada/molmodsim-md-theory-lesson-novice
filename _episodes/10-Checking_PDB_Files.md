---
title: "Checking PDB Files"
teaching: 20
exercises: 10
questions:
- "What problems are commonly found in PDB files?"
- "Why fixing errors in PDB files is essential for a simulation?"
objectives:
- "Understand why it is necessary to check PDB files before simulation."
- "Learn how to look for problems in PDB files."
- "Learn how to fix the common errors in PDB files."
keypoints:
- "Small errors in the input structure may cause MD simulations to became unstable or give unrealistic results."
---

In this lesson we will go through the steps of setting up a fully solvated protein system for simulation with GROMACS and AMBER/NAMD. While there are many commercial programs and interactive graphical interfaces designed to assist with the system preparation, we will describe system preparation using only free command line driven tools. The emphasis of this lesson is to explore the ways to create a reproducible molecular modeling workflow by automating preparation and simulation steps.

## Important Things to Check in a PDB File
Small errors in the input structure may cause MD simulations to become unstable or give unrealistic results. The most common problems in PDB files include:

- multiple chains
- alternate conformations
- disulfide bonds
- missing amino acids
- missing atoms
- mutations
- non-protein molecules (crystallographic waters, ligands, modified amino acids, etc.)
- clashes between atoms
- wrong assignment of the N and O atoms in the amide groups of ASN and GLN, and the N and C atoms in the imidazole ring of HIS

Some problems can be identified and corrected automatically (e.g. missing atoms) while other may have mutiple solutions (e.g. alternate conformations, several protein chains, non-protein molecules) and require researcher's decision.

#### Downloading a Protein Structure File from PDB
Let's start with creating a directory for MD simulation setup and downloading a protein structure file into it:
~~~
cd ~/scratch
mkdir md_system
cd md_system
wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .bash}

#### Checking a PDB File for Presence of Non-Protein Molecules
Any non-protein molecules present in a PDB file would require special treatment. Let's check if the downloaded file has any non-protein atoms:
~~~
grep "^HETATM " 1ERT.pdb | wc -l
~~~
{: .bash}
~~~
      46
~~~
{: .output}
We used the "grep" command to find all lines beginning with the word "HETATM" and then we sent these lines to the character counting command "wc". The output tells us that the downloaded PDB file contains 46 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules. In general PDB files may contain solvents, ions, lipid molecules, protein cofactors, e.t.c. In some cases, these extra components are essential for protein function and should be included in the simulation, while in other cases they were added to facilitate crystallization and are not important. In this introductory lesson, we will limit the simulation to standard protein residues.

Let's select only protein atoms from the downloaded PDB file and save them in the new file "protein.pdb". We will use he molecular visualization and analysis program [VMD](https://www.ks.uiuc.edu/Research/vmd/) to carry out this task. To make the program available we need to load its module. After the module is loaded we can start using the program:
~~~
module load nixpkgs/16.09 intel/2016.4 vmd/1.9.3
vmd
~~~
{: .bash}

~~~
mol new 1ERT.pdb
set s [atomselect top "protein"]
$s writepdb protein.pdb
quit
~~~
{: .vmd}


> ## Selecting Protein Atoms Using Shell Commands
> Use "grep" command to select protein atoms (the "ATOM" and the "TER" records) from a PDB file and save them in a new file.
>
> > ## Solution
> >~~~
> > grep "^ATOM\|^TER " 1ERT.pdb > protein.pdb
> >~~~
> >{: .bash}
> {: .solution}
{: .challenge}


#### Checking a PDB File for Alternate Conformations.

Some PDB files may contain alternate positions of residues. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs (for example pdb2gmx or pdb4amber) will automatically select the first conformation labeled "A" in the "altLoc" column. If you want to keep a different conformation, all conformations except the desired one must be removed from a PDB file.

Let's check if the downloaded file has any alternate conformations:
~~~
grep "^.\{16\}[A-Z]" 1ERT.pdb
~~~
{: .bash}

The output of the command tells us that residues 20, 43 and 90 have alternate conformations A and B.

Let's select conformations A for residues 43, 90 and conformation B for resid 20 using VMD:
~~~
mol new 1ERT.pdb
set s [atomselect top "(protein and altloc '') or (altloc A and resid 43 90) or (altloc B and resid 20)"]
$s writepdb protein.pdb
quit
~~~
{: .vmd}


> ## Preparing a Protein Coordinate File for Simulation
> 1. Make PDB file containing only conformations "B" from the file 1ERT.pdb
> 2. Make PDB file containing only conformations AHIS43, BASP20, BSER90 from the file 1ERT.pdb
> 3. Retrieve the coordinate file for barnase (PDB code 1BNI) and check if there are any alternate conformations in the file.
> 4. Generate structure file for a single molecule (there are 3 molecules in the file). Hint: use chain identifiers to select a molecule. Chain identifiers are found in the column 22 of PDB files.
>
> > ## Solution
> Make PDB file containing only conformations “B” from the file 1ERT.pdb:
> >~~~
> >grep "^.\{16\}[ B]" 1ERT.pdb | grep "^ATOM \|^TER " > protein.pdb
> >~~~
> >{: .bash}
>  Make PDB file containing conformations AHIS43, BASP20, BSER90 from the file 1ERT.pdb
>>~~~
>> grep -v "AHIS A  43" 1ERT.pdb | grep -v "ASER A  90" | grep -v "BASP A  20" | grep "^ATOM \|^TER " > protein.pdb
>>~~~
>>{: .bash}
> > Download PDB file 1BNI:
> >~~~
> >wget http://files.rcsb.org/view/1BNI.pdb
> >~~~
> >{: .bash}
> > Check for alternate conformations:
> >~~~
> >grep "^.\{16\}[A-Z]" 1BNI.pdb | grep ^ATOM
> >~~~
> >{: .bash}
> > Select chain A using shell commands:
> >~~~
> > grep "^.\{21\}A \| ^TER" 1BNI.pdb | grep "^ATOM\|^TER" > barnase.pdb
> >~~~
> > {: .bash}
>> Select chain A using VMD:
> >~~~
> >mol urlload pdb "http://files.rcsb.org/view/1BNI.pdb"
> >set s [atomselect top "chain A"]
> >$s writepdb barnase.pdb
> >quit
> >~~~
> >  {: .vmd}
> {: .solution}
{: .challenge}

#### Checking a PDB File for Disulfide Bonds.
Disulfide bonds are covalent bonds between the sulfur atoms of two cystein residues. They are very important for stabilization of protein structure.
Disulfide bonds are fairly easy to spot in PDB files with any visualisation program. For example [MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can identify disulfide bonds as well as many other problems in PDB files. *GROMACS pdb2gmx* module can add S-S bonds to the topology automatically based on the distance between sulfur atoms (option *-ss*). For MD simulations with *AMBER/NAMD* cross-linked cysteins must be renamed to "CYX" to distinguish them from normal cysteins.


#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue is quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to prepare a simulation for CHARMM program. CHARMM-GUI offers some other useful features, for example the [Membrane Builder](http://www.charmm-gui.org/?doc=input/membrane.bilayer) and the  [Multicomponent Assembler](http://www.charmm-gui.org/?doc=input/multicomp).
