---
title: "Checking PDB Files"
teaching: 30
exercises: 0
questions:
- "What problems are commonly found in PDB files?"
- "Why fixing errors in PDB files is essential for a simulation?"
- "How to find and fix problems in PDB files?"
objectives:
- "Explain why it is necessary to check PDB files before simulation"
- "Explain what things to look for in a PDB file"
- "Explain how to fix the common errors in PDB files"
keypoints:
- "Small errors in the input structure may cause MD simulations to became unstable or give unrealistic results"
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

#### Download a Protein Structure File from PDB
Let's create the md_system directory for the MD simulation setup and download 1ERT protein structure file into it.
~~~
cd ~/scratch
mkdir md_system
cd md_system
wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .bash}

#### Checking a PDB File for Presence of Non-Protein Molecules
Any non-protein molecules present in a PDB file would require special treatment. For this introductory lesson we will select protein only.

Let's check if the downloaded file has any non-protein atoms:
~~~
grep "^HETATM " 1ERT.pdb | wc -l
~~~
{: .bash}
~~~
      92
~~~
{: .output}
We used the **grep** command to find all lines beginning with the word "ATOM" and then we sent these lines to the character counting command **wc**.

The output of the command indicates that the file "1GOA.pdb" contains 92 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules. In general PDB files may contain solvents, ions, lipid molecules, protein cofactors, e.t.c. In some cases, these extra components are essential for protein function and should be included in the simulation, while in other cases they were added to facilitate crystallization and are not important. In this introductory lesson, we will limit the simulation to standard protein residues.

Let's select only protein atoms from the PDB file (the "ATOM" and the "TER" records) and save them in the new file **protein.pdb**:
~~~
grep "^ATOM\|^TER " 1ERT.pdb > protein.pdb
~~~
{: .bash}

We selected only protein atoms from the PDB file using the standard unix shell commands. We coud also use the molecular visualization and analysis program [VMD](https://www.ks.uiuc.edu/Research/vmd/) to carry out this task:

**VMD**
~~~
module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3
vmd
vmd> mol new 1ERT.pdb
vmd> set s [atomselect top "protein"]
vmd> $s writepdb protein.pdb
vmd> quit
~~~
{: .source}

#### Checking a PDB File for Alternate Conformations.

Some PDB files may contain alternate positions of residues. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs (for example pdb2gmx or pdb4amber) will automatically select the first conformation labeled by "A" in the altLoc column. If you want to keep a different conformation, all conformations except the desired one must be removed from the PDB file.

Let's check if the downloaded file has any alternate conformations:
~~~
grep "^.\{16\}[A-Z]" 1ERT.pdb
~~~
{: .bash}

The output of the command indicates that the residues  20, 43 and 90 have alternate conformations A and B.

Let's select the conformation B for all aminoacids:
~~~
grep "^.\{16\}[B]" 1ERT.pdb
~~~
{: .bash}

In addition to the conformation B we need to select the base conformation (space in the altLoc column).
~~~
grep "^.\{16\}[ B]" 1ERT.pdb
~~~
{: .bash}

Finally we can filter out everything except the protein atoms and save the result in the file **protein.pdb**:
~~~
grep "^.\{16\}[ B]" 1ERT.pdb | grep "^ATOM \|^TER " > protein.pdb
~~~
{: .bash}

We can also be more specific and select a specific alternate conformation for different aminoacids (for example AHIS43,  BASP20, and BSER90):
~~~
grep -v "AHIS A  43" 1ERT.pdb | grep -v "ASER A  90" | grep -v "BASP A  20" | grep "^ATOM \|^TER " > protein.pdb
~~~
{: .bash}

**Using VMD:**
~~~
module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3
vmd
vmd> mol new 1ERT.pdb
vmd> set s [atomselect top "(protein and altloc '') or (altloc A and resid 43 90) or (altloc B and resid 20)"]
vmd> $s writepdb protein.pdb
vmd> quit
~~~
{: .bash}

> ## Generating Protein Coordinate File
> 1. Retrieve the coordinate file for barnase (PDB code 1BNI) and check if there any alternate conformations in the file.
> 3. Generate structure file for a single molecule (there are 3 molecules in the file). Hint: use chain identifiers to select a molecule. Chain identifiers are found in the column 22 of PDB files.
>
> > ## Solution
> > Download PDB file:
> > >~~~
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
> >vmd> mol urlload pdb "http://files.rcsb.org/view/1BNI.pdb"
> >vmd> set s [atomselect top "chain A"]
> > vmd> $s writepdb barnase.pdb
> >vmd> quit
> >~~~
> >  {: .bash}
> {: .solution}
{: .challenge}

#### Checking a PDB File for Disulfide Bonds.

1DPX.pdb

#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) Web server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue is quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to prepare a simulation for CHARMM program. CHARMM-GUI offers some other useful features, for example the [Membrane Builder](http://www.charmm-gui.org/?doc=input/membrane.bilayer) and the  [Multicomponent Assembler](http://www.charmm-gui.org/?doc=input/multicomp).
