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

In this lesson we will go through the steps of setting up a fully solvated protein system for simulation with AMBER/NAMD and GROMACS. While there are many commercial programs and interactive graphical interfaces designed to assist with system preparation. While these tools are easy to use and don't require as much learning efforts as command line tools, they are offer only a limited functionality, and most importantly results obtained with WEB/GUI tools are not reproducible and prone to human error. Therefore, we will focus on system preparation using only scriptable command line driven tools. The emphasis of this lesson is to expose you to the various methods that can be used to create a reproducible molecular modeling workflow by automating preparation and simulation steps. One of the advantages of such approach is that once a workflow script have been developed it can be easily modified for other systems or conditions (for example if an updated version of pdb file is released, you can prepare a new simulation system with a single click).

## Important Things to Check in a PDB File
Small errors in the input structure may cause MD simulations to become unstable or give unrealistic results. The most common problems in PDB files include:

- missing side-chain atoms
- missing fragments
- clashes between atoms
- multiple chains
- alternate conformations
- non-protein molecules (crystallographic waters, ligands, modified amino acids, etc.)
- di-sulfide bonds
- wrong assignment of the N and O atoms in the amide groups of ASN and GLN, and the N and C atoms in the imidazole ring of HIS

Some problems can be identified and corrected automatically (e.g. missing atoms and some clashes) while other may have multiple solutions (e.g. alternate conformations, several protein chains, non-protein molecules, missing residues) and require your decision.

In this section we will learn how to identify and correct for multiple chains, alternate conformations, non-protein molecules, and disulphide bonds.  

#### Downloading a Protein Structure File from PDB
Let's start with creating a directory for MD simulation setup and downloading a protein structure file into it:
~~~
cd ~/scratch
mkdir -p workshop/pdb/1ERT
cd workshop/pdb/1ERT
wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .bash}

#### Checking a PDB File for Presence of Non-Protein Molecules
PDB files are just text files, they contain a lot useful information such as detaitls of the crystallographic experiment, secondary structure, missing residues ... etc. To setup a MD simulation system we will only need the coordinate section, the ATOM. HETATM and TER records. 

The lines beginning with "ATOM" present the atomic coordinates for standard amino acids and nucleotides. All other chemical compounds use the "HETATM" record type. Both of these record types use a simple fixed-column format described [here](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM). 

Any non-protein molecules present in a PDB file would require special treatment. Let's check if the downloaded file has any non-protein atoms. We could install some molecular visualization program for this task, but standard Linux text searching utility *grep* available everywhere is sufficient. Grep searches for a given patterns in the input files and prints out each line that matches a pattern.

~~~
grep "^HETATM " 1ERT.pdb | wc -l
~~~
{: .bash}
~~~
      46
~~~
{: .output}
The '^' expression matches beginning of line. We used the "grep" command to find all lines beginning with the word "HETATM" and then we sent these lines to the character counting command "wc". The output tells us that the downloaded PDB file contains 46 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules. In general PDB files may contain solvents, ions, lipid molecules, protein co-factors, e.t.c. In some cases, these extra components are essential for protein function and should be included in the simulation, while in other cases they were added to facilitate crystallization and are not important. In this introductory lesson, we will ignore non-polymer compounds.

Let's select only protein atoms from the downloaded PDB file and save them in the new file "protein.pdb". Let's  use the molecular visualization and analysis program [VMD](https://www.ks.uiuc.edu/Research/vmd/) to carry out this task. 

To make a program available we need to load its module. After the module is loaded we can start using the program:
~~~
module load StdEnv/2020 gcc vmd
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
The first line created a new molecule from the file 1ERT.pdb. Then we used the **atomselect** method to select all protein atoms from the top molecule. VMD has a powerful Atom Selection Language described [here](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html). Finally we saved the selection in the file protein.pdb.

> ## Selecting Protein Atoms Using Shell Commands
> The same task can be done with one Linux command. Try using *grep* command instead of *vmd* to select protein atoms (the "ATOM" and the "TER" records) from a PDB file and save them in a new file. Hint: the OR operator in *grep* is "\|". The output from *grep* can be redirected into a file using the output redirection operator ">": [command] > [filename].
>
> > ## Solution
> >~~~
> > grep "^ATOM\|^TER " 1ERT.pdb > protein.pdb
> >~~~
> >{: .bash}
> {: .solution}
{: .challenge}

#### Checking a PDB File for Alternate Conformations.
Some PDB files may contain alternate positions of residues. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs (for example pdb2gmx or pdb4amber) will automatically select the first conformation labeled "A" in the "altLoc" column (column 17). If you want to keep a different conformation, all conformations except the desired one must be removed from a PDB file.

Let's check if the downloaded file has any alternate conformations:
~~~
mol new 1ERT.pdb
set s [atomselect top "altloc A"]
$s get resid
set s [atomselect top "altloc B"]
$s get resid
$s get {resid resname name} 
set s [atomselect top "altloc C"]
$s get resid
quit
~~~
{: .vmd}
The output of the commands tells us that residues 20, 43 and 90 have alternate conformations A and B.

Let's select conformations A for residues 43, 90 and conformation B for resid 20 using VMD:
~~~
mol new 1ERT.pdb
set s [atomselect top "(protein and altloc '') or (altloc B and resid 20) or (altloc A and resid 43 90)"]
$s writepdb protein_20B_43A_90A.pdb
quit
~~~
{: .vmd}

You can also check for conformations using grep: 
~~~
grep "^ATOM" 1ERT.pdb | egrep "^.{16}[A-Z]" 
~~~
{: .bash}

First print lines starting with "ATOM", then send these lines into the second grep command. The second grep matches beginning of line, then matches any single character "." 16 times, then matches a single literal A-Z. The symbol "\|" is called "pipe". In Linux you can send output of one command into another command chaining them together. 

> ## Preparing a Protein Coordinate File for Simulation
> 1. Make PDB file containing only conformations "B" from the file 1ERT.pdb
> 2. Make PDB file containing conformations BASP20, AHIS43, BSER90 from the file 1ERT.pdb
> 3. Retrieve the coordinate file for barnase (PDB code 1BNI) and check if there are any alternate conformations in the file.
> 4. Generate structure file for a single molecule (there are 3 molecules in the file). Hint: use chain identifiers to select a molecule. In *vmd* chain A can be chosen using the selection "chain A". Chain identifiers are found in the column 22 of PDB files.
>
> > ## Solution
> >1.Make PDB file containing only conformations “B” from the file 1ERT.pdb:  
> >~~~
> >mol new 1ERT.pdb
> >set s [atomselect top "(protein and altloc '') or (altloc B)"]
> >$s writepdb protein_B.pdb
> >quit
> >~~~
> >{: .vmd}
> >~~~
> >egrep "^.{16}[ B]" 1ERT.pdb | grep "^ATOM"
> >~~~
> >{: .bash}
>  2.Make PDB file containing conformations B HIS43, A SER90, and B ASP20  from the file 1ERT.pdb. 
> >~~~
> >
> >mol new 1ERT.pdb
> >set s [atomselect top "(protein and altloc '') or (altloc A and resid 90) or (altloc B and resid 20 43)"]
> >$s writepdb protein_20B_43B_90A.pdb
> >quit
> >~~~
> >{: .vmd}
> >~~~
> > grep -v "AASP A  20" 1ERT.pdb | grep -v "AHIS A  43" | grep -v "BSER A  90" | grep "^ATOM " > protein_20B_43B_90A.pdb 
> >~~~
> >{: .bash}
> > 3.Download PDB file 1BNI:
> >~~~
> >mkdir ~/scratch/workshop/pdb/1BNI
> >cd  ~/scratch/workshop/pdb/1BNI
> >wget http://files.rcsb.org/view/1BNI.pdb
> >~~~
> >{: .bash}
> > Check for alternate conformations:
> >~~~
> >egrep "^.{16}[A-Z]" 1BNI.pdb | grep ^ATOM
> >~~~
> >{: .bash}
> > 4.Select chain A using shell commands:
> >~~~
> > egrep "^.{21}A" 1BNI.pdb | grep "^ATOM" > barnase.pdb
> >~~~
> > {: .bash}
>> Load molecule directly into *vmd* and select chain A:
> >~~~
> >mol urlload pdb "http://files.rcsb.org/view/1BNI.pdb"
> >set s [atomselect top "protein and chain A"]
> >$s writepdb barnase.pdb
> >quit
> >~~~
> >  {: .vmd}
> {: .solution}
{: .challenge}

#### Checking a PDB File for Disulfide Bonds.
Disulfide bonds are covalent bonds between the sulfur atoms of two cystein residues. They are very important for stabilization of protein structure.
Disulfide bonds are fairly easy to spot in PDB files with any visualisation program. For example [MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can identify disulfide bonds as well as many other problems in PDB files. *GROMACS pdb2gmx* module can add S-S bonds to the topology automatically based on the distance between sulfur atoms (option *-ss*). For MD simulations with *AMBER/NAMD* cross-linked cysteins must be renamed to "CYX" to distinguish them from normal cysteins.

#### BioBB *check_structure* utility
[check_structure](https://pypi.org/project/biobb-structure-checking/) is a command line utility from [BioBB project](https://github.com/bioexcel/biobb) performing [MDWeb](http://mmb.irbbarcelona.org/MDWeb2) structure checking. It includes some structure manipulation options like selecting models or chains, removing components of the system, completing missing atoms, and some quality checking as residue quirality, amide orientation, or vdw clashes.

Installation
~~~
~/scratch/workshop/scripts/install_check_structure.sh 
~~~
{: .bash}

Usage
~~~
module load StdEnv/2020 python
source ~/env-biobb/bin/activate

check_structure commands
check_structure -h
cd ~/scratch/workshop/pdb/1ERT
check_structure -i 1ERT.pdb checkall
check_structure -i 1ERT.pdb -o output.pdb altloc --select A20:A,A43:B,A90:B 
~~~
{: .bash}

#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue may be quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to prepare a simulation. CHARMM-GUI offers useful features, for example the [Membrane Builder](http://www.charmm-gui.org/?doc=input/membrane.bilayer) and the  [Multicomponent Assembler](http://www.charmm-gui.org/?doc=input/multicomp).
