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

- non-protein molecules (crystallographic waters, ligands, modified amino acids, etc.)
- alternate conformations
- missing side-chain atoms
- missing fragments
- clashes between atoms
- multiple copies of the same protein chains
- di-sulfide bonds
- wrong assignment of the N and O atoms in the amide groups of ASN and GLN, and the N and C atoms in the imidazole ring of HIS

Some problems can be identified and corrected automatically (e.g. missing atoms and some clashes) while other may have multiple solutions (e.g. alternate conformations, several protein chains, non-protein molecules, missing residues) and require your decision.

In this section we will learn how to identify and correct for multiple chains, alternate conformations, non-protein molecules, and disulphide bonds.  

#### Downloading a Protein Structure File from PDB
Let's start with downloading a protein structure file:
~~~
cd ~/scratch/workshop/pdb/1ERT
wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .bash}

#### Checking PDB Files for Presence of Non-Protein Molecules
PDB files are just text files. They contain helpful information such as details of the crystallographic experiment, secondary structure, missing residues, etc. To set up an MD simulation system, we will only need the coordinate section, including ATOM, TER, and HETATM records. 

The lines beginning with "ATOM" present the atomic coordinates for standard amino acids and nucleotides. The terminal atom of a protein chain is given in the "TER" record. "HETATM" record type is used for all other chemical compounds. Both of these record types use a simple fixed-column format described [here](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM). 


Any non-protein molecules present in a PDB file require special treatment. In general, PDB files may contain solvents, ions, lipid molecules, protein co-factors, e.t.c. In some cases, these extra components are essential for protein function and should be included in the simulation. Often, some compounds are added to facilitate crystallization and are not essential. In this introductory lesson, we will ignore such compounds.

 Let's select only protein atoms and save them in the new file "protein.pdb". We will use [VMD](https://www.ks.uiuc.edu/Research/vmd/) to carry out this task. 

Load vmd module and start the program:
~~~
module load StdEnv/2020 gcc vmd
vmd
~~~
{: .bash}

~~~
mol new 1ERT.pdb
set prot [atomselect top "protein"]
$prot writepdb protein.pdb
quit
~~~
{: .vmd}
The first line loaded a new molecule from the file 1ERT.pdb. Then we used the **atomselect** method to select all protein atoms from the top molecule. VMD has a powerful Atom Selection Language described [here](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html). Finally we saved the selection in the file "protein.pdb".

>## Selecting ATOM Records Using Shell Commands
>Standard Linux text searching utility `grep` can find and print all "ATOM" records from a PDB file. This is a good example of using Unix command line, and `grep` is very useful for many other purposes such as finding important things in log files. `Grep` searches for a given pattern in files and prints out each line that matches a pattern. 
>1. Check if a PDB file has "HETATM" records using `grep` command.
>2. Select only protein atoms from the file `1ERT.pdb` and save them in the new file `protein.pdb` using `grep` command to select protein atoms (the "ATOM" and the "TER" records). 
>
>Hint: the `OR` operator in *grep* is `\|`. The output from a command can be redirected into a file using the output redirection operator `>`.
>
>>## Solution
>> 1.
>>
>>~~~
>>grep "^HETATM " 1ERT.pdb | wc -l
>>~~~
>>{: .bash}
>>~~~
>>      46
>>~~~
>>{: .output}
>>The `^` expression matches beginning of line. We used the `grep` command to find all lines beginning with the word "HETATM" and then we sent these lines to the character counting command `wc`. The output tells us that the downloaded PDB file contains 46 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules.
>> 
>> 2.
>>~~~
>> grep "^ATOM\|^TER " 1ERT.pdb > protein.pdb
>>~~~
>>{: .bash}
>{: .solution}
{: .challenge}

#### Checking PDB Files for Alternate Conformations.
Some PDB files may contain alternate positions of residues. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs such as `pdb2gmx` or `pdb4amber` will automatically select the first conformation labeled "A" in the "altLoc" column (column 17). 

Sometimes you may want to compare simulations starting from different initial conformations. If you want to select a particular conformation, all conformations except the desired one must be removed from a PDB file.

>## Selecting Alternate Conformations with VMD
>1. Check if the file 1ERT.pdb has any alternate conformations. 
>2. Select conformation A for residues 43, 90. Select conformation B for residue 20. Save the selection in the file "protein_20B_43A_90A.pdb". 
>
>> ## Solution
>>1.
>>
>>~~~
>>mol new 1ERT.pdb
>>set s [atomselect top "altloc A"]
>>$s get resid
>>set s [atomselect top "altloc B"]
>>$s get resid
>>$s get {resid resname name} 
>>set s [atomselect top "altloc C"]
>>$s get resid
>>quit
>>~~~
>>{: .vmd}
>>The output of the commands tells us that residues 20, 43 and 90 have alternate conformations A and B.    
>>  
>>2.
>>  
>>~~~
>>mol new 1ERT.pdb
>>set s [atomselect top "(protein and altloc '') or (altloc B and resid 20) or >(altloc A and resid 43 90)"]
>>$s writepdb protein_20B_43A_90A.pdb
>>quit
>>~~~
>>{: .vmd}
>>
>{:.solution}
{:.challenge}

#### Checking PDB Files for cross-linked cysteines.
Disulfide bonds are covalent bonds between the sulfur atoms of two cysteine residues. They are very important for the stabilization of protein structure.
Disulfide bonds are easy to spot in PDB files with any visualization program. For example, [MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can identify disulfide bonds and many other problems in PDB files. GROMACS `pdb2gmx` utility can automatically add S-S bonds to the topology based on the distance between sulfur atoms (option *-ss*).  
For simulation preparation with the AMBER `tleap` program, cross-linked cysteines must be renamed from "CYS" to "CYX" to distinguish them from normal cysteines. Then the bond between them must be made manually using the `bond` command.

#### Check_structure utility from BioExcel building blocks project
[Check_structure](https://pypi.org/project/biobb-structure-checking/) is a command-line utility from [BioBB project](https://github.com/bioexcel/biobb) for exhaustive structure quality checking (residue chirality, amide orientation, vdw clashes, etc.).  Using this utility, you can perform manipulations with structures, such as selecting chains or conformations, removing components, mutating residues, adding missing atoms, adding hydrogens, etc.  

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
cd ~/scratch/workshop/pdb/1ERT
check_structure -i 1ERT.pdb checkall
check_structure -i 1ERT.pdb -o output.pdb altloc --select A20:A,A43:B,A90:B 
~~~
{: .bash}

#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue may be quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to generate input files for simulation with CHARMM force fields. CHARMM-GUI offers useful features, for example the "Membrane Builder" and the "Multicomponent Assembler".
