---
title: "Assigning Protonation States to Residues in a Protein"
teaching: 30
exercises: 10
questions:
- "Why titratable aminoacid residues can have different protonation states?"
- "How to determine protonation state of a residue in a protein?"
- "What are the weaknesses of fixed protonation state simulations?"
objectives:
- "Understand why it is necessary to assign the correct protonation state"
- "Learn how to determine protonation state of a protein"
- "Learn how to assign protonation state to a residue"
keypoints:
- "Assigning correct protonation states of aminoacids in proteins is crucial for realistic MD simulations"
- "Conformational changes in proteins may be accompanied by changes in protonation pattern."
---

The protonation states of titratable aminoacids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local microenvironment of the residue and the pH. For example, highly polar microenvironment will stabilize the charged form, while in less polar microenvironment the neutral form will predominate.

The protonation pattern of proteins is critical for their catalytic function and structural stabilization. In a classic MD simulation protonation states are fixed and therefore must be determined and assigned before simulation. Assigning correct protonation states is crucial for realistic MD simulations. Inappropriate charges can have drastic effects and invalidate all of the results. Numerous MD simulation studies demonstrated importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)].


## How to Detemine Protonation States of Residues in a Protein?

Several free web servers and standalone programs are available for prediction of the *pKa* values of protein residues:

[*H++*](http://biophysics.cs.vt.edu/index.php) server. *H++* calculations are based on the classical continuum electrostatic model (solvent molecules are replaced by a continuous medium with the average properties of the solvent). *H++* uses *AmberTools* modules to preprocess a PDB file and it is capable to generate basic topology and coordinate files for MD simulations in AMBER format. *H++* is a single-conformation method and without intervention from the user it selects, the "A" conformation. For detais of the methodology see [[ref]](https://doi.org/10.1093/nar/gks375).

[*PDB2PQR*](http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/) server. *PDB2PQR* calculations are based on the empirical *pKa* predicting program [*PROPKA3.0*](https://doi.org/10.1021/ct100578z).

[*MCCE*](https://sites.google.com/site/mccewiki/) program. For more rigorous calculations try *MCCE* program. *MCCE* uses the same classical continuum electrostatic approach as H++. In addition *MCCE* calculations account for protein conformational degrees of freedom giving a more accurate picture of coupled ionization and position changes. Taking into account conformational flexibility significantly improves *pKa* prediction. For detais of the methodology see [[ref]](https://doi.org/10.1002/jcc.21222).

[*PKAD*](http://compbio.clemson.edu/pkad) database. *PKAD* is the database of experimentally determined pKa values. Currently it includes *pKa* values for residues in 157 wild-type and 45 mutant proteins [[ref]](https://doi.org/10.1093/database/baz024).

None of the *pKa* prediction methods are perfect. While the average estimated *pKa* values are reasonably accurate, in some extreme cases deviation from the experimental values may be significant. It is advisable to compare results obtained with different methods and if available use the experimentally measured values.

> ## Calculating *pKa*'s
>1. Calculate *pKa*'s of residues in the PDB entry 1RGG using *H++* server.
>2. What protonation states of Asp79 and His53 are appropriate for simulation at *pH* 6?
>3. Repeat calculations using *PDB2PQR* server and compare the results.
>4. Compare calculated *pKa*'s with the experimental. How accurate are the predicted *pKa* values?
>
{: .challenge}
Solution.
if pKa > pH the probability that the residue is protonated is > 50% and we use the protonated form.
if pKa < pH the probability that the residue is protonated is < 50% and we use the deprotonated form.

ASP79 has pKa 7.2, it is protonated at pH 6 and we rename it to ASH
HIS53 has pKa 8.3, it is also protonated at pH 6 and we rename it to HIP

## How to select protonation state of a residue?

### Protonating residues with the *GROMACS pdb2gmx* module
Protonation states can be assigned by the *GROMACS pdb2gmx* program. By default, *pdb2gmx* will select charged forms of LYS, ASP or GLU. For HIS it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. Protonation states can be selected interactively by using options  -lys, -asp, -glu, -his.

The downside of this method is that it can not be scripted. The manual selection is cumbersome because *pdb2gmx* will be prompting to select protonation state for each of the residues in a PDB file. Residue names used by *pdb2gmx* in the output topology file are inconsistent with the output structure files (they are left unchanged in the .pdb and .gro files). This is problematic if you will want to rebuild topology file from PDB structure file or do some further processing of a structure file.

A more consistent and convenient way to select a desired form of aminoacid is to change its name in structure file before generating topology. The neutral forms of LYS, ASP, and GLU can be chosen by renaming them to LYN, ASH, and GLH respectively.  The appropriate form of HIS can be selected by renaming HIS to HIE (proton on NE1), HID (proton on NE2) or HIP (both protons).

Let's change ASP20 and ASP26 in the protein.pdb structure created on the previous step to the neutral form ASH using the *leap* module from the *AmberTools* or *VMD*.

### Selecting protonation states with the *AmberTools leap* module.
~~~
$ module load gcc/5.4.0 openmpi/2.1.1 amber/18
$ tleap -f leaprc.protein.ff14SB
> s = loadpdb protein.pdb
> set {s.20 s.26} name "ASH"
> savepdb s protonated.pdb
> quit
~~~
{: .bash}

### Selecting protonation states with *VMD*.
~~~
$ module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3
$ vmd
vmd> mol new protein.pdb
vmd> set s [atomselect top "resid 20 26"]
vmd> $s set resname ASH
vmd> set s [atomselect top all]
vmd> $s writepdb protonated.pdb
vmd> quit
~~~
{: .bash}


### Limitations of Fixed Protonation State Simulations
Molecular dynamics simulations employing constant protonation states have many drawbacks. In real systems, conformational changes are often accompanied by changes in protonation pattern. In fixed state molecular dynamics simulations, these processes are decoupled hindering understanding proton-coupled conformational dynamics. If proton-coupled dynamics is essential for your research consider using constant *pH* simulations. Constant pH MD is currently implemented in *AMBER* and *NAMD*. At the moment it is not officially implemented in *GROMACS*, but the modified *GROMACS* version is available [[6](https://pubs.acs.org/doi/10.1021/ct200061r)].


> ## Combining all structure preparation steps in one *VMD* script
> Combine all previous steps together and create *VMD* script to prepare MD simulation system for the hydrolaze PDB structure 1RGG. The script should perform the following steps:
>
> 1. Select molecule A
> 2. Remove non-protein molecules
> 3. Select location 'B' for residues 5, 54 and location 'A' for all other residues with alternative locations
> 4. Protonate Asp79 and His53
> 5. Rename CYS 7 and 96 into CYX (cross-linked cystein)
> 6. Save the resulting structure as 1RGG_chain_A_prot.pdb
>
>>## Solution
>> Save the following commands in a file,  e.g. prep_1RGG.vmd
>> ~~~
>># Load 1RGG.pdb into a new (top) molecule
>>mol urlload pdb "http://files.rcsb.org/view/1RGG.pdb"
>># Select and save all chain A protein atoms
>>set s [atomselect top "protein and chain A"]
>>$s writepdb 1RGG_chain_A.pdb
>># Delete the top molecule
>>mol delete top
>># Load chain A into a new molecule
>>mol new 1RGG_chain_A.pdb
>># Protonate ASP79
>>set s [atomselect top "resid 79"]
>>$s set resname ASH
>># Protonate HIS53
>>set s [atomselect top "resid 53"]
>>$s set resname HIP
>># Rename cross-linked cysteins
>>set s [atomselect top "resid 7 96"]
>>$s set resname CYX
>># Select the base and the alternate locations
>>set s [atomselect top "(altloc '') or (altloc A and resid 6 13 42 85 91) or (altloc B and resid 5 54)"]
>># Save the selection
>>$s writepdb 1RGG_chain_A_prot.pdb
>>quit
>>~~~
>>{: .bash}
> Execute the script: vmd -e prep_1RGG.vmd
> {: .solution}
{: .challenge}
