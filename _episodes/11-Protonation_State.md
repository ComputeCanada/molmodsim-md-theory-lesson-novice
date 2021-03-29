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


One of the important jobs of setting up a simulation system is assigning the protonation states and most likely tautomers of the HIS residues. The protonation states of titratable aminoacids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local microenvironment of the residue and the pH. For example, highly polar microenvironment will stabilize the charged form, while in less polar microenvironment the neutral form will predominate. While TYR, LYS, CYS, and ARG are almost always in their standard ptotonation states at physiological pH, GLU, ASP, and HIS can be in non-standard state. You should decide for each of these residues which state is most likely.
 
The protonation pattern of proteins is critical for their catalytic function and structural stabilization. In a classic MD simulation protonation states are fixed and therefore must be determined and assigned before simulation. Assigning correct protonation states is crucial for realistic MD simulations. Inappropriate charges can have drastic effects and invalidate all of the results. Numerous MD simulation studies demonstrated importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)].


## How to Detemine Protonation States of Residues in a Protein?
Several free web servers and standalone programs are available for prediction of the *pKa* values of protein residues. the underlying methods used by all ofthem fall into two categories: empirical (prokka) and Poisson-Boltzmann solvers. The drawbakc of the empirical method is that it works well for a protein similar to the ones it was trained with. It is not guaranteed to give good results for all proteins. Methods implementing Poisson-Boltzmann solution have a sound physical background. Most of them however are limited to sampling only a single protein conformation, and this is the major source of inaccuracy for this method. The exception is MCCE which samples sidechain rotamers giving a more accurate picture of coupled ionization and position changes. It is also limited becaue it does not sample backbone conformations. The most rigorous method is constant pH simulations where ionization states are sampled and dynamically adjusted in a course of MD simulation. This method is computationally very expensive, but recently a very efficient GPU implementation of the constant pH molecular dynamics has been developed and [implemented in AMBER](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00754).

[*H++*](http://biophysics.cs.vt.edu/index.php) server. *H++* calculations are based on the classical continuum electrostatic model (solvent molecules are replaced by a continuous medium with the average properties of the solvent). *H++* uses *AmberTools* modules to preprocess a PDB file and it is capable to generate basic topology and coordinate files for MD simulations in AMBER format. *H++* is a single-conformation method and without intervention from the user it selects, the "A" conformation. For detais of the methodology see [[ref]](https://doi.org/10.1093/nar/gks375).

[PlayMolecule-ProteinPrepare](https://www.playmolecule.org/proteinPrepare/) Preparation process includes the titration of the protonation states using PROPKA 3.1 and addition of missing atoms and overall optimization of the H-network using PDB2PQR 2.1. The first can be manually overriden by examining the protein in our modern 3D webGL viewer. 

[*PDB2PQR*](http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/) server. *PDB2PQR* solves Poisson-Boltzmann equation using the APBS solver. 

[*PROPKA3.0*](https://github.com/jensengroup/propka-3.0) is the empirical pKa prediction software. 

[*MCCE*](https://sites.google.com/site/mccewiki/) program. For more rigorous calculations try *MCCE* program. *MCCE* uses the same classical continuum electrostatic approach as H++. In addition *MCCE* calculations account for protein conformational degrees of freedom giving a more accurate picture of coupled ionization and position changes. Taking into account conformational flexibility significantly improves *pKa* prediction. For detais of the methodology see [[ref]](https://doi.org/10.1002/jcc.21222).

[*PKAD*](http://compbio.clemson.edu/pkad) database. *PKAD* is the database of experimentally determined pKa values. Currently it includes *pKa* values for residues in 157 wild-type and 45 mutant proteins [[ref]](https://doi.org/10.1093/database/baz024).

None of the *pKa* prediction methods are perfect. While the average estimated *pKa* values are reasonably accurate, in some extreme cases deviation from the experimental values may be significant. It is advisable to compare results obtained with different methods and if available use the experimentally measured values.

> ## Calculating *pKa*'s
>1. Calculate *pKa*'s of residues in the PDB entry 1RGG using *H++* server.
>2. What protonation states of Asp79 and His53 are appropriate for simulation at *pH* 6?
>3. Repeat calculations using *PDB2PQR* server and compare the results.
>4. Compare calculated *pKa*'s with the experimental. How accurate are the predicted *pKa* values?
>
>> ## Solution
>>If pKa > pH the probability that the residue is protonated is > 50%, and we use the protonated form.  
>>If pKa < pH the probability that the residue is protonated is < 50% and we use the deprotonated form.
>>
>>ASP79 has pKa 7.2 (experimental 7.37), it is protonated at pH 6 and we rename it to ASH  
>>HIS53 has pKa 8.3 (experimental 8.27), it is also protonated at pH 6 and we rename it to HIP
> {: .solution}
{: .challenge}

## How to select protonation state of a residue?

### Protonating residues with the *GROMACS pdb2gmx* module
Protonation states can be assigned by the *GROMACS pdb2gmx* program. By default, *pdb2gmx* will select charged forms of LYS, ASP or GLU. For HIS it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. Protonation states can be selected interactively by using options  -lys, -asp, -glu, -his.

The downside of this method is that it can not be scripted. The manual selection is cumbersome because *pdb2gmx* will be prompting to select protonation state for each of the residues in a PDB file. Residue names used by *pdb2gmx* in the output topology file are inconsistent with the output structure files (they are left unchanged in the .pdb and .gro files). This is problematic if you will want to rebuild topology file from PDB structure file or do some further processing of a structure file.

A more consistent and convenient way to select a desired form of aminoacid is to change its name in structure file before generating topology. The neutral forms of LYS, ASP, and GLU can be chosen by renaming them to LYN, ASH, and GLH respectively.  The appropriate form of HIS can be selected by renaming HIS to HIE (proton on NE1), HID (proton on NE2) or HIP (both protons).

Let's change ASP20 and ASP26 in the protein.pdb structure created previously from the file 1ERT.pdb to the neutral form ASH.  We can either the *leap* module from the *AmberTools* or *VMD*.

### Selecting protonation states with the *AmberTools leap* module.

~~~
cd  ~/scratch/workshop/pdb/1RGG/1ERT
module --force purge
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 ambertools
source $EBROOTAMBERTOOLS/amber.sh
tleap -f leaprc.protein.ff14SB
~~~
{: .bash}

~~~
s = loadpdb protein.pdb
set {s.20 s.26} name "ASH"
savepdb s protonated.pdb
quit
~~~
{: .leap}

### Selecting protonation states with *VMD*.

~~~
module --force purge
module load StdEnv/2020 intel vmd
vmd
~~~
{: .bash}

~~~
mol new protein.pdb
set s [atomselect top "resid 20 26"]
$s set resname ASH
set s [atomselect top all]
$s writepdb protonated.pdb
quit
~~~
{: .vmd}


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
>>~~~
>> mkdir ~/scratch/workshop/pdb/1RGG
>> cd ~/scratch/workshop/pdb/1RGG
>>~~~
>>{: .bash}
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
>>{: .vmd}
>> Execute the script
>>~~~
>> vmd -e prep_1RGG.vmd
>>~~~
>>{: .bash}
> {: .solution}
{: .challenge}
