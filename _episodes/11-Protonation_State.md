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
One of the critical jobs of setting up a simulation system is assigning the protonation states and most likely tautomers of the HIS residues. The protonation states of titratable amino acids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local micro-environment and the pH. For example, a highly polar micro-environment will stabilize the charged form, while in a less polar micro-environment, the neutral form will predominate. While TYR, LYS, CYS, and ARG are almost always in their standard protonation states at physiological pH, GLU, ASP, and HIS can be in non-standard forms. 
 
The protonation pattern of proteins is critical for their catalytic function and structural stabilization. In a classic MD simulation, protonation states are fixed and must be determined and assigned before simulation. Assigning correct protonation states is crucial for realistic MD simulations. Inappropriate charges can have drastic effects and invalidate the results. Numerous MD simulation studies demonstrated importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)].

## How to Determine Protonation States of Residues in a Protein?
Several web servers and standalone programs are available for the prediction of the *pKa* values of protein residues. The underlying methods used by all of them fall into two categories: empirical (propka) and Poisson-Boltzmann solvers. The drawback of the empirical method is that it works well for proteins similar to those in the training set. It is not guaranteed to give good results for all proteins. Programs implementing the Poisson-Boltzmann continuum electrostatics have a sound physical background. However, most of them are limited to sampling only a single protein conformation, which is the primary source of inaccuracy for this method. The exception is MCCE which samples sidechain rotamers giving a more accurate picture of coupled ionization and position changes. It also has limitations because it does not sample backbone conformations. The most rigorous method is constant pH simulations, where ionization states are dynamically adjusted in MD simulation. This method is computationally costly, but recently a very efficient GPU implementation of the constant pH molecular dynamics has been developed and [implemented in AMBER](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00754).

[*H++*](http://biophysics.cs.vt.edu/index.php) server. *H++* calculations are based on the classical continuum electrostatic model (a continuous medium with the average properties of the solvent replaces solvent molecules ). *H++* uses *AmberTools* modules to pre-process a PDB file, and it is capable of generating basic topology and coordinate files for MD simulations in AMBER format. *H++* is a single-conformation method, and without intervention from the user it selects, the "A" conformation. For details of the methodology, see [[ref]](https://doi.org/10.1093/nar/gks375).

[PlayMolecule-ProteinPrepare](https://www.playmolecule.org/proteinPrepare/) The preparation process includes the determination of the protonation states using PROPKA 3.1, the addition of missing atoms, and overall optimization of the H-bond network using PDB2PQR 2.1. 

[*PDB2PQR*](http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/) server. *PDB2PQR* solves Poisson-Boltzmann equation using the APBS solver. 

[*PROPKA3.0*](https://github.com/jensengroup/propka-3.0) is the empirical pKa prediction software. 

[*MCCE*](https://sites.google.com/site/mccewiki/) program. For more rigorous calculations, try *MCCE* program. *MCCE* uses the same classical continuum electrostatic approach as H++. Besides, *MCCE* calculations consider protein conformational degrees of freedom giving a more accurate picture of coupled ionization and position changes. Taking into account conformational flexibility improves *pKa* prediction significantly. For details of the methodology, see [[ref]](https://doi.org/10.1002/jcc.21222).

[*PKAD*](http://compbio.clemson.edu/pkad) database. *PKAD* is the database of experimentally determined pKa values. Currently it includes *pKa* values for residues in 157 wild-type and 45 mutant proteins [[ref]](https://doi.org/10.1093/database/baz024).

None of the *pKa* prediction methods are perfect. In some extreme cases, deviation from the experimental values may be significant. The best practice is to compare results obtained with different techniques and use the experimentally measured values if available.

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
### Assigning protonation states with the *GROMACS pdb2gmx* module
The *GROMACS pdb2gmx* program can assign protonation states. By default, *pdb2gmx* will select charged forms of LYS, ASP, or GLU. For HIS, it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. You can override the default behavior and select protonation manually using options  -lys, -asp, -glu, -his.

The downside of this method is that it can not be scripted. The manual selection is cumbersome because *pdb2gmx* will be prompting to select the protonation state for each of the residues in a PDB file. Besides, *pdb2gmx* changes residue names only in the output topology file. Residue names are left unchanged in the output .pdb and .gro files.  

A more consistent and convenient way to select the desired form of amino acid is to change its name in the structure file before generating a topology. The neutral states of LYS, ASP, and GLU can be chosen by renaming them  LYN, ASH, and GLH, respectively.  You can select the appropriate form of HIS by renaming HIS to HIE (proton on NE2), HID (proton on ND1), or HIP (both protons).

Let's change ASP20 and ASP26 in the "protein.pdb" file created previously from the "1ERT.pdb" to the neutral form ASH.  We can either the *leap* module from the *AmberTools* or *VMD*.

### Selecting protonation states with the *AmberTools leap* module.

~~~
cd  ~/scratch/workshop/pdb/1ERT
module load StdEnv/2020 gcc/9.3.0 amber
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
module load vmd
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
Molecular dynamics simulations employing constant protonation states have drawbacks. In natural systems, conformational changes are often accompanied by changes in protonation patterns. In fixed state molecular dynamics simulations, these processes are decoupled, hindering understanding of proton-coupled conformational dynamics. If proton-coupled dynamics is essential for your research, consider using constant *pH* simulations. Constant pH MD is currently implemented in *AMBER* and *NAMD*. At the moment, it is not officially implemented in *GROMACS*, but the modified *GROMACS* version is available [[6](https://pubs.acs.org/doi/10.1021/ct200061r)].

References

1. [Constant-pH Molecular Dynamics Simulations for Large Biomolecular Systems](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00875)

2. [ GPU-Accelerated Implementation of Continuous Constant pH Molecular Dynamics in Amber: pKa Predictions with Single-pH Simulations](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00754)

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
>> cd ~/scratch/workshop/pdb/1RGG
>>~~~
>>{: .bash}
>> Save the following commands in a file,  e.g. prep_1RGG.vmd
>> ~~~
>># Load 1RGG.pdb into a new (top) molecule
>>mol pdbload 1RGG
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
