---
title: "System Preparation"
teaching: 30
exercises: 0
questions:
- "How to clean up PDB structure file?"
- "How to assign protonation states"
- "How to generate molecular topology file?"
- "How to setup periodic box?"
- "How to solvate MD system?"
- "Why the simulation system should be neutralized"
- "How to add ions?"
objectives:
- "Explain why it is necessary to neutlalize the simulation system"
keypoints:
- "Periodic boundary conditions are used to approximate an infinitely large system"
- "Simulation system must be neutralized by adding counterions"
---

## Downloading and Preparing a PDB File for Simulation

This episode shows how to prepare a protein for simulation. While there are many commercial programs and interactive graphical interfaces designed to assist with the system preparation, we will describe system preparation using free tools that are scriptable, i.e. allow a complete automation of the process.

The process of molecular dynamics system setup can be automated by saving the  whole sequence of the commands into a text file. This file (shell script) can be subsequently executed to prepare simulation systems from other pdb structure files or it can be used to reproduce a simulation setup.

 The emphasis of this lesson is to explore the ways to combine all preparation and simulation steps into a reproducible modeling workflow.

#### Create directory for MD simulation system and download an experimental structure from PDB:
~~~
$ cd ~/scratch
$ mkdir md_system
$ cd md_system
$ wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .bash}

#### Discard Non-Protein Molecules.
Any non-protein molecules present in a PDB file would require special treatment. For this introductory lesson we will select protein only. Let's inspect the downloaded file for pesence of HETATM records:
~~~
$ grep "^HETATM " 1GOA.pdb | wc -l
~~~
{: .bash}
~~~
      92
~~~
{: .output}
The **grep** command finds all lines beginning with the word "ATOM" and sends them to the word counting command **wc**. The option "-l" instructs **wc** to count only lines. The output of the command indicates that the file "1GOA.pdb" contains 92 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules. In general PDB files may contain solvents, ions, lipid molecules, protein cofactors, e.t.c. In some cases, these extra components are essential for protein function and should be included in the simulation, while in other cases they were added to facilitate crystallization and are not important. In this introductory lesson, we will limit simulation to standard protein residues.

In a PDB file protein atoms are saved in lines beginning with the word "ATOM", while all other components are stored in the "HETATM" records. Let's select the lines starting with "ATOM" and save them in a new file. This task can be done using
~~~
$ grep "^ATOM \|^TER " 1GOA.pdb > protein.pdb
~~~
{: .bash}

#### Selecting Specific Aminoacid Alternate Conformations.

Some PDB files may contain several alternate conformations of residues. If an alternate locations are found in a PDB file, both **pdb2gmx** and **pdb4amber** programs will automatically select the first conformation labeled as "A" for each aminoacid. If you want to keep a different conformation all conformations except the desired one must be removed. This task can be done using vmd or grep.

Let's select alternate conformations B for all aminoacids:
~~~
grep "^................[ B]" 1ERT.pdb | grep "^ATOM \|^TER " > protein.pdb
~~~
{: .bash}

We can also be more specific, and select different alternate conformations for different aminoacids (AHIS43 and BASP20):
~~~
grep -v "AHIS A  43" 1ERT.pdb | grep -v "BASP A  20" | grep "^ATOM \|^TER " > protein.pdb
~~~
{: .bash}

> ## Using AMBER/GROMACS/VMD Utilities to Prepare a PDB File
>
> Download PDB file 1ERT and save only protein atoms using pdb4amber utility from the AMBER module
>
>
> > ## Solution
> > AMBER
> >~~~
> > $ module load gcc/5.4.0 openmpi/2.1.1 amber/18
> > $ pdb4amber --pdbid 1ERT --prot -o protein.pdb
> >~~~
> >{: .source}
> > VMD
> >~~~
> > $ module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3
> >$ vmd
> >vmd> mol new 1ERT.pdb
> >vmd> set s [atomselect top "protein and altloc "" or altloc A"]
> >vmd> $s writepdb protein.pdb
> >vmd> quit
> >~~~
> >{: .source}
> {: .solution}
{: .challenge}

## Assigning Protonation States
The protonation states of the 7 titratable aminoacids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local microenvironment of the residue and the pH. Highly polar microenvironment will stabilize the charged form, while in less polar microenvironment the neutral form will predominate. The protonation pattern of proteins is critical for their catalytic function and structural stabilization. Numerous MD simulation studies demonstrated importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)]. In a classic MD simulation protonation states are fixed and therefore must be determined and assigned before simulation. Assigning correct protonation states is crucial for realistic MD simulations because inappropriate charges can have drastic effects and invalidate all of the results.

### Detemining Protonation State of a Protein
- Use [PDB2PQR](http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/). PDB2PQR server can optionally clean up PDB file and optimize hydrogen bonding by sampling rotamers.
- Use local [PROPKA](https://github.com/jensengroup/propka-3.1) installation
- For more accurate multi-conformational calculations try [MCCE](https://sites.google.com/site/mccewiki/)
- Lookup in [PKAD](https://doi.org/10.1093/database/baz024) database of experimental pKa's.

~~~
sed -r -e 's/"|<|~|,/ /g'   WT_pka.csv  | sort -k5 -r -n  | grep ASP
~~~

1ERT ASP26 7.41 computed, 9.9 experimental (PKAD)

### Protonating residues with GROMACS (pdb2gmx)
Protonation states in GROMACS can be assigned by the pdb2gmx program. By default, pdb2gmx will select charged forms of LYS, ASP or GLU. For HIS it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. Non-standard protonation states can be selected interactively by using options  -lys, -asp, -glu, -his. Selecting protonation states with pdb2gmx can not be automated. Manual selection is cumbersome because pdb2gmx will be prompting to select a protonation state for each residue.

A better way to select a desired form of aminoacid is to change its name in the input PDB file. The neutral forms of LYS, ASP, and GLU can be selected by renaming them to LYN, ASH, and GLH respectively.  The correct form of HIS can be selected by renaming HIS to HIE (proton on NE1), HID (proton on NE2) or HIP (both protons).

### Protonating residues with AMBER (tleap)
~~~
$ module load gcc/5.4.0 openmpi/2.1.1 amber/18
$ tleap -f leaprc.protein.ff14SB
> s = loadpdb protein.pdb
> set {s.20 s.26} name "ASH"
> savepdb s protonated.pdb
> quit
~~~
{: .bash}

### Protonating residues with VMD
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
Molecular dynamics simulations employing constant protonation states have many drawbacks. In real systems, conformational changes are often accompanied by changes in protonation pattern. In fixed state molecular dynamics simulations, these processes are decoupled hindering understanding proton-coupled conformational dynamics. If proton-coupled dynamics is essential for your research consider using constant pH simulations. Constant pH MD is currently implemented in AMBER and NAMD. At the moment it is not officially implemented in GROMACS, but the modified GROMACS version is available [[6](https://pubs.acs.org/doi/10.1021/ct200061r)].


## Generating Molecular Topology for GROMACS
Several versions of GROMACS are installed on CC clusters. To start using the progrtam we need to load one of the GROMACS module files, for example:
~~~
$ module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
~~~
{: .bash}
When the module is loaded the environment variable EBROOTGROMACS will be set. This variable is pointing to the GROMACS installation directory. Knowing where the GROMACS installation is we can find out what force fields are available:
~~~
$ ls -d $EBROOTGROMACS/share/gromacs/top/*.ff | xargs -n1 basename | column -c 80
~~~
{: .bash}
~~~
amber03.ff		amber99sb-ildn.ff	gromos45a3.ff
amber94.ff		amberGS.ff		gromos53a5.ff
amber96.ff		charmm27.ff		gromos53a6.ff
amber99.ff		gromos43a1.ff		gromos54a7.ff
amber99sb.ff		gromos43a2.ff		oplsaa.ff
~~~
{: .output}
Let's generate the topology for the 1GOA_protein.pdb using AMBER ff99SBildn force field and spc/e water model:
~~~
gmx pdb2gmx -f protonated.pdb  -ff amber99sb-ildn -water spce
~~~
{: .bash}
By default the pdb2gmx program saves topology, gromacs formatted coordinates, and position restraints in the files **topol.top**, **conf.gro**, and  **posre.itp**, respectively. The names of the output files can be changed by using output options.

## Setting Up Periodic Box and Solvating a Macromolecule
Once the gromacs coordinate file is created we can add a periodic box to it:
~~~
gmx editconf -f conf.gro -o boxed.gro -c -d 1.0 -bt cubic
~~~
{: .bash}
The option '-c' positions solute in the middle of the box, the option -d specifies the distance (in nm) between the solute and the box.

~~~
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
~~~
{: .bash}
 Next, we need to create MD run parameters file. An empy file will be sufficient for now. Using the empty file mdrun.mdp we can generate the binary topology **.tpr** file:
~~~
$ touch mdrun.mdp
$ gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
~~~
{: .bash}
When the grompp program runs it generates a lot of diagnostic messages and prints out the net charge. We saved the output in the file grompp.log so that we can find out what is the total charge of the system:
~~~
grep "total charge" grompp.log
~~~
{: .bash}
~~~
System has non-zero total charge: 2.000000
~~~
{: .output}

## Adding Ions for Accurate Modeling Of Electrostatic Interactions in Solution
Under periodic boundary conditions and using grid-based methods to compute Coulomb energy the box interacts with the infinite number of its periodic images. That is if the system is charged the electrostatic energy will essentially add to infinity. To obtain the correct electrostatic energy during the simulation we need to neutralize the system by adding counterions and then add more ions to represent a desired ionic strength.

#### How Many Ions Do I Need to Add?
The GROMACS **genion** command adds ions to a simulation system. Several approaches can be used to determine the number of ions for your simulation:

1. To generate a system representing salt-free solution add enough counterions to neutralize a system (option **-neutral**)

2. To achieve the desired salt concentration add *N* cation/anion pairs and then add counterions to neutralize the system (options **-conc** [Mol/L] **-neutral**). The number of ion pairs is given by the formula: $$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$$. The problem with this approach is that it does not take into account the charge of solute molecules. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high.

3. For more accurate salt concentration you can calculate the number of ions corrected for screening effects in biomolecule using [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) web server and then pass these numbers to **genion** (options **-nn** and **-np**).

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


### Automating Simulation Setup and Making it Reproducible
The process of molecular dynamics system setup can be automated by saving the  whole sequence of the commands into a text file. This file (shell script) can be subsequently executed to prepare simulation systems from other pdb structure files or it can be used to reproduce a simulation setup. Here is the example shell script carrying out all system preparation steps that we have done manually can be downloaded from here: [setup_GROMACS.sh]({{ page.root }}/code/setup_GROMACS.sh). Make the file setup_GROMACS.sh executable and run it:

 It will download the molecular structure file from PDB and recreate the simulation system in the ~/scratch directory.

### References


#### Notes

Rational to develop our own MD tutorial:
- Can be customized to teach how to use CC clusters for system preparation and MD
- Can be focused on automation and reproducibility by  introducing batch mode and scripting.
- Can discuss alternative ways to do things, shortcomings and sthrengths of different programs/codes.
- Can integrate steps specific for different simulation packages in one tutorial.


AMBER:
changeProtState command in ParmEd
#### Prepare protein with HTMD
https://www.acellera.com/molecular-dynamics-playmolecule/
