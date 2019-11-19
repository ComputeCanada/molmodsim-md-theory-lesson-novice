---
title: "System Preparation"
teaching: 30
exercises: 0
questions:
- "How to clean up PDB structure file?"
- "How to generate molecular topology file?"
- "How to setup periodic box?"
- "How to solvate MD system?"
- "How to add ions?"
objectives:
- "Explain why and when periodic boundary conditions are used"
keypoints:
- "Periodic boundary conditions are used to approximate an infinitely large system"
- "Simulation system must be neutralized by adding counterions"
---

## Downloading and Cleaning Up PDB Files

 Create directory for MD simulation system and download PDB structure file of the *E. Coli* RNase HI:
~~~
$ cd ~/scratch
$ mkdir md_system
$ cd md_system
$ wget http://files.rcsb.org/view/1GOA.pdb
~~~
{: .bash}
~~~
--2019-11-13 09:19:34--  http://files.rcsb.org/view/1GOA.pdb
Resolving files.rcsb.org (files.rcsb.org)... 128.6.244.12
Connecting to files.rcsb.org (files.rcsb.org)|128.6.244.12|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [text/plain]
Saving to: ‘1GOA.pdb’

1GOA.pdb                [ <=>                ] 136.21K   889KB/s    in 0.2s

2019-11-13 09:19:34 (889 KB/s) - ‘1GOA.pdb’ saved [139482]
~~~
{: .output}
Any non-protein molecules possibly present in PDB files would require special treatment, so let's inspect the downloaded file for pesence of HETATM records:
~~~
$ grep ^HETATM 1GOA.pdb | wc -l
~~~
{: .bash}
~~~
      92
~~~
{: .output}
The **grep** command finds all lines beginning with the word "ATOM" and sends them to the word counting command **wc**. The option "-l" instructs **wc** to count only lines. The output of the command indicates that the file "1GOA.pdb" contains 92 non-protein atoms. In this case they are just oxygen atoms of the crystal water molecules. In general PDB files may contain solvents, ions, lipid molecules, protein cofactors, e.t.c. In some cases these extra components are essential for protein function and should be included in simulation, while in another cases they were simply added to facilitate crystal formation and are not important. In this introductory lesson we will limit simulation to standard protein residues saved in the "ATOM" records. Let's select the lines starting with "ATOM" and save them in the new file protein.pdb:
~~~
$ grep ^ATOM 1GOA.pdb > protein.pdb
~~~
{: .bash}

## Generating Molecular Tolopogy for GROMACS
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
gmx pdb2gmx -f protein.pdb  -ff amber99sb-ildn -water spce
~~~
{: .bash}
By default the pdb2gmx program saves topology, gromacs formatted coordinates, and position restraints in the files **topol.top**, **conf.gro**, and  **posre.itp**, respectively. The names of the output files can be changed by using output options.

## Setting Up Periodic Box and Solvating a System
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

#### How Do I Add Ions and How Many Ions Do I Need to Add?
The GROMACS **genion** command adds ions to a simulation system. Several approaches can be used to determine the number of ions for your simulation:

1. To generate a system representing salt-free solution add enough counterions to neutralize a system (option **-neutral**)

2. To achieve the desired salt concentration add *N* cation/anion pairs and then add counterions to neutralize the system (options **-conc** [Mol/L] **-neutral**). The number of ion pairs is given by the formula: $$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$$. The problem with this approach is that it does not take into account the charge of solute molecules. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with salt concentration that is too high.

3. For more accurate salt concentration you can calculate the number of ions corrected for screening effects in biomolecule using [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) webserver and then pass these numbers to **genion** (options **-nn** and **-np**).

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
As genion places ions randomly the generated system will be most likely in the completely dissociated, energetically unfavorable state. The random placement of charges is particularly problematic if the electric charge of the molecule is big (for example DNA) because ions form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time requred for the system equilibration and may affect structural stability. Much better approach is to place ions according to the electrostatic potential of the molecule. Such method is implemented in the **tleap** program from AMBERTOOLS package.

### Assigning Protonation States
The protonation states of the 7 titratable aminoacids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local microenvironment of the residue and the pH. Highly polar microenvironment will stabilize the charged form, while in less polar microenvironment the neutral form will predominate. The protonation pattern of proteins is critical for their catalytic function and structural stabilization. Numerous MD simulation studies demonstrated importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)]. In a classic MD simulation protonation states are fixed and therefore must be determined and assigned before simulation. Assigning correct protonation states is crucial for realistic MD simulations because inappropriate charge can have drastic effects and invalidate all of the results.

Protonation states in GROMACS can be assigned by the pdb2gmx program. By default pdb2gmx will select charged forms of LYS, ASP or GLU. For HIS it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. Non-standard protonation states can be selected interactively by using options  -lys, -asp, -glu, -his. However, selecting protonation with pdb2gmx can not be automated, and it is cumbersome because pdb2gmx will prompt to select protonation state for each residue.

A better way is to instruct pdb2gmx to use a desired form of aminoacid by changing its name in the input PDB file. The neutral forms of LYS, ASP, and GLU are selected by renaming them to LYN, ASH, and GLH respectively.  The correct form of HIS can be selected by renaming HIS to HIE (proton on NE1), HID (proton on NE2) or HIP (both protons).


#### Limitations of Fixed Protonation State MD Simulations
Molecular dynamics simulations employing constant protonation states have many drawbacks. In real systems conformational changes are often accompanied by changes in protonation pattern. In fixed state molecular dynamics simulations these processes are decoupled hindering understanding proton-coupled conformational dynamics. If proton-coupled dynamics is essentilal for your research consider using constant pH simulations. Constant pH MD is currently implemented in AMBER and NAMD. At the moment it is not officially implemented in GROMACS, but modified GROMACS version is available [[6](https://pubs.acs.org/doi/10.1021/ct200061r)].


### Automating Simulation Setup and Making it Reproducible
The process of molecular dynamics system setup can be automated by saving the  whole sequence of the commands into a text file. This file (shell script) can be subsequently executed to prepare simulation systems from other pdb structure files or it can be used to reproduce a simulation setup. You can download the script carrying out all system preparation steps that we have done manually from here: [setup_GROMACS.sh]({{ page.root }}/code/setup_GROMACS.sh). Make the file setup_GROMACS.sh executable and run it. It will recreate the simulation system in the ~/scratch directory.

### References


#### Notes

Rational for the own tutorial:
- Customize to teach how to use CC for system preparation and MD
- Focus on automation and reproducibility (scripting)
- Discuss alternative ways to do things, shortcomings and sthrengths of different codes.
- Integrate simulation package specific steps into one tutorial.

AMBER:
changeProtState command in ParmEd
#### Prepare protein with HTMD
https://www.acellera.com/molecular-dynamics-playmolecule/
