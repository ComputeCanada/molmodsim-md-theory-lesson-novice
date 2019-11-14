---
title: "System Preparation"
teaching: 10
exercises: 0
questions:
- "How to simlate a large system?"
objectives:
- "Explain why and when periodic boundary conditions are used"
keypoints:
- "Periodic boundary conditions are used to approximate an infinitely large system"
-  "Simulation system should be neutralized by adding counterions"
---
 Create directory for MD simulation system and download PDB structure file of the *E. Coli* RNase HI:
~~~
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
The file contains 92 non-protein atoms. In this case they are just oxygen atoms of the crystal water molecules. Some pdb files may contain solvents, ions, lipid millecules or protein cofactors. For this introductory lesson we select only protein atoms (the lines starting with "ATOM") and save them in the new file:
~~~
$ grep ^ATOM 1GOA.pdb > 1GOA_protein.pdb
~~~
{: .bash}

Several versions of GROMACS are available on CC clusters. To make GROMACS programs available in the login session we need to load one of the modulefiles:
~~~
$ module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
~~~
{: .bash}
When a GROMACS module is loaded the environment variable EBROOTGROMACS pointing to the GROMACS installation directory will be set. Now we can find out what force fields are available:
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
Let's generate the topology for the 1GOA_protein.pdb using AMBER ff99SBildn force field:
~~~
gmx pdb2gmx -f 1GOA_protein.pdb  -ff amber99sb-ildn -water none
~~~
{: .bash}
By default the pdb2gmx program saves three files: topology **topol.top**, gromacs coordinates **conf.gro** and position restraints **posre.itp**. Once the gromacs coordinate file is created we can add a periodic box to it:
~~~
gmx editconf -f conf.gro -o boxed.gro -c -d 1.0 -bt cubic
~~~
{: .bash}
The option '-c' positions solute in the middle of the box. The option -d specifies the distance between the solute and the box.

~~~
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
~~~
{: .bash}
 To move forward we need to create MD run parameters file. An empy file will be sufficient for now. Using the empty file we can generate the binary topology **.tpr** file:
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
We need to replace 2 water molecules with chloride ions to neutralize this total positive charge. The **genion** program can automatically determine the number of ions required to neutralize a system if the option **-neutral** us used:
~~~
$ echo "SOL" | gmx genion -s solvated.tpr -p topol.top -nname CL -neutral -o neutralized.pdb
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

[setup_GROMACS.sh]({{ page.root }}/code/setup_GROMACS.sh)

~~~
#!/bin/bash
module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
mkdir md_system
cd md_system
wget http://files.rcsb.org/view/1GOA.pdb
grep ^ATOM 1GOA.pdb > 1GOA_protein.pdb
gmx pdb2gmx -f 1GOA_protein.pdb  -ff amber99sb-ildn -water none
gmx editconf -f conf.gro -o boxed.gro -c -d 1.0 -bt cubic
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
touch mdrun.mdp
gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
echo "SOL" | gmx genion -s solvated.tpr -p topol.top -nname CL -neutral -o neutralized.gro
~~~
{: .bash}
