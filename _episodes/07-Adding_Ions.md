---
title: "Solvating a System and Adding Ions"
teaching: 20
exercises: 5
questions:
- "Why the simulation system should be neutralized?"
- "How to add water and ions to a simulation system?"
- "How to choose size and shape of a solvent box?"
objectives:
- "Explain why it is necessary to neutlalize the simulation system"
- "Learn how to neutralize a system"
- "Learn how to make a salt solution with a specific molarity"
- "Learn how to solvate a macromolecule"
keypoints:
- "Simulation system must be neutralized by adding counterions to obtain the correct electrostatic energy"
- "Ions are added to a simulations system to mimic composition of a local macromolecule environment"
- "Solvent box should be large enough to allow for unrestricted conformatioal dynamics of a macromolecule"
---
In this lesson we will learn how to add ions and solvent to a simulation system.

There are two main reasons to add ions to a simulation system:

1. Under periodic boundary conditions and using grid-based methods to compute Coulomb energy the box interacts with the infinite number of its periodic images. If simulation system is charged the electrostatic energy will essentially add to infinity. We need to neutralize the system by adding counterions to obtain the correct electrostatic energy during the simulation.

2. The conformations, dynamics and function of biological macromolecules are sensitive to salt concentration and composition of the local environment.

## Neutralizing a system

Fist we will add enough counterions to neutralize the system. The neutralized system will represent a salt-free solution. Ions can be added using two approaches:
1. Solvate the system and then replace random solvent molecules with ions.
2. Place ions according to the electrostatic potential of the macromolecule before solvation.

### Caveats and limitations of the random ion placement
Random placement of ions will generate a system in the completely dissociated, energetically unfavourable state. The random placement of ions is problematic if the electric charge of the macromolecule is big (for example DNA) because ions tend to form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time required for the system equilibration and may affect structural stability of a macromolecule. A better approach is to place ions according to the electrostatic potential of the macromolecule. Such method is implemented in the *leap* module of the *AmberTools*. The *addions* command adds ions to simulation cells near the minima of the solute's electrostatic potential field.

Let's neutralize 1RGG protein using the *leap* module. We will add ions prior to solvation so that the potential from unequilibrated water will not interfere with ion placement:

~~~
module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.1.1 amber/18
tleap
> ; Load parameters for water and ions
> source leaprc.water.spce
> ; Load parameters for protein
> source leaprc.protein.ff14SB
> s = loadpdb 1RGG_chain_A_prot.pdb
> ; Print total charge
> charge s
> ; Add soduim ions to compensate for the negative charge
> addions s Na+ 0
~~~
{: .bash}


## Adding Ions to Mimic the Macroscopic Salt Concentration
To mimic the macroscopic salt concentration in the environment being studied we will need to add more cation/anion pairs to the simulation system. The number of ion pairs can be estimated using the formula:

$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$

The drawback of this calculation is that it does not take into account the charge of a macromolecule. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high. For more accurate salt concentration you can calculate the number of ions corrected for screening effects using the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server.

As you can see from the equation above, to calculate the number of ions we need to know the number of water molecules in the simulation system. So we continue our *leap* session and solvate the simulation system.

### Solvating a system: what size/shape of a solvent box should I use?

Solvated macromolecules rotate during simulations. Furthermore macromolecules may undergo conformational changes. Often these changes are of major interest and shoud not be restricted in any way. If the molecule is not spherical and the box dimension is not large enough rotation will result in the interaction between copies. This artefactual interaction may influence the motions of the system and affect the outcome of the simulation. To avoid these problems the minimum box dimension should be larger than the largest dimension of the macromolecule plus at least 10 <span>&#8491;</span>.

 A cubic box is the most intuitive and common choice, but it is inefficient due to irrelevant water molecules in the corners. The extra water will make your simulation run slower. Ideally you need a sufficiently large sphere of water surrounding the macromolecule, but that's impossible because spheres can't be packed to fill space. A common alternatives that are closer to spherical are the dodecahedron or the truncated octahedron. These shapes work reasonably well for globular macromolecules, but if the solute is elongated there will be a large amount of the unnecessary water located far away from the solute. In this case you may consider constraining the rotational motion [[ref]](https://aip.scitation.org/doi/10.1063/1.480557) and using a smaller rectangular box. But be aware that the box shape itself may influence conformational dynamics by restricting motions in certain directions [[ref]](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20341). This effect may be significant when the amount of solvent is minimal.

In this lesson we will create a simple cubic solvent box. As we discussed above, a properly solvated simulation system should have at least 10 <span>&#8491;</span> distance between the solute and the edge of the periodic box after equilibration. Standard practice is to tile a pre-equilibrated solvent box across the system and eliminate solvent molecules which clash with the solute.

When water is added to the system in this way, some VDW voids at the macromolecule and the box interfaces are inevitable because packing is not perfect. In the process of equilibration the water molecules will move to fill the voids and minimize the interaction energy. The box will shrink and the distance between the solute and the edge of the periodic box wil become smaller. To compensate for this box shrinkage we need to start with a larger box size than the desired 10 <span>&#8491;</span> solvent padding.

We will use the *solvateBox* command to create the periodic solvent box around the macromolecule. The *solvateBox* command has many options. Let's create a cuboid water box around the 1RGG protein. We will use the pre-equilibrated box of SPCE water (SPCBOX), set the minimum distance between the solute and the edge of the box to 15<span>&#8491;</span>, and request an isometric (cubic) box:
~~~
> solvatebox s SPCBOX 15 iso
~~~
{: .bash}

~~~
  Solute vdw bounding box:              40.514 32.235 37.352
  Total bounding box for atom centers:  70.514 70.514 70.514
      (box expansion for 'iso' is  18.6%)
  Solvent unit box:                     18.774 18.774 18.774
  Volume: 399256.044 A^3
  Total mass 194369.824 amu,  Density 0.808 g/cc
  Added 10202 residues.
 ~~~
 {: .output}

Now that we know the number of water molecules in the simulation system, we can add salt to the desired concentration.

> ## Preparing Aqueous Salt Solution
> How many Na+ and Cl- ions do we need to add to the simulation box with 1RGG protein and 10202 water molecules to prepare 0.15 M salt solution?
> Calculate the number of ions using two methods: the formula above and the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server.
>
>>## Solution
>> 1. N_ions = 0.0187 x 0.15 x 10202 = 29. We need to add 35 Na+ and 29 Cl- ions
>> 2. *SLTCAP* calculation with the following input: (MW 11 KDa, 10202 water molecules, charge 6, 150 mM salt) yields 30 Na+ and 24 Cl- ions. You can calculate MW of a protein in kilodaltons [*here*](https://www.bioinformatics.org/sms/prot_mw.html).
>>
> {: .solution}
{: .challenge}

We already have the neutralized and solvated simulation system, and in the exercise above we determined that we need to add 24 ion pairs to prepare 150 mM salt solution. Let's replace 48 randomly selected water molecules with 24 Na+ and 24 Cl- ions:
~~~
> addionsrand s Na+ 24 Cl- 24
~~~
{: .bash}


## Generating Molecular Topology for Simulation with *AMBER* or *NAMD*
Setup of our simulation is almost complete. Our protein has cross-linked cysteine residues, so the last thing to do is to make disulfide bond between Cys7 and Cys96:

~~~
> bond s.7.SG s.96.SG
~~~
{: .bash}

We can now save the molecular topology (*parm7*) file and *AMBER* coordinates (*rst7*). To build *GROMACS* topology later we will also save the solvated system in PDB format:

~~~
> saveamberparm s 1RGG_chain_A.parm7 1RGG_chain_A.rst7
> savepdb s 1RGG_chain_A_solvated.pdb
~~~
{: .bash}


### Summary: script for LEaP to prepare the topology and the coordinate files
Save the following commands in a file, e.g. solvate_1RGG.leap
~~~
source leaprc.water.spce
source leaprc.protein.ff14SB
s = loadpdb 1RGG_chain_A_prot.pdb
addions s Na+ 0
solvatebox s SPCBOX 15 iso
addionsrand s Na+ 24 Cl- 24
bond s.7.SG s.96.SG
saveamberparm s prmtop inpcrd
savepdb s 1RGG_chain_A_solvated.pdb
quit
~~~
{: .bash}
Execute the script: tleap -f solvate_1RGG.leap


## Generating Molecular Topology for Simulation with *GROMACS*.
### What *GROMACS* force fields are available?
When the *GROMACS* module is loaded the environment variable *EBROOTGROMACS* will be set. This variable is pointing to the GROMACS installation directory. Knowing where the *GROMACS* installation is we can find out what force fields are available:
~~~
$ module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
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

### Generate *GROMACS* topology and coordinate files from the solvated system
We can generate gromacs topology from the complete simulation system prepared previously and saved in the file 1RGG_chain_A_solvated.pdb. For *pdb2gmx* to work correctly we need to rename ions to NA and CL and
rename CYX to CYS:

~~~
cat 1RGG_chain_A_solvated.pdb |\
sed s/"Cl-  Cl-  "/" CL  CL  B"/g |\
sed s/"Na+  Na+  "/" NA  NA  B"/g |\
sed s/CYX/CYS/g > 1RGG_chain_A_solvated_gro.pdb
~~~
{: .bash}

Let's make the topology using the *AMBER ff99SBildn* force field and the *spc/e* water model:
~~~
$ module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
$ gmx pdb2gmx -f 1RGG_chain_A_solvated_gro.pdb -ff amber99sb-ildn -water spce -ignh -chainsep id -ss << EOF >log
y
EOF
~~~
{: .bash}

By default *pdb2gmx* program saved topology, *GROMACS* - formatted coordinates, and position restraints in the files *topol.top*, *conf.gro*, and *posre.itp*, respectively. The names of the output files can be changed by using output options *-p*, *-o* and *-i*.

### Solvate protein and add ions using *GROMACS* module *pdb2gmx*.
To demonstate how to solvate protein and add ions using *gromacs* we can go back to the protein structure file 1RGG_chain_A_prot.pdb saved before solvation and repeat all system preparation steps using *pdb2gmx*. Note that in this case the neutralizing ions will be added in randomly selected positions.

First we generate the topology and the coordinate file using the *AMBER ff99SBildn* force field and the *spc/e* water model:
~~~
$ module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
$ gmx pdb2gmx -f 1RGG_chain_A_prot.pdb -ff amber99sb-ildn -water spce -ignh -chainsep id -ss << EOF >log
y
EOF
~~~
{: .bash}

Once the gromacs coordinate file *conf.gro* is created we can add a periodic box to it:
~~~
gmx editconf -f conf.gro -o boxed.gro -c -d 1.5 -bt cubic
~~~
{: .bash}
The option '-c' positions solute in the middle of the box, the option -d specifies the distance (in nm) between the solute and the box.

~~~
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
~~~
{: .bash}
 Next, we need to create MD run parameters file. An empy file will be sufficient for now. Using the empty file mdrun.mdp we can generate the binary topology *.tpr* file:
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
System has non-zero total charge: -5.000000
~~~
{: .output}

The *GROMACS genion* command replaces randomly chosen solvent molecule with ions. We will first add cation/anion pairs to a desired salt concentration and then neutralize the system by adding sodium ions (the options *-conc* [Mol/L] and *-neutral*):
~~~
$ echo "SOL" | gmx genion -s solvated.tpr -p topol.top -neutral -conc 0.15 -neutral -o neutralized.pdb
~~~
{: .bash}

By default genion uses Na+ and Cl- ions. Other ions can be chosen by selecting options *-pname* [positive ion] and *-nname* [negative ion]

The topology file *topol.top* has been updated to include 38 sodium and 33 chloride ions:
~~~
$ tail -n 4 topol.top
~~~
{: .bash}
~~~
Protein_chain_A     1
SOL         11482
NA               38
CL               33
~~~
{: .output}
