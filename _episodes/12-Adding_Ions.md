---
title: "Solvating a System, Adding Ions and Generating Input Files"
teaching: 30
exercises: 5
questions:
- "Why the simulation system should be neutralized?"
- "How to add water and ions to a simulation system?"
- "How to choose size and shape of a solvent box?"
objectives:
- "Understand why it is necessary to neutralize the simulation system."
- "Neutralize a system."
- "Solvate a macromolecule."
- "Add ions to a simulation system to mimic a salt solution."
- "Generate molecular topology for simulation with GROMACS and NAMD."
keypoints:
- "Simulation system must be neutralized by adding counter-ions to obtain the correct electrostatic energy."
- "Ions are added to a simulations system to mimic composition of a local macromolecule environment."
- "Solvent box should be large enough to allow for unrestricted conformational dynamics of a macromolecule."
---
There are two main reasons to add ions to a simulation system:

1. Under periodic boundary conditions and using grid-based methods to compute Coulomb energy a simulation box interacts with the infinite number of its periodic images. If simulation system is charged the electrostatic energy will essentially add to infinity. We need to neutralize the system by adding counter-ions to obtain the correct electrostatic energy during the simulation.

2. The conformations, dynamics and function of biological macromolecules are sensitive to salt concentration and composition of the local environment.

## Neutralizing a system

Fist we will add enough counter-ions to neutralize the system. The neutralized system will represent a salt-free solution. Ions can be added using two approaches:
1. Solvate the system and then replace random solvent molecules with ions.
2. Place ions according to the electrostatic potential of the macromolecule before solvation.

### Caveats and limitations of the random ion placement
Random placement of ions will generate a system in the completely dissociated, energetically unfavorable state. The random placement of ions is problematic if the electric charge of the macromolecule is big (for example DNA) because ions tend to form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time required for the system equilibration and may affect structural stability of a macromolecule. A better approach is to place ions according to the electrostatic potential of the macromolecule. Such method is implemented in the *leap* module of the *AmberTools*. The *addions* command adds ions to simulation cells near the minima of the solute's electrostatic potential field.

Let's neutralize 1RGG protein using the *leap* module. We will add ions prior to solvation so that the potential from un-equilibrated water will not interfere with ion placement:

~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER
module load StdEnv/2020 gcc amber
tleap
~~~
{: .bash}

~~~
source leaprc.water.tip3p
source leaprc.protein.ff14SB
s = loadpdb ../1RGG_chain_A_prot.pdb
charge s
addions s Na+ 0
~~~
{: .leap}

## Adding Ions to Mimic the Macroscopic Salt Concentration
To mimic the macroscopic salt concentration in the environment being studied we will need to add more cation/anion pairs to the simulation system. The number of ion pairs can be estimated using the formula:

$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$

The drawback of this calculation is that it does not take into account the charge of a macromolecule. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high. For more accurate salt concentration you can calculate the number of ions corrected for screening effects using the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server.

As you can see from the equation above, to calculate the number of ions we need to know the number of water molecules in the simulation system. So we continue our *leap* session and solvate the simulation system. In this lesson we will create a simple cubic solvent box. As we discussed in the episode "Periodic Boundary Conditions", a properly solvated simulation system should have at least 10 <span>&#8491;</span> distance between the solute and the edge of the periodic box after equilibration. Standard practice is to tile a pre-equilibrated solvent box across the system and eliminate solvent molecules which clash with the solute.

When water is added to the system in this way, some VDW voids at the macromolecule and the box interfaces are inevitable because packing is not perfect. In the process of equilibration the water molecules will move to fill the voids and minimize the interaction energy. The box will shrink and the distance between the solute and the edge of the periodic box will become smaller. To compensate for this box shrinkage we need to start with a larger box size than the desired. The rule of thumb is that you need to add at least 2.5 <span>&#8491;</span> to the box size.

We will use the *solvateBox* command to create the periodic solvent box around the macromolecule. The *solvateBox* command has many options. Let's create a cuboid water box around the 1RGG protein. We will use the pre-equilibrated box of SPCE water (SPCBOX), set the minimum distance between the solute and the edge of the box to 15 <span>&#8491;</span>, and request an isometric (cubic) box:
~~~
solvatebox s SPCBOX 15 iso
~~~
{: .leap}

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

> ## Preparing an Aqueous Salt Solution
> How many Na+ and Cl- ions do we need to add to the simulation box with 1RGG protein and 10202 water molecules to prepare 0.15 M salt solution?
> Calculate the number of ions using two methods: the formula above and the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server. For this calculation you need to know molecular weight of the protein. You can calculate it [here](https://www.bioinformatics.org/sms/prot_mw.html). 
>
>>## Solution
>> 1. N_ions = 0.0187 x 0.15 x 10202 = 29. We need to add 35 Na+ and 29 Cl- ions
>> 2. *SLTCAP* calculation with the following input: (MW 11 KDa, 10202 water molecules, charge 6, 150 mM salt) yields 30 Na+ and 24 Cl- ions.
>>
> {: .solution}
{: .challenge}

We already have the neutralized and solvated simulation system, and in the exercise above we determined that we need to add 24 ion pairs to prepare 150 mM salt solution. Let's replace 48 randomly selected water molecules with 24 Na+ and 24 Cl- ions:
~~~
addionsrand s Na+ 24 Cl- 24
~~~
{: .leap}


## Generating Molecular Topology for Simulation with *AMBER* or *NAMD*
Setup of our simulation is almost complete. Our protein has cross-linked cysteine residues, so the last thing to do is to make disulfide bond between Cys7 and Cys96:

~~~
bond s.7.SG s.96.SG
~~~
{: .leap}

We can now save the molecular topology (*parm7*) file and *AMBER* coordinates (*rst7*). To build *GROMACS* topology later we will also save the solvated system in PDB format:

~~~
saveamberparm s 1RGG_chain_A.parm7 1RGG_chain_A.rst7
savepdb s 1RGG_chain_A_solvated.pdb
quit
~~~
{: .leap}


### Summary: script for LEaP to prepare the topology and the coordinate files
Save the following commands in a file, e.g. solvate_1RGG.leap
~~~
source leaprc.water.tip3p
source leaprc.protein.ff14SB
s = loadpdb ../1RGG_chain_A_prot.pdb
addions s Na+ 0
solvatebox s SPCBOX 15 iso
addionsrand s Na+ 24 Cl- 24
bond s.7.SG s.96.SG
saveamberparm s 1RGG_chain_A.parm7 1RGG_chain_A.rst7
savepdb s 1RGG_chain_A_solvated.pdb
quit
~~~
{: .leap}
Execute the script: 
~~~
tleap -f solvate_1RGG.leap
~~~
{: .bash}

## Generating Input Files for Simulation with *GROMACS*.
>## What force fields are available in the loaded *GROMACS* module?
>When the *GROMACS* module is loaded the environment variable *EBROOTGROMACS* will be set. This variable is pointing to the GROMACS installation directory. Knowing where the *GROMACS* installation is we can find out what force fields are available:
>~~~
>module load gromacs
>ls -d $EBROOTGROMACS/share/gromacs/top/*.ff | xargs -n1 basename
>~~~
>{: .bash}
>~~~
>amber03.ff		amber99sb-ildn.ff	gromos45a3.ff
>amber94.ff		amberGS.ff		gromos53a5.ff
>amber96.ff		charmm27.ff		gromos53a6.ff
>amber99.ff		gromos43a1.ff		gromos54a7.ff
>amber99sb.ff		gromos43a2.ff		oplsaa.ff
>~~~
>{: .output}
{: .callout}

### Generate *GROMACS* Topology and Coordinate Files from the Solvated System.

~~~
cd ~/scratch/workshop/pdb/1RGG/GROMACS
~~~
{: .bash}

We can generate GROMACS topology from the complete simulation system prepared previously and saved in the file 1RGG_chain_A_solvated.pdb. For *pdb2gmx* to work correctly we need to rename ions from (Na+, Cl-) to (NA, CL), and rename CYX to CYS:

~~~
ATOM   1444 Na+  Na+    97      -5.058 -11.031  -0.206  1.00  0.00
ATOM   1450 Cl-  Cl-   103      19.451  -3.022   8.361  1.00  0.00
~~~
{: .file-content}

We will also assign ions to chain B.  

~~~
mol new ../1RGG_chain_A_solvated.pdb 
set s [atomselect top "resname 'Na+'"]
$s set resname "NA"
$s set name "NA"
$s set chain "B"
set s [atomselect top "resname 'Cl-'"]
$s set resname "CL"
$s set name "CL"
$s set chain "B"
set s [atomselect top "resname CYX"]
$s set resname "CYS"
set s [atomselect top all]
$s writepdb 1RGG_chain_A_solvated_gro.pdb
quit
~~~
{:.vmd}

- Do it using the global substitution function of the stream editor (sed).

~~~
cat ../1RGG_chain_A_solvated.pdb |\
sed s/"Cl-  Cl-  "/" CL  CL  B"/g |\
sed s/"Na+  Na+  "/" NA  NA  B"/g |\
sed s/CYX/CYS/g > 1RGG_chain_A_solvated_gro.pdb
~~~
{: .bash}


Let's make the topology using the *AMBER ff99SBildn* force field and the *tip3* water model:
~~~
gmx pdb2gmx -f 1RGG_chain_A_solvated_gro.pdb -ff amber99sb-ildn -water tip3 -ignh -chainsep id -ss << EOF 
y
EOF
~~~
{: .bash}

Option       |Value| Description                         
--------------|:---|:-----------------------
ignh          | -  | Ignore hydrogens in file|  
chainsep      |id  | Separate chains by chain ID. Since we assigned ions to chain B pbb2gmx will ignore TER records and put them in a separate chain  
ss            | -  | Interactive S-S bridge selection. Detect potential S-S bonds, and ask for confirmation.

The construct
~~~
<< EOF
y
EOF
~~~
{: .bash}

at the end of the command is to automatically confirm 'y' S-S bond. 

By default *pdb2gmx* program saved three output files:

Type                | Filename 
--------------------|---------- 
topology            | topol.top
coordinates         | conf.gro
position restraints | posre.itp
 
The names of the output files can be changed by using output options *-p*, *-o* and *-i*.

### Prepare the System Using *GROMACS* Module *pdb2gmx*.
To demonstrate how to solvate protein and add ions using *pdb2gmx* we can go back to the protein structure file 1RGG_chain_A_prot.pdb saved before solvation and repeat all system preparation steps with this GROMACS utility. Note that in this case the neutralizing ions will be added in randomly selected positions.

First we generate the topology and the coordinate file using the *AMBER ff99SBildn* force field and the *spc/e* water model:
~~~
gmx pdb2gmx -f ../1RGG_chain_A_prot.pdb -ff amber99sb-ildn -water tip3 -ignh -chainsep id -ss << EOF 
y
EOF
~~~
{: .bash}

Once the gromacs coordinate file *conf.gro* is created we add a periodic box to it:
~~~
gmx editconf -f conf.gro -o boxed.gro -c -d 1.5 -bt cubic
~~~
{: .bash}
The option '-c' positions solute in the middle of the box, the option -d specifies the distance (in nm) between the solute and the box.
Add water
~~~
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
~~~
{: .bash}
 Next, we need to create a binary topology file. For this we need a MD run parameters file. An empy file will be sufficient for now. Using the empty file *'mdrun.mdp'* we can generate the binary topology *'solvated.tpr'* file:
~~~
touch mdrun.mdp
gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
~~~
{: .bash}
When the grompp program runs it generates a lot of diagnostic messages and prints out the net charge. We saved the output in the file *'grompp.log'* so that we can find out what is the total charge of the system:
~~~
grep "total charge" grompp.log
~~~
{: .bash}
~~~
System has non-zero total charge: -5.000000
~~~
{: .output}

Finally we can use the *GROMACS genion* command to replace random solvent molecules with ions. We will first add cation/anion pairs to mimic a desired salt concentration and then neutralize the system by adding sodium ions (the options *-conc* [Mol/L] and *-neutral*). By default genion uses Na+ and Cl- ions. Other ions can be chosen by selecting options *-pname* [positive ion] and *-nname* [negative ion]. We also need to select a target group of solvent molecules to be replaced with ions. We will chose the 'SOL' group which is the default name of the solvent group in *GROMACS*:
~~~
$ echo "SOL" | gmx genion -s solvated.tpr -p topol.top -neutral -conc 0.15 -o neutralized.pdb
~~~
{: .bash}

Let's inspect the last section of the updated topology file:
~~~
tail -n 4 topol.top
~~~
{: .bash}
~~~
Protein_chain_A     1
SOL         11482
NA               38
CL               33
~~~
{: .output}

Create a binary topology file with the new topology and run parameters. 
~~~
gmx grompp -f minimization.mdp -p topol.top -c neutralized.pdb -o solvated_neutral.tpr
~~~
{: .bash}


You can see that 38 sodium and 33 chloride ions were added to the system.

>## Why PDB to GROMACS conversion fails?
>Download the structure file 2F4K from PDB and try to generate the molecular topology with *pdb2gmx*:
>~~~
>wget http://files.rcsb.org/view/2F4K.pdb
>gmx pdb2gmx -f 2F4K.pdb -ff amber99sb-ildn -water spce
>~~~
>{: .bash}
>Why this command fails and how to fix it?
>> ## Solution
>> The file contains two noncanonical norleucine aminoacids (NLE65 and 70). GROMACS does not have parameters for this aminoacid. You have two options: replace norleucine with leucine or add norleucine parameters. 
>>
> {: .solution}
{: .challenge}

### Automation and Reproducibility of a Simulation Setup
The process of molecular dynamics system setup can be automated by saving the  whole sequence of the commands into a text file. This file (shell script) can be easily modified to prepare simulation systems from other PDB structure files or used to reproduce your simulation setup. Using the script, all system preparation steps can be accomplished in seconds.

The shell script performing all preparation steps can be downloaded from [here]({{ page.root }}/code/run_setup.sh). This script will download the molecular structure file from PDB and generate input files for simulation with AMBER, NAMD and GROMACS.


## Running simulations with AMBER
### Energy minimization.
Before simulating a system we need to relax it. Any atomic clashes must be resolved, and potential energy minimized to avoid unphysically large forces that can crash a simulation. 

The general minimization strategy is first to restrict all solute atoms with the experimental coordinates and relax all atoms that were added. (solvent, ions and missing fragments). This will help to stabilize the native conformation. There are no strict rules defining how many minimization steps are necessary. The choice will depend on the composition of a simulation system. For a big systems with a large amount of missing residues it is safer to carry out several minimization steps gradually releasing restraints. For example, you can first relax only solvent and ions, then lipid bilayer (if this is a membrane protein), then added fragments, then the original protein side-chains. Having more steps may be unnecessary, but it will not cause any problems. 

We could do a two a two stage minimization. In the first stage we restrain all original atoms. In the second stage we restrain only the original backbone atoms. Our example protein is very small and we have limited time, so we skip the first step and restrain only protein backbone.

~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER/1_minimization
~~~
{: .bash}

MD programs can do a lot of different things, and every type of calculation has a number of parameters that allow us to control what will be done. To run a minimization we need to make an input file describing exactly what we want to do and how we want to do it:

- instruct a simulation program to minimize energy  
- choose a method of minimization  
- specify the maximum number of cycles of minimization  
- apply restraints to a subset of atoms (optionally)

AMBER MD programs read simulation parameters from an input file. Simulation parameters in AMBER are called "FLAGS". The Table below lists some important minimization FLAGS.   

#### AMBER minimization parameters

| Flag        | Value     | Description
|-------------|-----------|-----------
|imin         |     1     | Turn on minimization
|ntmin        | 0,1,2,3   | Flag for the method of minimization    
|maxcyc       | integer   | The maximum number of cycles of minimization  
|ncyc         | integer   | If NTMIN=1 switch from SD to CG after NCYC cycles  
|ntpr         | integer n | Print energies every n steps    
|ntr          |    1      | Use cartesian harmonic restraints   
|restraint_wt | float     | Restraint force kcal/mol   
|restraintmask| ambermask | Specifies restrained atoms   

#### Methods of minimization

|--|
|0 |Steepest descent+conjugate gradient. The first 4 cycles are steepest descent at the start of the run and after every non-bonded pair-list update.
|1 | For NCYC cycles the steepest descent method is used then conjugate gradient is switched on.
|2 | Steepest descent only
|3 | XMIN family methods. The default is LBFGS (Limited-memory Broyden-Fletcher-Goldfarb-Shanno). It is a popular algorithm in machine learning. The method incrementally learns from previous steps, so that it can make the next step more accurate. It converges considerably faster than CG, but requires more memory.

Minimization input file *min.in*:
~~~
Energy minimization 
&cntrl 
imin=1, ntmin=0, maxcyc=1000,   ! Minimization, method, number of cycles 
ntpr=5,                         ! Print energies every ntpr steps
ntr=1,                          ! Use harmonic cartesian restraints   
restraint_wt=10.0,              ! Restraint force kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

- There are several molecular dynamics programs in AMBER package: *sander, sander.MPI, pmemd, pmemd.MPI, pmemd.cuda*, and *pmemd.cuda.MPI*.  
- *Sander* is distributed free, *pmemd* is high performance binary which requires a license.

Submission script *submit.sh*:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=1000 --time=3:0:0 --ntasks=4
ml --force purge
ml arch/avx2 StdEnv/2020 gcc amber/20  
mpiexec pmemd.MPI -O -i min.in -p prmtop -c inpcrd -ref inpcrd -r minimized.nc -o mdout
~~~
{: .file-content}

This job runs about 30 sec on 4 CPUs.

- The option -O means: overwrite the output files if present.  
- The output from the minimization goes into the file *mdout*. The total energy of the system is printed in the lines beginning with "EAMBER =". 

If minimization is successful we expect to see large negative energies.

### Heating
~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER/2_heating
~~~
{: .bash}

#### Molecular dynamics parameters

Flag        | Value  | Description
-------------|--------|-----------
dt           | 0.001  | Time step, ps. Default 0.001 
ntf          | 2      | Force evaluation. Omit interactions involving H-atoms. Default 1 (complete interaction)
ntc          | 2      | Flag for SHAKE. Bonds involving hydrogen are constrained.
ntt          | 1      | Constant temperature, using the Berendsen weak-coupling algorithm.
tempi        | 150    | Initial temperature. The velocities are assigned from a Maxwellian distribution at TEMPI 
temp0        | 300    | Reference temperature at which the system is to be kept
tautp        | 1      | Time constant, in ps, for heat bath coupling, default is 1 ps. 
ntp          | 1      | Flag for constant pressure dynamics. 1 - MD with isotropic position scaling
barostat     | 1      | Berendsen (default)
pres0        | 1      | Reference pressure, default 1 bar
taup         | 4      | Pressure relaxation time (in ps), default 1 
ntwx         | 1000   | Every ntwx steps, the coordinates will be written to the mdcrd file
ntpr         | 100    | Print energies in the log every 100 steps, default 50 

In the examples bonds with hydrogens are not constrained and the default timestep of 1 fs in used. To turn on SHAKE use ntf=2 and ntc=2.

#### Heating input file:
~~~
Heating 
&cntrl 
ntt=1,                               ! Use Berendsen thermostat
tempi=150,temp0=300,tautp=1,         ! Initial and reference temperature, time constant
ntp=0,                               ! No barostat
ntpr=100,                            ! Print energies every ntpr steps
ntwx=1000,                           ! Write coordinates every ntws steps
nstlim=10000,                        ! Simulate nstlim steps
ntr=1,                               ! Use harmonic cartesian restraints 
restraint_wt=10,                     ! Restraint force kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

#### Heating submission script:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=4
module --force purge
module load arch/avx2 StdEnv/2020 gcc amber/20

mpiexec pmemd.MPI -O -i heat.in -p prmtop -c minimized.nc -ref inpcrd -r heated.nc -o mdout
~~~
{: .bash}

This job runs about 2 min on 4 CPUs.

### Equilibration
#### Constrained equilibration
~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER/3_equilibration
~~~
{: .bash}

- Turn on restart flag.

Flag         | Value  | Description
-------------|--------|-----------
ntx          | 5      | Coordinates and velocities will be read from a restart file
irest        | 1      | Restart simulations


Input file *equilibrate_1.in*:
~~~
&cntrl 
irest=1,ntx=5,                      ! Restart using coordinates and velocities
ntt=1,temp0=300,tautp=1,            ! Use Berendsen thermostat
ntp=1,barostat=1,pres0=1,taup=1,    ! Use Berendsen barostat
ntpr=100,                           ! Print energies every ntpr steps   
ntwx=1000,                          ! Save coordinates every ntwx steps
nstlim=5000,                        ! Simulate nstlim steps
ntr=1,                              ! Turn on restraints
restraint_wt=10,                    ! Restraint force, kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

Submission script *submit_1.sh*:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=4

module --force purge
module load arch/avx2 StdEnv/2020 gcc  amber/20

mpiexec pmemd.MPI -O -i equilibrate_1.in -p prmtop -c heated.nc -ref inpcrd -r equilibrated_1.nc -o equilibration_1.log
~~~
{: .bash}

This job runs about 3 min on 4 CPUs.

#### Unconstrained equilibration

- Use Langevin thermostat + Berendsen barostat.
- Simulate on GPU

Input file *equilibrate_2.in*:
~~~
Equilibration, Langevin thermostat + Berendsen barostat
&cntrl 
irest=1,ntx=5,                     ! Restart using coordinates and velocities
ntt=3,temp0=300,gamma_ln=1.0,      ! Langevin thermostat, collision frequency
ntp=1,barostat=1,pres0=1,taup=1.0, ! Berendsen barostat
ntpr=100,                          ! Print energies every ntpr steps   
ntwx=1000,                         ! Save coordinates every ntwx steps
nstlim=10000000,                   ! Simulate nstlim steps
&end
END
~~~
{: .file-content}

Submission script *submit_2.sh*
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --gres=gpu:1
module --force purge
module load arch/avx2 StdEnv/2020 gcc cuda amber

pmemd.cuda -O  -i equilibrate_2.in -p prmtop -c equilibrated_1.nc  -r equilibrated_2.nc -o equilibration_2.log
~~~
{:.bash}

### Analyzing simulation logs
#### Extract selected energy components from MD log and save in a table using *cpptraj*.

Use the script *~/bin/extract_energies.sh*:
~~~
extract_energies.sh equilibration_1.log
~~~
{: .bash}

The contents of the script *~/bin/extract_energies.sh*:
~~~
#!/bin/bash
echo "Usage: extract_energies simulation_log_file" 

log=$1
cpptraj << EOF
readdata $log
writedata energy.dat $log[Etot] $log[TEMP] $log[PRESS] $log[VOLUME] time 0.1
EOF
~~~
{: .bash}

- Modify the script to add or remove energy components.
- Plot data table with any plotting program.

#### Plot energy components with python

Get interactive allocation and start vncserver
Connect VNC to the node running vncserver.

On the node load the python module:
~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER/3_equilibration
ml StdEnv/2020 python scipy-stack
python
~~~
{: .bash}

Read table from the file *energies.dat* into pandas dataframe and plot it:

~~~
python ~/bin/plot_energies.py
~~~
{: .bash}

*plot_energies.py*:
~~~
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table('energy.dat', delim_whitespace=True)
df.columns=["Time", "Etot", "Temp", "Press", "Volume"]

df.plot(subplots=True, x="Time", figsize=(6, 8))
plt.legend(loc='best')
plt.show()
~~~
{: .python}

### Managing trajectories
You can remove from a trajectory all components that are not essential for analysis, for example water and ions. The following command will remove everything except residues from 1 to 96 and save every second frame.
~~~
cpptraj<<EOF
parm prmtop
trajin mdcrd.nc 1 last 2
strip !(:1-96)
trajout mdcrd_nowat.nc 
run
EOF
cpptraj<<EOF
parm prmtop
parmstrip !(:1-96)
parmwrite out prmtop_nowat.parm7
run
EOF
~~~
{: .bash}

## Transferring equilibrated system between simulation packages.
Simulation packages have different methods and performance. It is useful to be able to transfer a running simulation from one software to another. Imagine that you started your project with GROMACS, but later realized that you need to run a constant pH simulation. You need to switch to AMBER. Want to study conformational transitions? Gaussian accelerated MD is not available in GROMACS. Another reason to move to AMBER/NAMD. 
Want to apply custom forces - move to NAMD.

#### Moving simulation from AMBER to GROMACS.
To transfer simulation to GROMACS we need to convert topology and restart file .

~~~
cd ~/scratch/workshop/pdb/1RGG/AMBER_to_GROMACS
~~~

First convert AMBER topology to GROMACS
~~~
module --force purge
module load arch/avx2 StdEnv/2020 gcc amber gromacs
python
~~~
{: .bash}
~~~
import parmed as pmd
amber = pmd.load_file("prmtop.parm7", "inpcrd.rst7")
amber.save('topol.top')
amber.save('inpcrd.gro')
~~~
{: .python}

Make index file, the default index groups are OK:
~~~
gmx make_ndx -f inpcrd.gro <<EOF
q
EOF
~~~

Generate positional restraints file for the protein backbone.
~~~
gmx genrestr -f inpcrd.gro -fc 500.0 -n index.ndx -o backbone.itp<<EOF
Backbone
EOF
~~~
{: .bash}

Add definitions of the position restraints to the topology "gromacs.top". Use a text editor of your choice to insert the following lines at the end of the syste1 molecule block:
~~~
#ifdef POSRES
#include "backbone.itp"
#endif


[ moleculetype ]
; Name            nrexcl
Na+          3
~~~
{: .text}

Define position restraints on the input file min.mdp:

~~~
; Turn on position restraints
define = -D_POSRES
~~~
{: .text}

Convert AMBER restart to GROMACS restart.
~~~
import parmed as pmd

amber = pmd.load_file("prmtop.parm7", "rest.rst7")
ncrest=pmd.amber.Rst7("rest.rst7")
amber.velocities=ncrest.vels
amber.save("restart.gro")
~~~
{:.python}

Convert GROMACS restart to portable trajectory, and make binary topology
~~~
gmx trjconv -f restart.gro -o restart.trr
gmx grompp -p topol.top  -c restart.gro -t restart.trr -f gromacs_production.mdp
~~~
{: .bash}


