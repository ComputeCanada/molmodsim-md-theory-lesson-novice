---
title: "Running Simulations"
teaching: 30
exercises: 5
questions:
- "How to minimize energy of a system?"
- "How to run molecular dynamics from cold start?"
objectives:
- "?"
keypoints:
- "?"
---

## 1. Running simulations with AMBER
### 1.1 Energy minimization.
Before simulating a system we need to relax it. Any atomic clashes must be resolved, and potential energy minimized to avoid unphysically large forces that can crash a simulation. 

Let's check our model for clashes. 
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/setup/
source ~/env-biobb/bin/activate
check_structure -i inpcrd_noWAT.pdb checkall
~~~
{: .bash}
~~~
...
8 Steric severe clashes detected
 LYS  124.CA  G    877.N2     0.747 A
 ARG  277.NE  U    878.P      0.845 A
 LYS  525.CE  C    888.N4     0.961 A
 THR  556.CA  A3   898.C4'    0.614 A
 PRO  557.C   A    897.N3     0.487 A
 GLN  558.N   A    897.C4     0.751 A
 THR  559.CA  A3   898.N3     0.435 A
 HIP  822.CD2 C    888.OP1    0.786 A
 ...
~~~
{: .output}

As the RNA model was built without protein, it is expected that the added RNA residues may be placed too close to some aminoacid residues. You can inspect severe clashes between the protein and the RNA residues visually to ensure that there are no severe problems that can not be fixed automatically (such as overlapping ring structures). 

There is nothing too serious that may crash simulation, the clashes will be resolved in the process of energy minimization.  As we want to keep our initial simulation structure as close to the experimental structure as possible we first allow energy minimizer to move freely only new added residues, and restrain all original residues. So we need a list of all original atoms to restrain them. 

In the simulation system residues are renumbered. Chain identifiers are not used. There is a single list of all atoms where atoms and residues are numbered sequentially starting form 1 without any gaps. To make a list of restrained atoms we need to convert PDB residue numbers to simulation residue numbers. Residue number mapping between the original PDB file and the simulation is as follows:

Chain       | Original | Shift | Simulation |
------------|----------|-------|------------|
Protein   A | 1-859    |  -    | 1-859      |
RNA chain C | 1-21     | 859   | 860-880    |
RNA chain D | 1-18     | 898   | 881-898    |
MG ions     |    -     |  -    | 899-901    |

>## Selecting atoms in a simulation system
>1. Create AMBERMASK representing all residues in the simulation system that are given in the original pdb entry 6N4O.   
>2. Modify the AMBERMASK that you created to narrow selection to only backbone atoms. 
>
>Use the following definition of backbone atoms/: CA, N, O, P for protein, and  P, C4', O3' for nucleic. Consult the [AMBER mask selection manual](https://amber-md.github.io/pytraj/latest/atom_mask_selection).
>
>>## Solution
>>1. :22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896   
>>2. (:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896)&(@CA,N,O,P,C4',O3')
>{: .solution}
{: .challenge}

The general minimization strategy is first to restrict all solute atoms with the experimental coordinates and relax all atoms that were added. (solvent, ions and missing fragments). This will help to stabilize the native conformation. There are no strict rules defining how many minimization steps are necessary. The choice will depend on the composition of a simulation system. For a big systems with a large amount of missing residues it is safer to carry out several minimization steps gradually releasing restraints. For example, you can first relax only solvent and ions, then lipid bilayer (if this is a membrane protein), then added fragments, then the original protein side-chains. Having more steps may be unnecessary, but it will not cause any problems. 

Let's do a two a two stage minimization. In the first stage we restrain all original atoms. In the second stage we restrain only the original backbone atoms. 

~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_pmemd/1-minimization
~~~
{: .bash}

Simulation programs can do a lot of different things, and every type of calculation has a number of parameters that allow us to control what will be done. To run a minimization we need to make an input file describing exactly what we want to do and how we want to do it:

- instruct a simulation program to minimize energy  
- choose a method of minimization  
- specify the maximum number of cycles of minimization  
- apply restraints to a subset of atoms (optionally)

A simulation program reads simulation parameters from an input file. Simulation parameters in AMBER are called "FLAGS". The Table below lists some important minimization FLAGS.   

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

Example minimization input file *min_1.in*
~~~
Title line. Energy minimization, stage 1.  
&cntrl 
imin=1, ntmin=0, maxcyc=200,    ! Minimization, method, number of cycles 
ntpr=5,                         ! Print energies every ntpr steps
ntr=1,                          ! Use harmonic cartesian restraints   
restraint_wt=100.0,             ! Restraint force kcal/mol
restraintmask="(:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896)",
&end
END
~~~
{: .file-content}

Allocate resources. The workshop resources are limited, do not ask for more than 8 tasks.
~~~
salloc --time=2:0:0 --mem-per-cpu=2000 --ntasks=8
~~~
{: .bash}
Load AMBER module and run minimization
~~~
module load StdEnv/2020 gcc ambertools
source $EBROOTAMBERTOOLS/amber.sh 
mpiexec sander.MPI -O -i min_1.in -p prmtop.parm7 -c inpcrd.rst7 -ref inpcrd.rst7 -r minimized_1.nc -o mdout_1&
~~~
{: .bash}

The option -O means: overwrite the output files if present.  
The output from the minimization goes into the file *mdout*. The total energy of the system is printed in the lines beginning with "EAMBER =". If minimization is successful we expect to see large negative energies.

>## Why minimization fails?
>When you run minimization with the input file *min1.in* as described above the program crashes after 4 cycles. Try to understand why minimization is unstable, and how to fix the problem. 
>>## Solution
>>With ntmin=0  the minimization method is switched from SD to CG only 4 SD cycles. As initial system has very high energy, 4 SD cycles are not sufficient to relax it well enough for CG to work.   
Try to increase the number of CG cycles to 30 (ntmin=1, ncyc=30).   
Examine mdout. How many SD steps are needed for total energy to become negative?     
Try using an alternative minimization method (SD only or LBFGS). What method converges faster?
>{: .solution}
{: .challenge}

In the second round of minimization constrain only backbone atoms of all original residues. The mask for this selection will be
~~~
(:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896)&(@CA,N,O,P,C4',O3')
~~~
Continue minimization from the restart coordinates:
~~~
mpiexec sander.MPI -O -i min_2.in -p prmtop.parm7 -c minimized_1.nc -ref inpcrd.rst7 -r minimized_2.nc -o mdout_2&
~~~
{: .bash}

### 1.2 Heating
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_pmemd/2-heating
~~~
{: .bash}

Try in the allocated interactive shell: 
~~~
mpiexec sander.MPI -O -i heat.in -p prmtop.parm7 -c minimized_2.nc -ref inpcrd.rst7 -r heated.nc -o mdout&
~~~
{: .bash}

Submit to the queue to run simulation with GPU accelerated pmemd.cuda.
~~~
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --gres=gpu:v100:1 --partition=all_gpus
ml StdEnv/2020 gcc/8.4.0 cuda/10.2 openmpi/4.0.3 amber
pmemd.cuda -O -O -i heat.in -p prmtop.parm7 -c minimized_2.nc -ref inpcrd.rst7 -r heated.nc -o heating.log
~~~
{: .bash}

GPU version runs less than 1 min.

 Flag        | Value  | Description
-------------|--------|-----------
dt           | 0.001  | Time step, ps. Default 0.001 
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

#### Plotting energy components.

Extract selected energy components from MD log and save in a table.
~~~
cpptraj << EOF
readdata heating.log
writedata energy.dat heating.log[Etot] heating.log[TEMP] heating.log[PRESS] heating.log[VOLUME] time 0.1
EOF
~~~
{: .bash}

Read table into pandas dataframe and plot
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

### 1.3 Equilibration
#### Constrained equilibration
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_pmemd/3-equilibration
~~~
{: .bash}

1. Turn on restart flag 
2. Shorten BerendsenPressureRelaxationTime to 1000
3. Decrease restraint force to 10 kcal/mol 
4. Run for 2 ns

Flag         | Value  | Description
-------------|--------|-----------
ntx          | 5      | Coordinates and velocities will be read from a restart file
irest        | 1      | Restart simulations

#### Unconstrained equilibration

1. Switch to Landevin dynamics  
2. Run for 2 ns.

Download log file and examine energy.

Ready for production.

### 1.4 Production

~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/4-production
~~~
{: .bash}

Run for 10 ns

~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --gres=gpu:v100:1 --partition=all_gpus
ml StdEnv/2020 gcc/8.4.0 cuda/10.2 openmpi/4.0.3 amber

pmemd.cuda -O -O -i md.in -p prmtop.parm7 -c equilibrated_2.nc -r rest.nc -o md.log
~~~

#### Managing trajectories
You can remove everything not essential for processing for example water and ions. The following command will remove everything except residues from 1 to 901 and save every second frame in the file mdcrd_nowat.nc
~~~
cpptraj<<EOF
parm prmtop.parm7
trajin mdcrd.nc 1 last 2
strip !(:1-901)
trajout mdcrd_nowat.nc 
parmwrite out prmtop-nowat.parm7
run
EOF
~~~
{: .bash}

## 2. Running simulations with NAMD
### 2.1 Energy minimization
There are only two minimization methods in NAMD, conjugate gradient and simple velocity quenching. All input and output related parameters are configured in input files. NAMD takes only one command line argument, the name of the input file.

NAMD reads constraints from a specially prepared pdb file describing constraint force for each atom in a system. Constraint forces can be given in x,y,z, occupancy or beta columns. 

 Parameter        | Value     | Description
------------------|-----------|-----------
 minimization     | on        | Perform conjugate gradient energy minimization
 velocityQuenching| on        | Perform energy minimization using a simple quenching scheme. 
 numsteps         | integer   | The number of minimization steps to be performed
 constraints      | on        | Activate cartesian harmonic restraints   
 conskfile        | path      | PDB file containing force constant values
 conskcol         | X,Y,Z,O,B | Column of the PDB file to use for the position restraint force constant
 consref          | path      | PDB file containing restraint reference positions

~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/1-minimization
~~~
{: .bash}

As we did in the previous section, in the first round of minimization we restrain all original atoms. Let's prepare PDB file with constraint forces in the occupancy field for this run. Such files can be prepared with VMD. 
~~~
module load vmd
vmd
~~~
{: .bash}
~~~
mol new prmtop.parm7
mol addfile inpcrd.rst7
set sel [atomselect top "all"]
$sel set occupancy 0.0
set sel [atomselect top "noh resid 22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858 860 to 868 870 to 877 879 to 885 889 to 896"]
$sel set occupancy 100.0
set sel [atomselect top "all"]
$sel writepdb constrain_all.pdb
quit
~~~
{: .vmd}

Minimization input file *min_1.in*. In addition to input, output, constraints, and general minimization parameters this file describes force field and PME calculations.
~~~
# Minimization, restrained backbone  
# Input
parmfile       prmtop.parm7  
ambercoor      inpcrd.rst7
# Output 
outputname          minimized_1
numsteps 400
# Constraints
constraints on
conskFile constrain_all.pdb
conskcol O
consref constrain_all.pdb
# Integrator
minimization   on
# AMBER FF settings 
amber on
cutoff         9.0 
pairlistdist   11.0 
switching      off 
exclude        scaled1-4 
readexclusions yes 
1-4scaling     0.83333333 
scnb           2.0 
ljcorrection   on
# PME 
PME                 on 
PMEGridSizeX        140  
PMEGridSizeY        140
PMEGridSizeZ        140 
# Periodic cell 
cellBasisVector1   137  0.0 0.0 
cellBasisVector2   0.0 137  0.0 
cellBasisVector3   0.0 0.0  137  
~~~
{: .file-content}

Run 500 steps of energy minimization:
~~~
module load StdEnv/2020 intel/2020.1.217 namd-multicore/2.14
charmrun ++local +p 8 namd2 min_1.in >&log&
~~~
{:.bash}

In the second round of minimization constrain only backbone atoms of all original residues.   
Prepare force constants file:

~~~
mol new prmtop.parm7
mol addfile inpcrd.rst7
set sel [atomselect top "all"]
$sel set occupancy 0.0
set sel [atomselect top "name CA N O P C4' O3' and resid 22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858 860 to 868 870 to 877 879 to 885 889 to 896"]
$sel set occupancy 10.0
set sel [atomselect top "all"]
$sel writepdb constrain_all_backbone_f10.pdb
quit
~~~
{: .vmd}

Run 1000 minimization steps.
~~~
charmrun ++local +p 8 namd2 min_2.in >&log&
~~~
{:.bash}

### 2.2 Heating 
After energy minimization we have the optimized coordinates that are ready for MD simulation.
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/2-heating
~~~
{: .bash}

| Parameter                      | Value     | Description
|--------------------------------|-----------|-----------
|temperature                     | 150       | Generate velocities at 150 K
|tCouple                         | on        | Use Berendsen thermostat 
|tCoupleTemp                     | 300       | Temperatute of the heat bath
|BerendsenPressure               | on        | Use Berendsen barostat
|BerendsenPressureRelaxationTime | 4000      | Use long relaxation time to slow down box rescaling  
|numsteps                        | 20000     | The number of simulation steps

### 2.3. Equilibration
####  Constrained equilibration
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/3-equilibration
~~~
{: .bash}

Read velocities and box from the restart files  
Shorten BerendsenPressureRelaxationTime to 1000   
Run for 2 ns.  

#### Unconstrained equilibration
Switch to Landevin dynamics  
Run for 2 ns.  

## 2.4 Production
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/4-production
~~~
{: .bash}
Run for 10 ns.  

## 3. Transferring equilibrated system between simulation packages.
Simulation packages have different methods and performance. It is useful to be able to transfer a running simulation from one software to another. Imagine that you started your project with GROMACS, but later realized that you need to run a constant pH simulation. You need to switch to AMBER. Want to study conformational transitions? Gaussian accelerated MD is not available in GROMACS. Another reason to move to AMBER/NAMD. 
Want to apply custom forces - move to NAMD.


#### 3.1. Moving simulation from NAMD to AMBER.
NAMD saves coordinates, velocity and periodic box in separate files. In AMBER all information required to restart simulation is in one text (.rst7) or netcdf (.ncrst) file. To prepare AMBER restart file we first convert namd binary files to amber text restart files with VMD:

~~~
cd ~/scratch/Ago-RNA_sim/sim_pmemd/2-production
cp ~/scratch/Ago-RNA_sim/sim_namd/5-equilibration-unconstrained/{equilibration.coor,equilibration.vel,equilibration.xsc} .
module load intel vmd
vmd
~~~
{: .bash}

~~~
# Convert velocities
mol new equilibration.vel type namdbin
set sel [atomselect top all]
$sel writerst7 vel.rst7
# Convert coordinates
mol new equilibration.coor
set sel [atomselect top "all"]
$sel writerst7 coor.rst7
quit
~~~
{: .vmd}

Now we can read these files along with the extended system (.xsc) file and prepare the complete rst7 file. We will do it with ParmEd, the program for editing AMBER parameter and coordinate files.

NAMD controls pressure differently from AMBER, so when we transfer simulation to AMBER pressure will be too high. It will relax on its own, but we can slightly rescale coordinates and periodic box to speed up pressure equilibration.

~~~
module purge
module load StdEnv/2020 gcc ambertools python scipy-stack
source $EBROOTAMBERTOOLS/amber.sh
python
~~~
{: .bash}

~~~
import parmed as pmd
import numpy as np

sc=1.002
xsc = np.loadtxt('equilibration.xsc')
vel_rst7 = pmd.load_file('vel.rst7')
coor_rst7 = pmd.load_file('coor.rst7')
new_rst7 = pmd.amber.Rst7(natom=vel_rst7.natom)
new_rst7.vels = vel_rst7.coordinates*20.455
new_rst7.coordinates = coor_rst7.coordinates*sc
new_rst7.box = [xsc[1]*sc, xsc[5]*sc, xsc[9]*sc, 90, 90, 90]
new_rst7.title = "Created with ParmEd"
new_rst7.write("restart.rst7")
quit()
~~~
{: .python}

We converted NAMD restart files to AMBER restart and we can continue simulation with AMBER. AMBER suite includes several simulation codes: sander, sander.MPI, pmemd, pmemd.MPI, pmemd.cuda. Sander is free version, pmemd is commercial. Sander and pmemd are serial (CPU only) programs; sander.MPI and pmemd.MPI are parallel (CPU only); and pmemd.cuda is GPU version.

Submitting pmemd.cuda on Siku:
~~~
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --gres=gpu:v100:1 --partition=all_gpus
ml StdEnv/2020 gcc/8.4.0 cuda/10.2 openmpi/4.0.3 amber

pmemd.cuda -O -i pmemd_prod.in -o production.log -p ../../prmtop.parm7 -c restart.rst7
~~~
{:.bash}

PMEMD is highly optimized to do all computations in one GPU, and it runs exceptionally fast. It CANNOT be used efficiently on more than one GPU because of the overhead from moving data between GPUs.

References:
[Amber file formats](https://ambermd.org/FileFormats.php#restart)


#### 3.2 Moving simulation from AMBER to GROMACS.
To transfer simulation to GROMACS in addition to converting restart file we need to convert topology.

First convert AMBER topology to GROMACS
~~~
module load StdEnv/2020 gcc ambertools
source $EBROOTAMBERTOOLS/amber.sh
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_gromacs/0-setup
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

Make index files with groups of atoms that we want to restrain, (one for the original residues and another for backbone of the original residues).

To recap, the original residues are
~~~
:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896
~~~

Position restraints are defined within molecule blocks, they must be included within the correct `moleculetype` block after all atoms are defined (after the `atoms` block). The end of the `moleculetype` block  is a good place.

The atom numbers in position restraints `.itp` files must match the atom numbers in a corresponding `moleculetype` block. Our topology file `gromacs.top`  includes several molecules. The protein is the first, `system1` molecule. The nucleic acids are is the second `system2` and the third `system3` molecules.

WARNING: the genrestr program can generate the correct restraints file only for the first molecule. If you need to restrain the second molecule this tool will not work. To generate a valid constraint file you need to shift indexes in the posre.itp file by the number of atoms in the preceding molecule(s). This can be done by generating separate topology files for each of the RNA chains.

This is a long job, so for now we will restrain only protein. Let's prepare position restraint files for `system1`.
~~~
gmx make_ndx -f inpcrd.gro <<EOF
del5-36
del6-7
r22-120|r126-185|r190-246|r251-272|r276-295|r303-819|r838-858
name 6 Orig_prot
6&aCA
6&aN
6&aO
7|8|9
del 7-9
name 7 Orig_prot_backbone
q
EOF
~~~

Check groups:
~~~
gmx make_ndx  -n index.ndx 
~~~
{: .bash}

Generate positional restraints files, one for all original protein atoms, another for the backbone of the original protein residues.
~~~
gmx genrestr -f inpcrd.gro -fc 500.0 -n index.ndx -o orig_prot.itp<<EOF
Orig_prot
EOF
gmx genrestr -f inpcrd.gro -fc 50.0 -n index.ndx -o orig_prot_backbone.itp<<EOF
Orig_prot_backbone
EOF
~~~
{: .bash}

Add definitions of the position restraints to the topology "gromacs.top". Use a text editor of your choice to insert the following lines at the end of the system1 molecule block:
~~~
#ifdef ORIG_PROT_POSRES
#include "orig_prot.itp"
#endif
#ifdef ORIG_PROT_BACKBONE
#include "orig_prot_backbone.itp"
#endif
~~~
{: .text}

Now we can include any of these two files from the minimization input file by defining a corresponding variable
~~~
; Turn on position restraints
define = -DORIG_PROT_POSRES
; Run parameters
integrator              = steep     
nsteps                  = 400     
; Output control
nstxout                 = 0       
nstvout                 = 0     
nstfout                 = 0        
nstenergy               = 5    
nstlog                  = 5    
nstxout-compressed      = 2000   
; Bond parameters
continuation            = no   
constraint_algorithm    = shake    
constraints             = h-bonds  
; Neighbor-searching
cutoff-scheme           = Verlet   
nstlist                 = 10  
rcoulomb                = 0.8    
rvdw                    = 0.8
DispCorr                = Ener ; anaytic VDW correction 
; Electrostatics
coulombtype             = PME    
pme_order               = 4        
fourier-nx              = 140
fourier-ny              = 140
fourier-nz              = 140
~~~
{: .file-content}

Make binary topology 
~~~
gmx grompp -f min.mdp -p gromacs.top -c inpcrd.gro -r inpcrd.gro -o input.tpr<<EOF
q
EOF
~~~



~~~
import parmed as pmd
amber = pmd.load_file('../prmtop.parm7', '../sim_pmemd/2-production/restart.rst7')
amber.save('gromacs.top')
~~~
{: .python}

Then convert velocities and coordinates:

Amber operates in kcal/mol units for energy, amu for masses,
and Angstoms for distances. For convenience when calculating KE from
velocity, the velocities have a conversion factor built in; as a result the Amber unit of time is (1/20.455) ps.
So to convert Amber velocities from internal units to Ang/ps multiply by 20.455. The number itself is derived from sqrt(1 / ((AMU_TO_KG * NA) / (1000 * CAL_TO_J))).

[AMBER constants](https://github.com/Amber-MD/cpptraj/blob/master/src/Constants.h)
~~~
vel_rst7 = pmd.load_file('../sim_pmemd/2-production/vel.rst7')
amber.velocities = vel_rst7.coordinates*20.455
amber.save('restart.gro')
~~~
{: .python}

~~~
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 gromacs
gmx trjconv -f restart.gro -o restart.trr
gmx make_ndx -f restart.gro
gmx grompp -p gromacs.top  -c restart.gro -t restart.trr -f gromacs_production.mdp
~~~
{: .bash}

Running simulation

~~~
#SBATCH --mem-per-cpu=4000 --time=10:0:0 -c16
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 gromacs
gmx mdrun -s input.tpr
~~~
{: .bash}
