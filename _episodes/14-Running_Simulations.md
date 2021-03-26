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

### 5. Energy minimization.
Before simulationg molecular dynamics we need to relax the system. Any atomic clashes must be resolved, and potential energy minimized to avoid unphysically large forces that can crash simulation. 
Let's check our model for clashes. 
~~~
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

As the RNA model was built without protein, it is expected that the added RNA residues may be placed too close to some aminoacid residues. You can inspect severe clashes between the protein and the RNA residues visually to ensure that there are no severe problems such as overlapping rings that can not be fixed automatically. 

There is nothing too serious that may crash simulation, the clashes will be resolved in the process of energy minimization.  As we want to keep our initial simulation structure as close to the experimental structure as possible we first allow energy minimizer to move freely only new added residues, and restrain all original residues. So we need a list of all original atoms to restrain them. 

In the simulation system resisues are renumbered. All residues in simulation systems are numbered sequentially starting with 1, and chain identifiers do not exist. Thus, the residue number mapping between the original pdb file and the simulation is as follows:

Chain       | Original | Shift | Simulation |
------------|----------|-------|------------|
Protein   A | 1-859    |  -    | 1-859      |
RNA chain C | 1-21     | 859   | 860-880    |
RNA chain D | 1-18     | 898   | 881-898    |
MG ions     |    -     |  -    | 899-901    |

Considering this mapping, the AMBERMASK selecting all original residues is:
~~~
:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896
~~~

#### 5.1 Energy minimization with AMBER

Minimization parameters

| Flag        | Value     | Description
|-------------|-----------|-----------
|imin         |     1     | Turn on minimization
|ntmin        | 0 - 4     | Flag for the method of minimization    
|maxcyc       | integer   | The maximum number of cycles of minimization  
|ncyc         | integer   | If NTMIN=1 switch from SD to CG after NCYC cycles  
|ntpr         | integer n | Print energies every n steps    
|ntr          |    1      | Use cartesian harmonic restraints   
|restraint_wt | float     | Restraint force kcal/mol   
|restraintmask| ambermask | Specifies restrained atoms   

Methods of minimization

|--|
|0 |Steepest descent+conjugate gradient. The first 4 cycles are steepest descent at the start of the run and after every nonbonded pairlist update.
|1 | For NCYC cycles the steepest descent method is used then conjugate gradient is switched on.
|2 | Steepest descent only
|3 | XMIN family methods. The default is LBFGS (Limited-memory Broyden-Fletcher-Goldfarb-Shanno). It is a popular algorithm in machine learning. The method incrementally learns from previous steps, so that it can make the next step more accurate. It converges considerably faster than CG, but requires more memory.

Minimization input file *min1.in*
~~~
Energy minimization
&cntrl
imin=1, ntmin=0, maxcyc=400,
ntpr=5,
ntr=1,
restraint_wt=100,
restraintmask="(:22-120,126-185,190-246,251-272,276-295,303-819,838-858,860-868,870-877,879-885,889-896)",
&end
END
~~~
{: .file-content}

For convenience make links in the working directory pointing to the topology and the initial coordinates.
~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_pmemd/1-minimization
ln -s ../../setup/prmtop.parm7 prmtop.parm7
ln -s ../../setup/inpcrd.rst7 inpcrd.rst7
~~~
{: .bash}

Allocate resources
~~~
salloc --time=2:0:0 --mem-per-cpu=2000 --ntasks=8
~~~
{: .bash}
Load AMBER module and run minimization
~~~
module load StdEnv/2020 gcc ambertools
mpiexec sander.MPI -O -i min1.in -p prmtop.parm7 -c inpcrd.rst7  -ref inpcrd.rst7 -r minimized_1.nc&
~~~
{: .bash}

The option -O means: overwrite the output files if present.  
The output from the minimization goes into the file *mdout*. The total energy of the system is printed in lines beginning with "EAMBER =". If minimization is successful we expect to see large negative energies.

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
mpiexec sander.MPI -O -i min2.in -p prmtop.parm7 -c minimized_1.nc -ref inpcrd.rst7 -r minimized_2.nc&
~~~

#### 5.2 Energy minimization with NAMD

There are only two minimization methods in NAMD, conjugate gradient and simple velocity quenching. All input and output related parameters are configured in input files. NAMD takes only one command line argument, the name of the input file.

NAMD reads constraints from a specially prepared pdb file describing constraint force for each atom in a system. Constraint forces can be given in either occupancy or beta columns. 

| Parameter        | Value     | Description
|------------------|-----------|-----------
| minimization     | on        | Perform conjugate gradient energy minimization
| velocityQuenching| on        | Perform energy minimization using a simple quenching scheme. 
| numsteps         | integer   | The number of minimization steps to be performed
| constraints      | on        | Activate cartesian harmonic restraints   
| conskfile        |  path     | PDB file containing force constant values
| conskcol         | X,Y,Z,O,B | Column of the PDB file to use for the position restraint force constant
| consref          | path      | PDB file containing restraint reference positions

~~~
cd ~/scratch/workshop/pdb/6N4O/simulation/sim_namd/1-minimization
ln -s ../../setup/prmtop.parm7 prmtop.parm7
ln -s ../../setup/inpcrd.rst7 inpcrd.rst7
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
$sel writepdb constrain_all_6n4o_residues.pdb
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
conskFile constrain_all_6n4o_residues.pdb
conskcol O
consref constrain_all_6n4o_residues.pdb
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
PMEGridSizeX        128  
PMEGridSizeY        128
PMEGridSizeZ        128 
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
$sel writepdb constrain_backbone_all_6n4o_residues_f10.pdb
quit
~~~
{: .vmd}

Run 1000 steps of minimization.
~~~
charmrun ++local +p 8 namd2 min_2.in >&log&
~~~
{:.bash}

### 6. Heating and equilibration.

After energy minimization we have the optimized coordinates that are ready to use for MD simulation.

#### 6.1. Heating

Generate the initial velocities at 150 K

~~~
temperature 150
~~~
{: .file-content}

Use Berendsen thermostat and barostat

~~~
tCouple on
tCoupleTemp 300
BerendsenPressure on
~~~
{: .file-content}

Long pressure relaxation time to prevent box from changing too fast

~~~
BerendsenPressureRelaxationTime 4000
~~~
{: .file-content}

Run heating for 20 ps.

#### 6.2 Constrained equilibration

Read velocities and box from restart Files

Shorten Berendsen Pressure Relaxation Time

Run for 2 ns.

Download and examine energy.

#### 6.3 Unconstrained equilibration

Switch to Landevin dynamics

Run for 2 ns.

Download and examine energy.

Ready for production.



Sbatch file for running simulation on a single (GPU) node on Siku:

~~~
#!/bin/bash
#SBATCH -c8 --mem-per-cpu=4000 --time=3:0:0 --gres=gpu:v100:1 --partition=all_gpus
module load StdEnv/2020 intel/2020.1.217 cuda/11.0 namd-multicore

charmrun ++local +p $SLURM_CPUS_PER_TASK namd2 namd_min.in
~~~
{: .file-content}

Sbatch file for running simulation on multiple nodes (CPU):

~~~
#!/bin/bash
#SBATCH --ntasks=8 --mem-per-cpu=4000 --time=3:0:0

module load StdEnv/2020  intel/2020.1.217 namd-ucx/2.14
srun --mpi=pmi2 namd2 namd.in
~~~
{: .file-content}

### 7. Transferring equilibrated system between simulation packages.
Simulation packages have different methods and performance. It is useful to be able to transfer a running simulation from one software to another.

```
Examples:
Constant pH witn replica exchange - Amber
Targeted molecular dynamics - namd
Custom forces - namd
```

#### 7.1. Moving simulation from NAMD to AMBER.

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
ml StdEnv/2016.4 nixpkgs/16.09  gcc/7.3.0  cuda/9.2.148  openmpi/3.1.2 amber/18.10-18.11 scipy-stack

pmemd.cuda -O -i pmemd_prod.in -o production.log -p ../../prmtop.parm7 -c restart.rst7
~~~
{:.bash}

PMEMD is highly optimized to do all computations in one GPU, and it runs exceptionally fast. It CANNOT be used efficiently on more than one GPU because of the overhead from moving data between GPUs.

- Note: pmemd.cuda is unstable with nfft = 144, but stable with nfft = 128 or 256.

References:
[Amber file formats](https://ambermd.org/FileFormats.php#restart)


#### 5.2 Moving simulation from AMBER to GROMACS.

To tansfer simulation to GROMACS in addition to converting restart file we need to convert topology.

First convert AMBER topology to GROMACS

~~~
module --force purge
module load StdEnv/2020 gcc ambertools python scipy-stack
source $EBROOTAMBERTOOLS/amber.sh
python
~~~
{: .bash}

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
gmx mdrun
~~~
{: .bash}



