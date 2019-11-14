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


