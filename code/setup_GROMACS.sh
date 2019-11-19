#!/bin/bash

PDB_file=1GOA.pdb
neutral_GLU="32 154"
neutral_ASP="10"
neutral_LYS="99"
neutral_ARG="none"
delta_HIS="none"
epsilon_HIS="none"
protonated_HIS="none"

mkdir md_system
cd md_system
wget http://files.rcsb.org/view/1GOA.pdb
grep ^ATOM $PDB_file > protein.pdb

module purge
module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3

echo "mol new protein.pdb" > mutate.vmd
if [ "$neutral_GLU" != "none" ]; then
 echo "set sel [atomselect top \"resid ${neutral_GLU}\"]" >> mutate.vmd
 echo "\$sel set resname GLH" >> mutate.vmd
fi
if [ "$neutral_ASP" != "none" ]; then
 echo "set sel [atomselect top \"resid ${neutral_ASP}\"]" >> mutate.vmd
 echo "\$sel set resname ASH" >> mutate.vmd
fi
if [ "$neutral_LYS" != "none" ]; then
 echo "set sel [atomselect top \"resid ${neutral_LYS}\"]" >> mutate.vmd
 echo "\$sel set resname LYN" >> mutate.vmd
fi
if [ "$neutral_ARG" != "none" ]; then
 echo "set sel [atomselect top \"resid ${neutral_ARG}\"]" >> mutate.vmd
 echo "\$sel set resname ARN" >> mutate.vmd
fi
if [ "$delta_HIS" != "none" ]; then
 echo "set sel [atomselect top \"resid ${delta_HIS}\"]" >> mutate.vmd
 echo "\$sel set resname HID" >> mutate.vmd
fi
if [ "$epsilon_HIS" != "none" ]; then
 echo "set sel [atomselect top \"resid ${epsilon_HIS}\"]" >> mutate.vmd
 echo "\$sel set resname HIE" >> mutate.vmd
fi
if [ "$protonated_HIS" != "none" ]; then
 echo "set sel [atomselect top \"resid ${protonated_HIS}\"]" >> mutate.vmd
 echo "\$sel set resname HIP" >> mutate.vmd
fi
cat << EOF >> mutate.vmd
set sel [atomselect top all]
\$sel writepdb protonated.pdb
quit
EOF
vmd -dispdev text -e mutate.vmd

module purge
module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
gmx pdb2gmx -f protonated.pdb  -ff amber99sb-ildn -water spce
gmx editconf -f conf.gro -o boxed.gro -c -d 1.0 -bt cubic
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
touch mdrun.mdp
gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
echo "SOL" | gmx genion -s solvated.tpr -p topol.top -nname CL -neutral -o neutralized.gro
