#!/bin/bash

PDB_file=1GOA.pdb
neutral_GLU=(155)
neutral_ASP=(10 70)
neutral_LYS=(100)
delta_HIS=(84)
epsilon_HIS=(62)
protonated_HIS=(115)

mkdir md_system
cd md_system
wget http://files.rcsb.org/view/1GOA.pdb
grep ^ATOM $PDB_file > protein.pdb

#module purge
#module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3

#echo "mol new protein.pdb" > mutate.vmd
#if [ "${neutral_GLU[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${neutral_GLU[@]}\"]" >> mutate.vmd
# echo "\$sel set resname GLH" >> mutate.vmd
#fi
#if [ "${neutral_ASP[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${neutral_ASP[@]}\"]" >> mutate.vmd
# echo "\$sel set resname ASH" >> mutate.vmd
#fi
#if [ "${neutral_LYS[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${neutral_LYS[@]}\"]" >> mutate.vmd
# echo "\$sel set resname LYN" >> mutate.vmd
#fi
#if [ "${neutral_ARG[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${neutral_ARG[@]}\"]" >> mutate.vmd
# echo "\$sel set resname ARN" >> mutate.vmd
#fi
#if [ "${delta_HIS[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${delta_HIS[@]}\"]" >> mutate.vmd
# echo "\$sel set resname HID" >> mutate.vmd
#fi
#if [ "${epsilon_HIS[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${epsilon_HIS[@]}\"]" >> mutate.vmd
# echo "\$sel set resname HIE" >> mutate.vmd
#fi
#if [ "${protonated_HIS[@]}" != "none" ]; then
# echo "set sel [atomselect top \"resid ${protonated_HIS[@]}\"]" >> mutate.vmd
# echo "\$sel set resname HIP" >> mutate.vmd
#fi
#cat << EOF >> mutate.vmd
#set sel [atomselect top all]
#$sel writepdb protonated.pdb
#quit
#EOF
#vmd -dispdev text -e mutate.vmd

# Prepend sys. to all elements of the arrays
neutral_GLU=( "${neutral_GLU[@]/#/sys.}" )
neutral_ASP=( "${neutral_ASP[@]/#/sys.}" )
neutral_LYS=( "${neutral_LYS[@]/#/sys.}" )
delta_HIS=( "${delta_HIS[@]/#/sys.}" )
epsilon_HIS=( "${epsilon_HIS[@]/#/sys.}" )
protonated_HIS=( "${protonated_HIS[@]/#/sys.}" )
# Load PDB file
echo "sys = loadpdb protein.pdb" > mutate.leap
# Rename titratable residues if needed
if (( ${#neutral_GLU[@]} )); then
 echo  "set {${neutral_GLU[@]}} name \"GLH\"" >> mutate.leap
fi
if (( ${#neutral_ASP[@]} )); then
 echo  "set {${neutral_ASP[@]}} name \"ASH\"" >> mutate.leap
fi
if (( ${#neutral_LYS[@]} )); then
 echo  "set {${neutral_LYS[@]}} name \"LYN\"" >> mutate.leap
fi
if (( ${#delta_HIS[@]} )); then
 echo  "set {${delta_HIS[@]}} name \"HID\"" >> mutate.leap
fi
if (( ${#epsilon_HIS[@]} )); then
 echo  "set {${epsilon_HIS[@]}} name \"HIE\"" >> mutate.leap
fi
if (( ${#protonated_HIS[@]} )); then
 echo  "set {${protonated_HIS[@]}} name \"HIP\"" >> mutate.leap
fi
echo savepdb sys protonated.leap.pdb >> mutate.leap
echo quit >> mutate.leap
module purge
module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.1.1 amber/18
tleap -f leaprc.protein.ff14SB -f mutate.leap


module purge
module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
gmx pdb2gmx -f protonated.leap.pdb  -ff amber99sb-ildn -ignh -water spce
#gmx editconf -f conf.gro -o boxed.gro -c -d 1.0 -bt cubic
#gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
#touch mdrun.mdp
#gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
#echo "SOL" | gmx genion -s solvated.tpr -p topol.top -nname CL -neutral -o neutralized.gro
