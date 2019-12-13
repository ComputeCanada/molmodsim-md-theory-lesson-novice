#!/bin/bash

# >>>>> Clean up and protonate PDB structure with VMD <<<<<
mkdir amber
cd amber
module purge
module load nixpkgs/16.09  intel/2016.4 vmd/1.9.3
# Save VMD commands in the variable 'script'
script="$(cat << EOF
# Load 1RGG.pdb into a new (top) molecule
mol urlload pdb "http://files.rcsb.org/view/1RGG.pdb"
# Select and save all chain A protein atoms
set s [atomselect top "protein and chain A"]
\$s writepdb 1RGG_chain_A.pdb
# Delete the top molecule
mol delete top
# Load chain A into a new molecule
mol new 1RGG_chain_A.pdb
# Protonate ASP79
set s [atomselect top "resid 79"]
\$s set resname ASH
# Protonate HIS53
set s [atomselect top "resid 53"]
\$s set resname HIP
# Rename cross-linked cysteins
set s [atomselect top "resid 7 96"]
\$s set resname CYX
# Select the base and the alternate locations
set s [atomselect top "(altloc '') or (altloc A and resid 6 13 42 85 91) or (altloc B and resid 5 54)"]
# Save the selection
\$s writepdb 1RGG_chain_A_prot.pdb
quit
EOF
)"
# Execute VMD commands
vmd -e <<<"${script}" /dev/stdin

# >>>>>  Prepare the system using LEaP + pbd2gmx <<<<<
module purge
module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.1.1 amber/18
# Save LeaP commands in the variable 'script'
script="$(cat << EOF
# Load spce water model
source leaprc.water.spce
# Load Force Field parameters
source leaprc.protein.ff14SB
# Load PDB structure file
s = loadpdb 1RGG_chain_A_prot.pdb
# Neutralize
addions s Na+ 0
# Solvate
solvatebox s SPCBOX 15 iso
# Add salt
addionsrand s Na+ 24 Cl- 24
# Make SS bond
bond s.7.SG s.96.SG
# Save molecular parameters for AMBER/NAMD
saveamberparm s prmtop inpcrd
savepdb s 1RGG_chain_A_solvated.pdb
quit
EOF
)"
# Execute LEaP commands
tleap -f  <<<"${script}"  /dev/stdin
cd ..
# Generate molecular parameters for GROMACS
mkdir gromacs_A
cp amber/1RGG_chain_A_solvated.pdb gromacs_A
cd gromacs_A
module purge
module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
# Rename ions for GROMACS and assign them to chain B
sed s/"Cl-  Cl-  "/" CL  CL  B"/g 1RGG_chain_A_solvated.pdb | sed s/"Na+  Na+  "/" NA  NA  B"/g > 1RGG_chain_A_solvated_gro.pdb
# Save GROMACS topology and coordinate files
gmx pdb2gmx -f 1RGG_chain_A_solvated_gro.pdb -ff amber99sb-ildn -water spce -ignh -chainsep id -ss << EOF > log
y
EOF
cd ..

# >>>>>  Prepare the system using GROMACS pdb2gmx <<<<<
mkdir gromacs_B
cp amber/1RGG_chain_A_prot.pdb gromacs_B
cd gromacs_B
module purge
module load gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
# Generate GROMACS topology and coordinate files
gmx pdb2gmx -f 1RGG_chain_A_prot.pdb -ff amber99sb-ildn -water spce -ignh -chainsep id -ss << EOF >log
y
EOF
# Add periodic box to conf.gro
gmx editconf -f conf.gro -o boxed.gro -c -d 1.5 -bt cubic
# Solvate the box
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
# Create MD run parameters file
touch mdrun.mdp
# Generate the binary topology 'solvated.tpr'
gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
# Add salt and neutralize the system
echo "SOL" | gmx genion -s solvated.tpr -p topol.top -conc 0.15 -neutral -o neutralized.pdb
