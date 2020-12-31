---
title: "Preparing a complex protein-RNA system for simulation"
teaching: 30
exercises: 5
questions:
- "How to add missing segments to a protein?"
- "How to add missing segments to a nucleic acid?"
- "How to align molecules with VMD?"
objectives:
- "?"
keypoints:
- "?"
---
For this workshop, we have chosen a complex of human argonaute-2 (hAgo2) protein with a micro RNA (miRNA) bound to a target messenger RNA (mRNA). miRNAs are short non-coding RNAs that are critical for regulating gene expression and the defense against viruses. miRNAs regulate a wide variety of human genes. They can control the production of proteins by targeting and inhibiting mRNAs. miRNAs can specifically regulate individual proteins' expression, and their selectivity is based on sequence complementarity between miRNAs and mRNAs. miRNAs that target mRNAs encoding oncoproteins can serve as selective tumor suppressors. They can inhibit tumor cells without a negative impact on all other types of cells. The discovery of this function of miRNAs has made miRNAs attractive tools for new therapeutic approaches. However, it is challenging to identify the most efficient miRNAs that can be targeted for medicinal purposes. To regulate protein synthesis miRNAs interact with hAgo2 protein forming the RNA-induced silencing complex that recognizes and inhibits the target mRNAs by slicing them. Therefore, elucidating the structural basis of the molecular recognition between hAgo2 and mRNA is crucial for understanding miRNA functions and developing new therapeutics for diseases.

Create working directory:
~~~
mkdir ~/scratch/workshop
cd ~/scratch/workshop
~~~
{: .bash}

### 1. Preparing protein for MD simulation
#### 1.1 Adding missing residues to protein structure files.

Almost all protein and nucleic acid crystallographic structure files have missing residues. The reason for it is that the most flexible parts of biopolymers are disordered in crystals, and therefore, the positions cannot be accurately determined. These atoms, however, may be crucial for MD simulations (e.g., loops connecting functional domains, nucleic acid chains, incomplete amino acid side chains ... etc.). For realistic simulation, we need to build a model contacting all atoms.

##### 1.1.1 Adding missing residues using SWISS MODEL

If you have bash and vmd on your computer, you can download and prepare pdb and sequence files locally. Otherwise do it on a CC system and download the results to your computer for insertion of missing fragments.

Download structure and sequence files from PDB database:
~~~
$ wget https://files.rcsb.org/download/6n4o.pdb
$ wget https://www.rcsb.org/fasta/entry/6N4O/download -O 6n4o.fasta
~~~
{: .bash}

Extract the full sequence of chain A:
~~~
$ grep -A1 "Chain A" 6n4o.fasta > 6n4o_chain_A.fasta
~~~
{: .bash}

Extract chain A from 6n4o.pdb using VMD. First, load the VMD module and launch VMD:
~~~
$ module load vmd
$ vmd
~~~
{: .bash}

In the VMD prompt, execute the following commands:
~~~
mol new 6n4o.pdb
set sel [atomselect top "chain A"]
$sel writepdb 6n4o_chain_A.pdb
quit
~~~
{: .vmd}

Download 6n4o_chain_A.pdb and 6n4o_chain_A.fasta to your computer for homology modeling with SWISS-MODEL.

- In a browser on your local computer navigate to [SWISS-MODEL](https://swissmodel.expasy.org) website.
- Click Start Modelling,
- Click User Template,
- Paste the full sequence of your protein or upload 6n4o_chain_A.fasta,
- Upload your structure file 6n4o_chain_A.pdb missing residues,
- Click Build Model,
- Download the homology model, and rename it to 6N4O_SWISS_PROT_model_chainA.pdb.

> ## Inspecting the model
> ~~~
> Compare the model with your original structure. Were all missing residues added?
> ~~~
{: .challenge}


##### 1.1.2 Adding missing residues using i-TASSER
The limitation of SWISS-MODEL server is that it is not capable of modeling long terminal fragments. Another homology modeling server [i-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/) (Iterative Threading ASSEmbly Refinement) uses the advanced protocol and is capable of predicting folding without any structural input. The downside of i-TASSER is that the process is much longer (about 60 hours for protein like 6n4o). In addition, i-TASSER optimizes, positions of all atoms, which is great, but sometimes not desirable. We can not wait for i-TASSER modeling to complete, but the result is available in the workshop data tarball.

##### 1.1.3. Preparing working directory.
Login to one of the CC systems, download data and unpack it in ~/scratch directory:
~~~
cd ~/scratch
wget md_workshop_data.tar.gz
tar -xf md_workshop_data.tar.gz
~~~
{: .bash}

Upload protein models from your computer:
~~~
scp 6N4O_SWISS_PROT_model_chainA.pdb \
6N4O_i-TASSER_model_chainA.pdb \
someuser@graham.computecanada.ca:scratch/workshop
~~~
{: .bash}
~~~
6N4O_SWISS_PROT_model_chainA.pdb              100%  516KB   1.2MB/s   00:00
6N4O_i-TASSER_model_chainA.pdb                100% 1082KB   5.3MB/s   00:00
~~~
{: .output}
In this command, "\\" is the line continuation character. Ensure that there is no whitespace characters after it.

#### 1.2. Aligning protein models.
i-TASSER procedure changes the orientation of the protein and slightly optimizes the positions of all atoms. We will keep the original atom positions and take only the terminal end from the i-TASSER model. To combine the i-TASSER model with the actual 6n4o coordinates, we need to align these two structures.
For alignment we can use only residues present in 6n4o.pdb. Missing residues are specified in the PDB file header in section "REMARK 465"
~~~
grep "REMARK 465" 6n4o.pdb
~~~
{: .bash}
~~~
REMARK 465 MISSING RESIDUES
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
REMARK 465
REMARK 465   M RES C SSSEQI
REMARK 465     MET A     1
REMARK 465     TYR A     2
REMARK 465     SER A     3
REMARK 465     GLY A     4
REMARK 465     ALA A     5
REMARK 465     GLY A     6
REMARK 465     PRO A     7
...
~~~
{: .output}

Using this information we can select all residues present in 6N4O.pdb:

~~~
atomselect top "resid 22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858"
~~~
{: .vmd}

You don't need to run this command right now, will use it later for alignment of structures and constrained energy minimization.

To begin the alignment process navigate to the working directory that you created on graham, and ensure that you have two files in the working directory:
~~~
cd  ~/scratch/workshop
ls
~~~
{: .bash}
~~~
6N4O_i-TASSER_model_chainA.pdb	6N4O_SWISS_PROT_model_chainA.pdb
~~~
{: .output}

Ensure that both files have only protein atoms (chain A).

Launch VMD and load two pdb files, they will be loaded as molecules 0 and 1:
~~~
mol new 6N4O_SWISS_PROT_model_chainA.pdb
mol new 6N4O_i-TASSER_model_chainA.pdb
~~~
{: .vmd}
Save the list of all residues present in 6N4O.pdb in the variable 6n4o_residues:
~~~
set 6n4o_residues "22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858"
~~~
{: .vmd}
Select residues defined in the variable *6n4o_residues* from both models, and save them in the variables *swissmodel* and *itasser*
~~~
set swissmodel [atomselect 0 "backbone and resid $6n4o_residues"]
set itasser [atomselect 1 "backbone and resid $6n4o_residues"]
~~~
{: .vmd}
Compute the transformation matrix *TransMat*
~~~
set TransMat [measure fit $itasser $swissmodel]
~~~
{: .vmd}
Select all residies of molecule 1 and apply the transformation matrix to the selection
~~~
echo rmsd before fit = [measure rmsd $itasser $swissmodel]
set itasser_all [atomselect 1 "all"]
$itasser_all move $TransMat
echo rmsd after fit = [measure rmsd $itasser $swissmodel]
~~~
{: .vmd}
Select residues 1-21 from molecule 1 and save them in the file 6n4o_resid_1-21.pdb
~~~
set terminal [atomselect 1 "noh resid 1 to 21"]
$terminal writepdb 6n4o_resid_1-21.pdb
quit
~~~
{: .vmd}
Combine the i-TASSER model of residues 1-21 with the SWISS-MODEL.
~~~
grep -h ATOM 6n4o_resid_1-21.pdb 6N4O_SWISS_PROT_model_chainA.pdb > 6n4o_chain_A_complete.pdb
~~~
{: .bash}

#### 1.3. Mutating residues
PDB entry 6N4O is the structure of the catalytically inactive hAgo2 mutant D669A. To construct the active form, we need to revert this mutation.

To accomplish this, we need to delete from ALA669 all atoms that are not present in ASP. Then change the residue name of ALA669 to ASP.

Let's begin by checking what atoms are present in residue 669:
~~~
grep 'A 669' 6n4o_chain_A_complete.pdb
~~~
{: .bash}
~~~
ATOM   5161  N   ALA A 669     -19.332  25.617 -27.862  1.00  0.97           N
ATOM   5162  CA  ALA A 669     -18.951  24.227 -27.916  1.00  0.97           C
ATOM   5163  C   ALA A 669     -17.435  24.057 -28.043  1.00  0.97           C
ATOM   5164  O   ALA A 669     -16.661  25.018 -28.091  1.00  0.97           O
ATOM   5165  CB  ALA A 669     -19.720  23.530 -29.053  1.00  0.97           C
~~~
{: .output}

If you are familiar with aminoacid structures, you remember that the alanine sidechain is made of only one beta carbon atom (CB). All amino acids except glycine have beta carbon as well. So there is nothing to delete. All we need to do is to change the resName of all five ALA669 atoms to ASP.
You can do it using stream editor:
~~~
sed 's/ALA A 669/ASP A 669/g' 6n4o_chain_A_complete.pdb > 6n4o_chain_A_complete_A669D.pdb
~~~
{: .bash}

Verify the result:
~~~
grep 'A 669' 6n4o_chain_A_complete_A669D.pdb
~~~
{: .bash}

~~~
ATOM   5161  N   ASP A 669     -19.332  25.617 -27.862  1.00  0.97           N
ATOM   5162  CA  ASP A 669     -18.951  24.227 -27.916  1.00  0.97           C
ATOM   5163  C   ASP A 669     -17.435  24.057 -28.043  1.00  0.97           C
ATOM   5164  O   ASP A 669     -16.661  25.018 -28.091  1.00  0.97           O
ATOM   5165  CB  ASP A 669     -19.720  23.530 -29.053  1.00  0.97           C
~~~
{: .output}

#### 1.4. Adding functionally important ions.
The catalytic site of hAgo2 is comprised of the three amino acids D597, E637, and D669. It is known that hAgo2 requires a divalent metal ion near the catalytic site to slice mRNA. The 6n4o PDB file does not have this ion, but another hAgo2 structure, 4w5o, does. We can align these two structures as we did in section 3 and then copy the Mg2+ ion located near the catalytic site from 4W5O to our model.

Download 4w5o.pdb
~~~
wget https://files.rcsb.org/download/4w5o.pdb
~~~
{: .bash}
Align 4w5o with 6n4o and save MG ions. For the alighnment we will use residues closest to the MG ions.
~~~
mol new 6n4o.pdb
mol new 4w5o.pdb
set 6n4o [atomselect 0 "backbone and resid 597 637 669 807"]
set 4w5o [atomselect 1 "backbone and resid 597 637 669 807"]
set TransMat [measure fit $4w5o $6n4o]
echo rmsd before fit = [measure rmsd $6n4o $4w5o]
set 4w5o_all [atomselect 1 "all"]
$4w5o_all move $TransMat
echo rmsd after fit = [measure rmsd $6n4o $4w5o]
set mg [atomselect 1 "resname MG"]
$mg set resid [$mg get residue]
$mg writepdb 4w5o_MG_ions.pdb
quit
~~~
{: .vmd}


#### 1.5. Assigning Protonation States to Residues
One of the important jobs of setting up a simulation system is assigning the protonation states and most likely tautomers of the HIS residues. TYR, LYN, CYS, and ARG are almost always in their standard ptotonation states at physiological pH. But you should decide for each GLU, ASP, and HIS which state is most likely.

You can use the H++ server to calculate pKa of titratable sites and select protonation states as described in Episode 6, "Assigning Protonation States to Residues in a Protein".

Force fields used by H++:
1. AMBER Lipid 2011 Force Field, lipid11.dat
2. PARM99 + frcmod.ff99SB + frcmod.parmbsc0 + OL3 for RNA, parm10.dat + frcmod.ff14SB
3. AMBER force-field parameters for phosphorylated amino acids, frcmod.phosaa10
4. Ions, see Leap log for more details

It does not yet supports 5' phosphorylated AA, so 5' phosphates must be removed.
TER records must be added.

Still does not work with RNA:
```
FATAL:  Atom .R<A 869>.A<HO3' 34> does not have a
type.
```

Does not make connection P-O3' because they are too far apart?
Temporary pdb file is created incorrectly with O3' of some added residues protonated because they are too far apart from P.

It looks like simrna does not insert fragments properly.
To Do:
Can RNA model be loaded in Leap?

### 2. Adding missing segments to RNA structure files.

First, we need to create a PDB file containing all RNA atoms placed at proper positions. At this initial step, we are not particularly concerned with the quality of the 3D structure because we will refine it afterward.

We can insert the missing residues using the freely available [ModeRNA server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone ModeRNA software. The automatic process used by ModeRNA server moves residues adjacent to the inserted fragment. Besides, modeRNA server offers a limited set of options, and the output PDB files will need more processing steps afterward. Changing atomic positions is not desirable because we want to keep all experimental coordinates. For these reasons, we will use the standalone modeRNA package. You can install modeRNA on CC systems, or if you are comfortable with installation of Python, you can install it on your computer.

#### 2.1. Installing ModeRNA on your computer

- Install python/2.7.14, numpy/1.11.3, biopython/1.58
- [Download ModerRNA](http://genesilico.pl/moderna/download/) and follow [installation instructions](http://genesilico.pl/moderna/installing/).

#### 2.3. Installing ModeRNA on CC systems.
~~~
module load StdEnv/2016.4 python/2.7.14
virtualenv ~/e27
source ~/e27/bin/activate
pip install numpy==1.11.3 biopython==1.58 ModeRNA==1.7.1
~~~
{: .bash}

Installation is required only once. When you login into your account next time you only need to activate the environment:

~~~
module load StdEnv/2016.4 python/2.7.14
source ~/e27/bin/activate
~~~
{: .bash}

As ModeRNA can only model a single RNA strand, we will model chains C and D separately. For this task we will need to prepare structural templates and sequence alignment files.

#### 2.4. Preparing structural templates for chains C and D.

Download 6n4o.pdb
~~~
wget https://files.rcsb.org/download/6n4o.pdb
~~~
{: .bash}
Residue 6 in chain D has only phosphate atoms and thus can not be used as a template. To prevent modeRNA from using it, we need to delete residue 6. Leaving it in the PDB file will lead to unwanted results.

~~~
mol new 6n4o.pdb
set sel [atomselect top "chain C or (chain D and not resid 6)"]
$sel writepdb 6n4o_chains_CD.pdb
quit
~~~
{: .vmd}

We created the file 6n4o_chains_CD.pdb suitable for use as a structural template.


#### 2.5. Preparing sequence alignment files for chains C and D.
Use a text editor of your choice (for example, nano or vi) to create two sequence alignment files, 6n4o_C.fasta and 6n4o_D.fasta.  The first file is for chain C and the latter is for chain D. Each file must contain two sequences, the sequence of the model to be built and the template sequence. The contents of the files is shown below.

Sequence alignment file for chain C:
~~~
cat 6n4o_C.fasta
~~~
{: .bash}
~~~
>Model
UGGAGUGUGACAAUGGUGUUU
>Template
UGGAGUGUG-CAAUGGUG-UU
~~~
{: .output}

Sequence alignment file for chain D:
~~~
cat 6n4o_D.fasta
~~~
{: .bash}
~~~
>Model
CAUUGUCACACUCCAAAA
>Template
CAUUG---CACUCCAA--
~~~
{: .output}

#### 2.6. Inserting missing segments.
Once you install modeRNA program, you will be able to use all its functions. Below are the commands needed to insert missing fragments in chains C and D.  Description of all commands is available [here](http://genesilico.pl/moderna/commands/).
~~~
from moderna import *
# Model chain C
tC = load_template('6n4o_chains_CD.pdb', 'C')
aC = load_alignment('6n4o_C.fasta')
mC = create_model(model_chain_name = 'A')
apply_alignment(tC, aC, mC)
apply_indel(mC, '9', '11', 'A')
apply_indel(mC, '18', '20', 'U')
apply_missing_ends(aC, mC)
renumber_chain(mC, '1')
write_model(mC, 'chain_C_model_A.pdb')
# Model chain D
tD = load_template('6n4o_chains_CD.pdb', 'D')
aD = load_alignment('6n4o_D.fasta')
mD = create_model(model_chain_name = 'B')
apply_alignment(tD, aD, mD)
apply_indel(mD, '5', '9', 'UCA')
apply_missing_ends(aD, mD)
renumber_chain(mD, '1')
write_model(mD, 'chain_D_model_B.pdb')
~~~
{: .python}

You can start python and execute the commands interactively:
~~~
module load StdEnv/2016.4 python/2.7.14
source ~/e27/bin/activate
python
~~~
{: .bash}
Or save these commands in the file make_models.py and run it non-iteractively:
~~~
module load StdEnv/2016.4 python/2.7.14
source ~/e27/bin/activate
python make_models.py
~~~
{: .bash}
The modeRNA program will create two model files, chain_C_model_A.pdb, and chain_D_model_B.pdb. Ensure that the program added all missing residues and did not move any residues present in the original structure file.

> ## How will ModeRNA server do the same task?
> Try using ModeRNA server. Compare automatically generated ModeRNA models with the original chains C and D. Did server only added missing residues without moving any other atoms?
>
{: .challenge}

#### 2.7. Pitfalls with structure generated automatically by ModeRNA server.
ModeRNA server inserts missing residues with the same residue number as the residue before insertion. It uses an insertion code (iCode) to mark inserted residues. For example, if residues 6-8 are missing, the program will insert them as residues 5A, 5B, and 5C. It is not possible to instruct ModeRNA webserver to renumber residues, and it is not possible to change chainID either. So you will need to take care of residue numbering and renaming chains manually.

**References:**
1. [ModeRNA: a tool for comparative modeling of RNA 3D structure](https://doi.org/10.1093/nar/gkq1320)
2. [PDB file format](
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)

### 3. Optimizing double stranded RNA.
As we built two complementary RNA strands independently of each other, they may clash when two strands are combined. To resolve clashes and obtain a good duplex structure, we need to minimize the energy of the whole RNA duplex.

We will minimize energy of the double stranded RNA using SimRNA. It is available as [SimRNAweb server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone  [SimRNA software] (https://genesilico.pl/SimRNAweb). This program allows for RNA 3D structure modeling with optional restraints.

SimRNA features:
- Coarse-grained representation (5 atoms per residue)
- Monte Carlo method for sampling the conformational space
- The energy function is composed of statistical potential terms, derived from the observed frequencies of occurrence of various proximate structural patterns.

#### 3.1. Installation of the SimRNA binary package:
Simulation submitted to public SimRNAweb server may wait up to a few days in the queue, while on a local computer, you can do it in a couple of minutes. SimRNA is available as a binary distribution, so no installation is required. You only need to download and unpack the package:
~~~
cd ~
wget https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz --no-check-certificate
tar -xf SimRNA_64bitIntel_Linux.tgz
~~~
{: .bash}

Then go back to the working directory.

#### 3.2. Preparing input files for simulation with SimRNAweb server.
If you will be using standalone SimRNA program you can skip this sections and proceed to the next one.

SimRNAweb server requires the user to provide RNA sequence, a list of residues not allowed to move in simulation and RNA structure in PDB format.

Sequence:
~~~
UGGAGUGUGACAAUGGUGUUU CCAUUGUCACACUCCAAA
~~~
{: .string}
The convention is that the first and the second sequences correspond to chains A and B in the PDB structure file.

List of residues not allowed to move (we don't want the program to move atoms resolved in the experimental structure).
~~~
 A:1-9,11-18,20-21;B:1-5,9-16
~~~
{: .string}
PDB file matching the sequence. All atoms must be present, and chains must be named A and B. To prepare this file combine chains A and B:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

You can use chains_CD_model_AB.pdb only for SimRNAweb simulation. For simulation with standalone SimRNA program you will need to modify the PDB structure file as described in the next section.

#### 3.3. Preparing structure file for simulation with standalone SimRNA program.

Command-line SimRNA program does not take a list of frozen atoms as a separate input. Instead, you need to insert this information into the PDB structure file. But before we flag frozen atoms the missing phosphate must be added to the 5' terminal residue of chain D.

We begin with merging chains A and B if you have not done this yet:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

Next we apply two modifications to this file. First, we need to add phosphate to the 5' terminal residue of chain D. SimRNA expects all residues to have a P atom. SimRNAweb will add P automatically, but for simulation with a standalone SimRNA program, we need to do it manually. There are several options to add the phosphate.

##### 3.3.1. Renaming O5' atom to P
The most straightforward fix is to rename O5' atom to P. if you chose to do this, save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb, and skip the next step.

##### 3.3.2. Adding 5' monophosphate with AmberTools/20.
First, rename the phosphorylated 5' terminal nucleotide according to AMBER convention. The names of phosphorylated terminals in AMBER are A5, C5, G5, U5, DA5, DC5, DG5, and DT5. Libraries of phosphorylated 5' terminal nucleotides are in the file '$AMBERHOME/dat/leap/lib/terminal_monophosphate.lib'.

Load the AmberTools module:
~~~
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 ambertools/20
source $EBROOTAMBERTOOLS/amber.sh
~~~
{: .bash}

Watch for the message printed on screen when the module is loaded. The message instructs how to use the module.

Launch Leap and load RNA force field:
~~~
tleap -f leaprc.RNA.OL3 -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}
In the Leap promt execute the commands:
~~~
loadoff terminal_monophosphate.lib
chainD = loadpdb chain_D_model_B.pdb
savepdb chainD chain_D5P.pdb
quit
~~~
{: .leap}
These commands will load libraries of phosphorylated 5' terminal nucleotides and chain D PDB file. Leap will automatically add all missing atoms based on library entries and save chainD in PDB file chain_D5P.pdb:

We don't want to use the PDB file prepared with Leap for SimRNA because AMBER has different aminoacid naming conventions. So we copy phosphate atoms from chain_D5P.pdb and paste them into chains_CD_model_AB.pdb. We then edit chain ID, residue ID, and residue name. Save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb

##### 3.3.3. Adding 5' monophosphate with [CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader).
You can also add phosphate using CHARMM-GUI. Beware that CHARMM-GUI changes residue names to the old-style RNA 3-letter names and changes chain ID to "R".

##### 3.3.4 Define frozen atoms.
Standalone SimRNA program accepts PDB file where frozen atoms have occupancy 0.0 and completely free have occupancy 1.0. You can change values of the occupancy with the following VMD commands:
~~~
mol new chains_CD_model_AB_5P.pdb
set sel [atomselect top all]
$sel set occupancy 0
set sel [atomselect top "chain A and resid 10 19"]
$sel set occupancy 1
set sel [atomselect top "chain B and resid 6 7 8 17 18"]
$sel set occupancy 1
set sel [atomselect top all]
$sel writepdb chains_CD_model_AB_5P_frozen.pdb
quit
~~~
{: .vmd}

#### 3.4. Running simulation
SimRNA needs two files in the working directory:
'chains_CD_model_AB_5P_frozen.pdb' and 'config', the file with SimRNA simulation parameters.

Example config file:
~~~
NUMBER_OF_ITERATIONS 1000000
TRA_WRITE_IN_EVERY_N_ITERATIONS 10000
INIT_TEMP 1.35
FINAL_TEMP 0.90

BONDS_WEIGHT 1.0
ANGLES_WEIGHT 1.0
TORS_ANGLES_WEIGHT 0.0
ETA_THETA_WEIGHT 0.40
~~~
{: .string}

In the working directory, make a symbolic link to the 'data' directory located in SimRNA distribution. Assuming that you installed SimRNA in $HOME the link command is:
~~~
ln -s data ~/SimRNA_64bitIntel_Linux/data
~~~
{: .bash}
Then run the simulation:
~~~
srun -A <desired account> -c10 --mem-per-cpu=1000 --time=30:0 \
~/SimRNA_64bitIntel_Linux/SimRNA \
-P chains_CD_model_AB_5P_frozen.pdb \
-c config -E 10
~~~
{: .bash}

The option -E \<number of replicas> turns on replica exchange mode.
Replica exchange mode is parallelized with OMP. Each replica simulation can run on its own CPU independently of others,  so for the optimal performance allocate the same number of cores (option -c) as the number of replicas (option -E).

The simulation will run for about two minutes and produce trajectory file *.trafl for each replica.

#### 3.5. Processing simulation trajectory
The simplest way of processing the trajectory files is obtaining the lowest energy structure. Generally, better results can be obtained by using clustering. Clustering tool is included with the distribution, but using clustering is outside the scope of this workshop.

Extract the lowest energy frame from the trajectory of the first replica
~~~
~/SimRNA_64bitIntel_Linux/trafl_extract_lowestE_frame.py \
chains_CD_model_AB_5P_frozen.pdb_01.trafl
~~~
{: .bash}
Convert the lowest energy frame to PDB format
~~~
~/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs chains_CD_model_AB_5P.pdb \
chains_CD_model_AB_5P_frozen.pdb_01_minE.trafl 1 AA
~~~
{: .bash}
This command will create PDB file of the lowest energy structure from trajectory of replica 1:
chains_CD_model_AB_5P_frozen.pdb_01_minE-000001_AA.pdb
We will use this relaxed structure for simulation. Rename it into a shorter name, for example chains_CD_minimized.pdb


**References**
1. [SimRNA: a coarse-grained method for RNA folding simulations and 3D structure prediction](https://doi.org/10.1093/nar/gkv1479)
2. [SimRNA manual](https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_UserManual_v3_20_20141002.pdf)
3. [VMD TCL commands](https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node121.html)


### 4. Preparing simulation system.

#### 4.1 Determine the number of water and salt molecules needed to prepare solvated system.
To prepare solution with the desired ionic strength we will use SLTCAP server. For this calulation we need to know the molecular weight of the macromolecules, their charge, and the number of water molecules in the simulation system.

The molecular weight of hAgo2 is 97,208 Da, and the MW of our nucleic acids is 12.5 KDa [[calculate MW of the RNA]](http://www.encorbio.com/protocols/Nuc-MW.htm). Thus, the total MW is 110 KDa.

To determine the number of water molecules we will solvate the system in a cubic box extending 13 A from the solute.

Launch Leap and load protein and RNA forcefields:
~~~
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
source $EBROOTAMBERTOOLS/amber.sh
tleap -f leaprc.RNA.OL3 -f leaprc.protein.ff14SB -f leaprc.water.tip3p -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}

Load pdb files into leap, combine them, and solvate the system
~~~
loadoff terminal_monophosphate.lib
rna = loadpdb chains_CD_minimized.pdb
prot = loadpdb 6n4o_chain_A_complete_A669D.pdb
mg = loadpdb 4w5o_MG_ions.pdb
sys = combine {prot,rna,mg}
solvatebox sys TIP3PBOX 13 iso
charge sys
quit
~~~
{: .leap}

~~~
  Solute vdw bounding box:              110.730 72.051 85.804
  Total bounding box for atom centers:  136.730 136.730 136.730
      (box expansion for 'iso' is  70.5%)
  Solvent unit box:                     18.774 18.774 18.774
  Volume: 2740120.355 A^3
  Total mass 1460917.373 amu,  Density 0.885 g/cc
  Added 74991 residues.
> charge sys
  Total unperturbed charge:  -4.000000
  Total perturbed charge:    -4.000000
~~~
{: .output}

Using this information (MW 110 KDa, charge -4.0, 75000 water molecules) as an input to [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server we obtain the number of ions: 188.64 anions and 192.64 cations.

#### 4.2 Determine protonation states of titratable sites.
For processing with H++ server we need to merge protein, nucleic acids and ions into one PDB file.

~~~
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
source $EBROOTAMBERTOOLS/amber.sh
mkdir ~/scratch/Ago-RNA_sim/HPP
cd ~/scratch/Ago-RNA_sim/HPP
tleap -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}

H++ server does not have library entries for phosphorylated 5' terminals. To workaround we simply remove 5' phosphate atoms.

~~~
source leaprc.RNA.OL3
source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadoff terminal_monophosphate.lib
rna=loadpdb ../prep_system/chains_CD_minimized.pdb
remove rna rna.1@P
remove rna rna.1@OP1
remove rna rna.1@OP2
remove rna rna.1@OP3
remove rna rna.1@HP3
remove rna rna.22@P
remove rna rna.22@OP1
remove rna rna.22@OP2
remove rna rna.22@OP3
remove rna rna.22@HP3
prot=loadpdb ../prep_system/6n4o_chain_A_complete_A669D.pdb
mg=loadpdb ../prep_system/4w5o_MG_ions.pdb
sys=combine {prot rna mg}
savepdb sys 6n4o_Hpp.pdb
quit
~~~
{: .leap}

Process 6n4o_Hpp.pdb with H++ server. Uncheck 'Correct orientation' in caclulation setup. When calculation completes download the list of computed pKs (0.15_80_10_pH6.5_6n4o_Hpp.pkout.txt)

Examine HIS, ASP, GLU.
~~~
grep ^HI 0.15_80_10_pH6.5_6n4o_Hpp.pkout.txt
grep ^AS 0.15_80_10_pH6.5_6n4o_Hpp.pkout.txt
grep ^GL 0.15_80_10_pH6.5_6n4o_Hpp.pkout.txt
~~~
{: .bash}

> ## Determine protonation states from the list of pKs.
> What titratable sites in 6n4o are in non-standard protonation state at pH 6.5?
>> ## Solution
>>Histidines 77, 766, 822, and 829 are protonated (HIP)
> {: .solution}
{: .challenge}

#### 4.3 Preparing the complete simulation system

Finally, we are ready to prepare the complete simulation system. We can run all the commands interactively, or save them in a file and then execute it.

Leap was designed to read commands from a file (-f option). This means that we need two scripts: one with the leap commands, and another with commands to run leap itself.

Taking advantage of shell flexibility, we can eliminate two files' necessity by creating a multiline variable holding all commands and then passing this variable instead of file to leap.

Load phosphorylated aminoacids forcefield:
frcmod.phosaa14SB

OP-P -OH        90.0    108.23  from O2-P -OH ! raised to prevent simulation instabilities
HO-OH-P         100.0   113.28  from HO-OH-P  ! raised to prevent simulation instabilities, equilibrium from G03 calculations


As Leap does not support input from STDIN we will use the <(echo "$inputData") syntax which provides a way to pass the output of a command (echo "$inputData") to a program that can not use pipeline.

~~~
#!/bin/bash
# FILE <<< prep_system.leap >>>
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
source $EBROOTAMBERTOOLS/amber.sh

inputData=$(cat << EOF
loadamberparams frcmod.monophosphate
loadoff terminal_monophosphate.lib
rna = loadpdb chains_CD_minimized.pdb
prot = loadpdb 6n4o_chain_A_complete_A669D.pdb
mg = loadpdb 4w5o_MG_ions.pdb
sys = combine {prot,rna,mg}
set {sys.77 sys.766 sys.822 sys.829} name "HIP"
addions sys Na+ 0
solvatebox sys TIP3PBOX 13 iso
addionsrand sys Na+ 189 Cl- 189
saveamberparm sys prmtop.parm7 inpcrd.rst7
savepdb sys inpcrd.pdb
quit
EOF)

tleap -f leaprc.RNA.OL3 -f leaprc.protein.ff14SB -f leaprc.water.tip3p -I $EBROOTAMBERTOOLS/dat/leap/lib/ -f <(echo "$inputData")
~~~
{:.file-content}


### 5. Energy minimization.

First we need to optimize positions of atoms. To restrain residues present in 6n4o we need to select all residues that have coordinates in the pdb file, but as in the simulation system resisues are renumbered we can not use the original numbers. All residues in simulation systems are counted sequentially starting from first to last.

Chain       | Original | Shift | Simulation |
------------|----------|-------|------------|
Protein     | 1:859    |  -    | 1:859      |
RNA chain C | 1:21     | 859   | 860:880    |
RNA chain D | 1:18     | 898   | 881:898    |
MG ions     |    -     |  -    | 899:901    |

Considering this,  selection command in VMD will be:
~~~
atomselect 0 "backbone and resid 22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858 860 to 868 870 to 877 879 to 885 889 to 896"
~~~
{: .vmd}

#### 5.1 Energy minimization with AMBER

~~~
sander -O -i min.in -p ../prmtop.parm7 -c ../inpcrd.rst7  -ref ../inpcrd.rst7
~~~
{: .bash}
~~~
...
LINMIN FAILURE
...
~~~
{: .output}

Minimization fails.


#### 5.2 Energy minimization with  NAMD

In the first roud of minimization we will not allow all original atoms to move. The only exception is ARG814 which is in close contact with the added RNA segment. Prepare constraints file for this run:

~~~
mol new prmtop.parm7
mol addfile inpcrd.rst7
set sel [atomselect top "all"]
$sel set occupancy 0.0
set sel [atomselect top "noh resid 22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858 860 to 868 870 to 877 879 to 885 889 to 896"]
$sel set occupancy 999.9
set sel [atomselect top "resid 814"]
$sel set occupancy 0.0
set sel [atomselect top "all"]
$sel writepdb constrain_all_6n4o_residues.pdb
quit
~~~
{: .vmd}

Run 1000 steps of energy minimization:
~~~
module load StdEnv/2020 intel/2020.1.217 namd-multicore/2.14
charmrun ++local +p 8 namd2 namd_min1.in >&log&
~~~
{:.bash}

In the second round of minimization we will constrain only backbone atoms of all original residues using the minimized coordinates as reference. Prepare the reference coordinates:

~~~
module load StdEnv/2020 intel vmd
cd ~/scratch/Ago-RNA_sim/sim_namd/1-minimization
vmd
~~~
{:.bash}

~~~
mol new ../../prmtop.parm7
mol addfile minimized.coor
set sel [atomselect top "all"]
$sel writepdb ../../minimized.pdb
quit
~~~
{: .vmd}

Prepare force constants file:

~~~
cd ~/scratch/Ago-RNA_sim
module purge
module load StdEnv/2020 intel vmd
vmd
~~~
{: .bash}

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
module load StdEnv/2020 intel/2020.1.217 namd-multicore/2.14
charmrun ++local +p 8 namd2 namd_min2.in >&log&
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
module load StdEnv/2016.4 nixpkgs/16.09  gcc/7.3.0  cuda/9.2.148  openmpi/3.1.2 amber/18.10-18.11 scipy-stack

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
module purge
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
Velocities in text .rst7 files are in 1/ps. For restart files they need to be scaled by 20.455

~~~
vel_rst7 = pmd.load_file('../sim_pmemd/2-production/vel.rst7')
amber.velocities = vel_rst7.coordinates*20.455
amber.save('restart.gro')
~~~
{: .python}

~~~
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gromacs
gmx trjconv -f restart.gro -o restart.trr
gmx make_ndx -f restart.gro
gmx grompp -p gromacs.top  -c restart.gro -t restart.trr -f gromacs_production.mdp
~~~
{: .bash}

Running simulation
~~~
#SBATCH --mem-per-cpu=4000 --time=10:0:0 -c16
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gromacs
gmx mdrun
~~~

References:
 [Lessons learned from comparing molecular dynamics engines on the SAMPL5 dataset](https://link.springer.com/article/10.1007%2Fs10822-016-9977-1)

### 6. Benchmarking

#### GROMACS

~~~
gmx convert-tpr -s topol.tpr -nsteps 10000 -o next.tpr
gmx mdrun -s next.tpr -cpi state.cpt
~~~
{: .bash}

~~~
#SBATCH --mem-per-cpu=4000 --time=10:0:0 -c4 --ntasks=2
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gromacs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
srun gmx mdrun -s next.tpr -cpi state.cpi
~~~
{: .file-content}

gmx_mpi, gromacs/2020.4, avx512, single
CPU | Tasks | Threads | ns/day
----|-------|---------|-------
16  |  16   |   1     | 7.13
32  |  32   |   1     | 10.99
40  |  40   |   1     | 13.64
64  |  64   |   1     | 20.76
128 |  128  |   1     | 36.56

gmx_mpi, gromacs/2020.4, avx2
CPU | Tasks | Threads | ns/day
----|-------|---------|-------
16  |  16   |   1     | 6.37
32  |  32   |   1     | 10.27
64  |  64   |   1     | 17.03
128 | 128   |   1     | 33.32

#### AMBER

pmemd.cuda
41.77 ns/day

#### NAMD

namd-ucx
CPU | ns/day |
----|--------|
80  |  4.67  |
160 |  9.06  |

namd-multicore-cuda
CPU | GPU | ns/day
----|-----|--------
 8  | 1   | 5.92
16  | 2   | 8.85
40  | 2   | 11.17

namd-multicore
CPU |  ns/day
----|---------
1   | 0.08163
2   | 0.15828
4   | 0.30568
8   | 0.59241
16  | 1.14872
32  | 2.17942

#### Accelerating simulation

It is possible to increase time step to 4 fs.

Increase step to 4 fs with Hydrogen mass repartitioning.

parmed prmtop.parm7

~~~
hmassrepartition
outparm prmtop_hmass.parm7
quit
~~~

Reference
[Long-Time-Step Molecular Dynamics through Hydrogen Mass Repartitioning](https://pubs.acs.org/doi/abs/10.1021/ct5010406)
