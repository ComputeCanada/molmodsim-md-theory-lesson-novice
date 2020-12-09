---
title: "Preparing a protein-RNA complex for MD simulation"
teaching: 30
exercises: 5
questions:
- "?"
objectives:
- "?"
keypoints:
- "?"
---
For this workshop, we have chosen a complex of human argonaute-2 (hAgo2) protein with a micro RNA (miRNA) bound to a target messenger RNA (mRNA). miRNAs are short non-coding RNAs that are critical for regulating gene expression and the defense against viruses. miRNAs regulate a wide variety of human genes. They can control the production of proteins by targeting and inhibiting mRNAs. miRNAs can specifically regulate individual proteins' expression, and their selectivity is based on sequence complementarity between miRNAs and mRNAs. miRNAs that target mRNAs encoding oncoproteins can serve as selective tumor suppressors. They can inhibit tumor cells without a negative impact on all other types of cells. The discovery of this function of miRNAs has made miRNAs attractive tools for new therapeutic approaches. However, it is challenging to identify the most efficient miRNAs that can be targeted for medicinal purposes. To regulate protein synthesis miRNAs interact with hAgo2 protein forming the RNA-induced silencing complex that recognizes and inhibits the target mRNAs by slicing them. Therefore, elucidating the structural basis of the molecular recognition between hAgo2 and mRNA is crucial for understanding miRNA functions and developing new therapeutics for diseases.

Create working directory:

~~~
$ mkdir ~/scratch/workshop
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
vmd > mol new 6n4o.pdb
vmd > set sel [atomselect top "chain A"]
vmd > $sel writepdb 6n4o_chain_A.pdb
vmd > quit
~~~
{: .bash}

Download 6n4o_chain_A.pdb and 6n4o_chain_A.fasta to your computer for homology modeling with SWISS-MODEL and i-TASSER.

&nbsp;1. On your local computer navigate to [SWISS-MODEL](https://swissmodel.expasy.org)

&nbsp;2. Click Start Modelling,

&nbsp;3. Click User Template,

&nbsp;4. Paste the full sequence of your protein or upload 6n4o_chain_A.fasta,

&nbsp;5. Upload your structure file 6n4o_chain_A.pdb missing residues,

&nbsp;6. Click Build Model,

&nbsp;7. Download the homology model, and save it in the file '6N4O_SWISS_PROT_model_chainA.pdb'.

&nbsp;8. Compare the model with your original structure. Were all missing residues added?

In the following sections we will assume that the SWISS-MODEL is saved in the file **6N4O_SWISS_PROT_model_chainA.pdb**

##### 1.1.2 Other homology modeling servers
SWISS-MODEL server does not add terminal fragments. Another homology modeling server [i-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/) (Iterative Threading ASSEmbly Refinement) uses the advanced protocol and is capable of threading terminal fragments. The downside of i-TASSER is that the process is much longer (about 60 hours for protein like 6n4o). i-TASSER optimizes, positions of all atoms, which is not always desirable.

Save the result of the i-TASSER modeling in the file **6N4O_i-TASSER_model_chainA.pdb**

##### 1.1.3. Preparing working directory
Login to one of the CC systems
Create working directory and descend into it:
~~~
$ mkdir ~/scratch/workshop
$ cd ~/scratch/workshop
~~~
{: .bash}
Upload protein models from your computer:
~~~
$ scp 6N4O_SWISS_PROT_model_chainA.pdb \
6N4O_i-TASSER_model_chainA.pdb\
someuser@graham.computecanada.ca:scratch/workshop
~~~
{: .bash}

#### 1.2. Aligning protein models.
i-TASSER procedure changes the orientation of the protein and slightly optimizes the positions of all atoms. We will keep the original atom positions and take only the terminal end from the i-TASSER model. To combine the i-TASSER model with the actual 6n4o coordinates, we need to align these two structures.

Navigate to the working directory that you created on graham. Ensure that you have two files in the working directory:
~~~
[svassili@gra-login3 workshop]$ cd  ~/scratch/workshop
[svassili@gra-login3 workshop]$ ls
~~~
{: .bash}

~~~
6N4O_i-TASSER_model_chainA.pdb	6N4O_SWISS_PROT_model_chainA.pdb
~~~
{: .output}

 Launch VMD and in VMD prompt execute the following commands:
~~~
vmd > mol new 6N4O_SWISS_PROT_model_chainA.pdb
vmd > mol new 6N4O_i-TASSER_model_chainA.pdb
vmd > set 6n4o_residues "22 to 120 126 to 185 190 to 246 251 to 272 276 to 295 303 to 819 838 to 858"
vmd > set swissmodel [atomselect 0 "backbone and resid $6n4o_residues"]
vmd > set itasser [atomselect 1 "backbone and resid $6n4o_residues"]
vmd > set RotMat [measure fit $itasser $swissmodel]
vmd > echo rmsd before fit = [measure rmsd $itasser $swissmodel]
vmd > set itasser_all [atomselect 1 "all"]
vmd > $itasser_all move $RotMat
vmd > echo rmsd after fit = [measure rmsd $itasser $swissmodel]
vmd > set terminal [atomselect 1 "noh resid 1 to 21"]
vmd > $terminal writepdb 6n4o_resid_1-21.pdb
vmd > quit
~~~
{: .bash}
These commands will align the i-TASSER model  with the SWISS-MODEL. Combine the i-TASSER model of residues 1-21 and the SWISS-MODEL.
~~~
grep -h ATOM 6n4o_resid_1-21.pdb 6N4O_SWISS_PROT_model_chainA.pdb > 6n4o_chain_A_complete.pdb
~~~
{: .bash}

#### 1.3. Mutating residues
PDB entry 6N4O is the structure of the catalytically inactive hAgo2 mutant D669A. To construct the active form, we need to revert this mutation.

To accomplish this, we need to delete from ALA669 all atoms that are not present in ASP. Then change the residue name of ALA669 to ASP. Let's  first check what atoms are in residue 669:
~~~
$ grep 'A 669' 6n4o_chain_A_complete.pdb
~~~
{: .bash}

If you are familiar with aminoacid structures, you remember that the alanine sidechain is made of only one beta carbon atom (CB). All amino acids except glycine have beta carbon as well. So there is nothing to delete. All we need to do is to change the resName of all five ALA669 atoms to ASP.
You can do it using stream editor:
~~~
$ sed 's/ALA A 669/ASP A 669/g' 6n4o_chain_A_complete.pdb > 6n4o_chain_A_complete_A669D.pdb
~~~
{: .bash}

Verify the result:
~~~
$ grep 'A 669' 6n4o_chain_A_complete_A669D.pdb
~~~
{: .bash}

#### 1.4. Adding functionally important Mg2+ ion.
The catalytic site of hAgo2 is comprised of the three acidic amino acids D597, E637, and D669. It is known that hAgo2 requires a divalent metal ion near its catalytic site to slice mRNA. The 6N4O PDB file does not have this ion, but another hAgo2 structure, 4W5O, does. We can align these two structures as we did in section 3 and then copy the Mg2+ ion located near the catalytic site from 4W5O to our model.

#### 1.5. Assigning Protonation States to Residues
Use H++ server to calculate pKa of titratable sites and select protonation states as described in Episode 6 - "Assigning Protonation States to Residues in a Protein".

### 2. Adding missing segments to RNA structure files.

First, we need to create a PDB file containing all RNA atoms placed at proper positions. At this initial step, we are not particularly concerned with the quality of the 3D structure because we will refine it afterward.

We can insert the missing residues using the freely available [ModeRNA server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone ModeRNA software. The automatic process used by ModeRNA server moves residues adjacent to the inserted fragment. Besides, modeRNA server offers a limited set of options, and the output PDB files will need more processing steps afterward. Changing atomic positions is not desirable because we want to keep all experimental coordinates. For these reasons, we will use the standalone modeRNA package. You can install modeRNA on CC systems, or if you are comfortable with installation of Python, you can install it on your computer.

#### 2.1. Installing ModeRNA on your computer

1. Install python/2.7.14, numpy/1.11.3, biopython/1.58
2. [Download ModerRNA](http://genesilico.pl/moderna/download/) and follow [installation instructions](http://genesilico.pl/moderna/installing/).

#### 2.3. Installing ModeRNA on CC systems.
~~~
$ module load StdEnv/2016.4 python/2.7.14
$ virtualenv e27
$ source e27/bin/activate
$ pip install numpy==1.11.3 biopython==1.58 ModeRNA==1.7.1
~~~
{: .bash}

As ModeRNA can only model a single RNA strand, we will model chains C and D separately. For this job, we will need to prepare several files.

#### 2.4. Preparing structural templates for chains C and D.

Download 6n4o.pdb
~~~
$ wget https://files.rcsb.org/download/6n4o.pdb
~~~
{: .bash}
Residue 6 in chain D has only phosphate atoms and thus can not be used as a template. To prevent modeRNA from using it, we need to delete residue 6. Leaving it in the PDB file will lead to unwanted results.

Execute the following VMD commands to save chain C and chain D without residue 6:
~~~
vmd> mol new 6n4o.pdb
vmd> set sel [atomselect top "chain C or (chain D and not resid 6)"]
vmd> $sel writepdb 6n4o_chains_CD.pdb
~~~
{: .bash}

We created the file 6n4o_chains_CD.pdb suitable for use as a structural template.


#### 2.5. Preparing sequence alignment files for chains C and D.
Prepare two sequence alignment files for chains C and D. Each file should contain two sequences, the sequence of the model to be built and the template sequence.

Sequence alignment file for chain C, 6n4o_C.fasta:
~~~
>Model
UGGAGUGUGACAAUGGUGUUU
>Template
UGGAGUGUG-CAAUGGUG-UU
~~~
{: .source}

Sequence alignment file for chain D, 6n4o_D.fasta:
~~~
>Model
CAUUGUCACACUCCAAAA
>Template
CAUUG---CACUCCAA--
~~~
{: .source}

#### 2.6. Inserting missing segments.
Once you install modeRNA program, you will be able to use all functions. Below are commands needed to build chains C and D.  Description of all commands is available [here](http://genesilico.pl/moderna/commands/).
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
{: .source}

Save these commands in the file make_models.py and run the program:
~~~
$ module load StdEnv/2016.4 python/2.7.14
$ source e27/bin/activate
$ python make_models.py
~~~
{: .source}
The modeRNA program will run and create two model files, chain_C_model_A.pdb, and chain_D_model_B.pdb. Ensure that the program added all missing residues and did not move any residues present in the original structure file.

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
$ wget https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz --no-check-certificate
$ tar -xf SimRNA_64bitIntel_Linux.tgz
~~~
{: .bash}

#### 3.2. Preparing input files for simulation with SimRNAweb server.
If you will be using standalone SimRNA program you can skip this sections and proceed to the next one.

SimRNAweb server requires the user to provide RNA sequence, a list of residues not allowed to move in simulation and RNA structure in PDB format.

Sequence:
~~~
UGGAGUGUGACAAUGGUGUUU CCAUUGUCACACUCCAAA
~~~
{: .bash}
The convention is that the first and the second sequences correspond to chains A and B in the PDB structure file.

List of residues not allowed to move (we don't want the program to move atoms resolved in the experimental structure).
~~~
 A:1-9,11-18,20-21;B:1-5,9-16
~~~
{: .bash}
PDB file matching the sequence. All atoms must be present, and chains must be named A and B. To prepare this file combine chains A and B:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

You can use chains_CD_model_AB.pdb for SimRNAweb simulation. You will need to modify the PDB structure file as described in the next section for simulation with standalone SimRNA program.

#### 3.3. Preparing structure file for simulation with standalone SimRNA program.

Command-line SimRNA program does not need sequence and list of frozen atoms. You need to incorporate this information into the PDB structure file.

Begin with combining chains A and B if you have not done this yet:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

Next apply two modifications to this file. First, you need to add phosphate to the 5' terminal residue of chain D. SimRNA expects all residues to have a P atom. SimRNAweb will add P automatically, but for simulation with a standalone SimRNA program, you need to do it manually. There are several options to add phosphate.

##### 3.3.1. Renaming O5' atom to P
The most straightforward fix is to rename O5' atom to P. if you chose to do this, save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb, and skip the next step.

##### 3.3.2. Adding 5' monophosphate with AmberTools/20.
First, rename the phosphorylated 5' terminal nucleotide according to AMBER convention. The names of phosphorylated terminals in AMBER are A5, C5, G5, U5, DA5, DC5, DG5, DT5. Libraries of phosphorylated 5' terminal nucleotides are in the file 'terminal_monophosphate.lib'.

Launch Leap and load RNA force field:
~~~

$ tleap -f leaprc.RNA.OL3 -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}
In Leap promt execute the commands:
~~~
> loadoff terminal_monophosphate.lib
> chainD = loadpdb chain_D_model_B.pdb
> savepdb chainD chain_D5P.pdb
~~~
{: .bash}
These commands will load libraries of phosphorylated 5' terminal nucleotides and chain D PDB file. Leap will automatically add all missing atoms based on library entries and save chainD in PDB file chain_D5P.pdb:

We don't want to use the PDB file prepared with Leap for SimRNA because AMBER has different aminoacid naming conventions. So we copy phosphate atoms from chain_D5P.pdb and paste them into chains_CD_model_AB.pdb. We then edit chain ID, residue ID, and residue name. Save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb

##### 3.3.3. Adding 5' monophosphate with [CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader).
You can also add phosphate using CHARMM-GUI. Beware that CHARMM-GUI changes residue names to the old-style RNA 3-letter names and changes chain ID to "R".

##### 3.3.4 Define frozen atoms.
Standalone SimRNA program accepts PDB file where frozen atoms have occupancy 0.0 and completely free have occupancy 1.0. You can change values of the occupancy with the following VMD commands:
~~~
vmd> mol new chains_CD_model_AB_5P.pdb
vmd> set sel [atomselect top all]
vmd> $sel set occupancy 0
vmd> set sel [atomselect top "chain A and resid 10 19"]
vmd> $sel set occupancy 1
vmd> set sel [atomselect top "chain B and resid 6 7 8 17 18"]
vmd> $sel set occupancy 1
vmd> set sel [atomselect top all]
vmd> $sel writepdb chains_CD_model_AB_5P_frozen.pdb
~~~
{: .bash}

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
{: .source}

In the working directory, make a symbolic link to the 'data' directory located in SimRNA distribution. Assuming that you installed SimRNA in $HOME the link command is:
~~~
$ ln -s data ~/SimRNA_64bitIntel_Linux/data
~~~
{: .bash}
Then run the simulation:
~~~
$ srun -A def-someuser -c10 --mem-per-cpu=1000 --time=30:0 \
~/SimRNA_64bitIntel_Linux/SimRNA \
-P chains_CD_model_AB_5P_frozen.pdb \
-c config -E 10
~~~
{: .bash}

Option -E \<number of replicas> turns on replica exchange mode.
Replica exchange mode is parallelized with OMP.

The simulation will run for about two minutes and produce trajectory file *.trafl for each replica.

#### 3.5. Processing simulation trajectory
The simplest way of processing the trajectory files is obtaining the lowest energy structure. Generally, better results can be obtained by using clustering. Clustering tool is included with the distribution, but using clustering is outside the scope of this workshop.

Extract the lowest energy frame from the trajectory of the first replica
~~~
$ ~/SimRNA_64bitIntel_Linux/trafl_extract_lowestE_frame.py \
chains_CD_model_AB_5P_frozen.pdb_01.trafl
~~~
{: .bash}
Convert the lowest energy frame to PDB format
~~~
$ ~/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs chains_CD_model_AB_5P.pdb \
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


### 4. Preparing simulation system
Launch Leap and load protein and RNA forcefields:
~~~
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
source $EBROOTAMBERTOOLS/amber.sh
tleap -f leaprc.RNA.OL3 -f leaprc.protein.ff14SB
~~~
{: .bash}
Load protein and RNA. Then combine them inot one unit.
~~~
rna = loadpdb chains_CD_minimized.pdb
prot = loadpdb 6n4o_chain_A_complete_A669D.pdb
sys = combine {prot,rna}
~~~
{: .bash}
After this follow Episode 7 "Solvating a System, Adding Ions and Generating Input Files".
