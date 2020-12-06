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

### Adding missing residues to protein structure files.

Almost all protein and nucleic acid crystallographic structure files have missing residues. The reason for it is that the most flexible parts of biopolymers are disordered in crystals, and therefore, the positions cannot be accurately determined. These atoms, however, may be crucial for MD simulations (e.g., loops connecting functional domains, nucleic acid chains, incomplete amino acid side chains ... etc.). For realistic simulation, we need to build a model contacting all atoms.

#### Adding missing residues using SWISS MODEL

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

Navigate to [SWISS-MODEL](https://swissmodel.expasy.org)

Click Start Modelling,

Click User Template,

Paste the full sequence of your protein or upload 6n4o_chain_A.fasta,

Upload your structure file 6n4o_chain_A.pdb missing residues,

Click Build Model,

Download the homology model, and save it in the file '6N4O_SWISS_PROT_model_chainA.pdb'. Compare it with your original structure. Were all missing residues added?

In the following sections we will assume that the SWISS-MODEL is saved in the file 6N4O_SWISS_PROT_model_chainA.pdb.

#### Other homology modeling servers
SWISS-MODEL server does not add terminal fragments. Another homology modeling server [i-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/) (Iterative Threading ASSEmbly Refinement) uses the advanced protocol and is capable of threading terminal fragments. The downside of i-TASSER is that the process is much longer (about 60 hours for protein like 6n4o). i-TASSER optimizes, positions of all atoms, which is not always desirable.

The result of the i-TASSER modeling is in file 6N4O_i-TASSER_model_chainA.pdb.

### Aligning protein models.
i-TASSER procedure changes the orientation of the protein and slightly optimizes the positions of all atoms. We will keep the original atom positions and take only the terminal end from the i-TASSER model. To combine the i-TASSER model with the actual 6n4o coordinates, we need to align these two structures.

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
~~~
{: .bash}

Now we can combine the i-TASSER model of residues 1-21 and the SWISS-MODEL.
~~~
grep -h ATOM 6n4o_resid_1-21.pdb 6N4O_SWISS_PROT_model_chainA.pdb > 6n4o_chain_A_complete.pdb
~~~
{: .bash}

### Mutating residues
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

### Adding functionally important Mg2+ ion
The catalytic site of hAgo2 is comprised of the three acidic amino acids D597, E637, and D669. It is known that hAgo2 requires a divalent metal ion near its catalytic site to slice mRNA. The 6N4O PDB file does not have this ion, but another hAgo2 structure, 4W5O, does. We can align these two structures and then copy the Mg2+ ion located near the catalytic site from 4W5O to our model.

### Adding missing residues to RNA structure files.

First, we need to create a PDB file containing all RNA atoms placed at proper positions. At this initial step, we are not particularly concerned with the quality of the 3D structure because we will refine it afterward.

We can insert the missing residues using the freely available [ModeRNA server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone ModeRNA software. The automatic process used by ModeRNA server moves residues adjacent to the inserted fragment. Besides, modeRNA server offers a limited set of options, and the output PDB files will need more processing steps afterward. Changing atomic positions is not desirable because we want to keep all experimental coordinates. For these reasons, we will use the standalone modeRNA package.

#### Installing on your computer

1. Install python/2.7.14, numpy/1.11.3, biopython/1.58
2. [Download ModerRNA](http://genesilico.pl/moderna/download/) and follow [installation instructions](http://genesilico.pl/moderna/installing/).

#### Installing [ModeRNA](http://genesilico.pl/moderna/) on CC systems.

~~~
$ module load StdEnv/2016.4 python/2.7.14
$ virtualenv e27
$ source e27/bin/activate
$ pip install numpy==1.11.3 biopython==1.58 ModeRNA==1.7.1
~~~
{: .bash}

*\*Note: one of the tests failed on Graham - could be some weird problem in the environment of my account.*

As ModeRNA can only model a single RNA strand, we will model chains C and D separately. For this job, we will need to prepare several files.

#### Preparing structural templates for chains C and D.

Download 6n4o.pdb:
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


#### Sequence alignment files for chains C and D.
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

#### Python script for inserting missing residues.
File make_models.py:
~~~
from moderna import *
# Model chain C
tC = load_template('6n4o_chains_CD.pdb', 'C')
aC = load_alignment('6n4o_C.fasta')
mC = create_model(model_chain_name = 'A')
apply_alignment(tC, aC, mC)
apply_indel(mC, '9', '11', 'A')
apply_indel(mC, '18', '20', 'A')
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

#### Running ModeRNA
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

#### Pitfalls with structure generated automatically by ModeRNA.
ModeRNA webserver inserts missing residues with the same residue number as the residue preceding gap. It uses an insertion code (iCode) to differentiate inserted residues. For example, if residues 6-8 are missing, the program will insert them as residues 5A, 5B, 5C. It is not possible to instruct ModeRNA webserver to renumber residues, and it is not possible to change model chainID either. So you will need to take care of residue numbering and renaming chains manually.

**References:**
1. [ModeRNA: a tool for comparative modeling of RNA 3D structure](https://doi.org/10.1093/nar/gkq1320)
2. [PDB file format](
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)

### Refining double stranded RNA.
As modeRNA modeled RNA strands independently of each other, they may clash when two chains are combined. To obtain a good duplex structure, we need to minimize the energy of the whole RNA duplex.

Refining of the double stranded RNA can be done using freely available [SimRNAweb server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone  [SimRNA software] (https://genesilico.pl/SimRNAweb). This program allows for RNA 3D structure modeling with optional restraints.

SimRNA features:
- Coarse-grained representation (5 atoms per residue)
- Monte Carlo method for sampling the conformational space
- The energy function is composed of statistical potential terms, derived from the observed frequencies of occurrence of various proximate structural patterns.

#### Installation of the SimRNA binary package:
Simulation submitted to public SimRNAweb server may wait up to a few days in the queue, while on a local computer, you can do it in a couple of minutes. SimRNA is available as a binary distribution, so no installation is required. You only need to download and unpack the package:
~~~
$ wget https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz --no-check-certificate
$ tar -xf SimRNA_64bitIntel_Linux.tgz
~~~
{: .bash}

#### Preparing input files for simulation with SimRNAweb server.
If you will be using standalone SimRNA program you can skip this sections.

SimRNAweb server requires the user to submit RNA sequence, RNA structure in PDB format, and a list of residues to freeze in simulation.

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
PDB file matching the sequence. All atoms must be present, and chains must be named A and B.

Combine chains A and B:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

You can now use chains_CD_model_AB.pdb for SimRNAweb simulation. You will need to modify the PDB structure file as described in the next section for simulation with standalone SimRNA program.

#### Preparing structure file for simulation with standalone SimRNA program.

Command-line SimRNA program does not need sequence and list of frozen atoms. You must incorporate this information into the PDB structure file.

Combine chains A and B if you have not done this yet:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

You will need to apply two modifications. First, you need to add  phosphate to the 5' terminal of chain D. SimRNA expects all residues to have a P atom. SimRNAweb will add P automatically, but for simulation with a standalone SimRNA program, you need to do it manually. There are several options to add phosphate.

##### Rename O5' atom to P
The most straightforward fix is to rename O5' atom to P. if you chose to do this, save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb, and skip the next step.

##### Add 5' monophosphate with AmberTools/20.
First, rename the phosphorylated 5' terminal nucleotide according to AMBER convention. The names of phosphorylated terminals in AMBER are A5, C5, G5, U5, DA5, DC5, DG5, DT5. Libraries of phosphorylated 5' terminal nucleotides are in the file 'terminal_monophosphate.lib'. We need to load the RNA force field and this file.
~~~
$ module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
$ tleap -f leaprc.RNA.OL3 -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}
Then in Leap execute the commands:
~~~
> loadoff terminal_monophosphate.lib
> chainD = loadpdb chain_D_model_B.pdb
> savepdb chainD chain_D5P.pdb
~~~
{: .bash}
We don't want to use the PDB file prepared with Leap for SimRNA because AMBER has different aminoacid naming conventions. So we copy phosphate atoms from chain_D5P.pdb and paste them into chains_CD_model_AB.pdb. We then edit chain ID, residue ID, and residue name. Save the edited file chains_CD_model_AB.pdb as chains_CD_model_AB_5P.pdb

**Adding 5' monophosphate with [CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader).**
Charmm gui changes residue names to 3-letter code and changes chain ID to "R".

2. Flag frozen atoms in the occupancy field.

**Changing occupancy values in PDB files with VMD**
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


#### Running simulation

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

The simulation will run for about a minute and produce trajectory file *.trafl for each replica.


#### Extracting a structure from a simulation trajectory

Extract the lowest energy frame from the trajectory of the first replica
~~~
$ ~/SimRNA_64bitIntel_Linux/trafl_extract_lowestE_frame.py \ chains_CD_model_AB_5P_frozen.pdb_01.trafl
~~~
{: .bash}
Convert the lowest energy frame to PDB format
~~~
$ ~/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs chains_CD_model_AB_5P.pdb \ chains_CD_model_AB_5P_frozen.pdb_01_minE.trafl 1 AA
~~~
{: .bash}

#### References
1. [SimRNA: a coarse-grained method for RNA folding simulations and 3D structure prediction](https://doi.org/10.1093/nar/gkv1479)
2. [SimRNA manual](https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_UserManual_v3_20_20141002.pdf)
3. [VMD TCL commands](https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node121.html)
