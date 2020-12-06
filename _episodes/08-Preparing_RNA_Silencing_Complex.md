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

For this workshop we have chosen a complex of human argonaute-2 (hAgo2) protein with a micro RNA (miRNA) bound to a target messenger RNA (mRNA).

miRNAs are short noncoding RNAs that are critical for the regulation of gene expression and the defence against viruses. miRNAs guide argonaute proteins to mRNAs targeted for repression.

miRNAs regulate a wide variety of human genes. It is estimated that they regulate over 50% of all genes. So there is enormous number of different miRNA types. miRNAS play essential roles in both normal physiological functions and progression of disease such as cancer.

So miRNAs can control production of proteins. They do it by targeting mRNAs (RNAs that are used as blueprints for synthesis of all proteins). miRNAS can very specifically control expression of individual proteins and their selectivity is based on sequence complementarity between miRNAs and mRNAs.

miRNAs which target mRNAs encoding oncoproteins can serve as very selective tumour suppressors. They can inhibit tumor cells without negative impact on all other types of cells in a body. Discovery of this function of miRNAs have made miRNAs attractive tools for new therapeutic approaches. However, it is challenging to identify the most efficient miRNAs that can be targeted for therapeutic purposes. That's why understanding of the structural basis of the molecular mechanism by which miRNAs control production of proteins is very important.

To regulate protein synthesis a miRNA has to be loaded into an Ago2 protein forming the RNA-induced silencing complex that recognizes and inhibits the target messenger RNA  with high sequence specificity.

Therefore elucidating the structural basis of the molecular recognition between hAgo2 and mRNA is crucial for both the in-depth understanding of miRNA functions and the development of new therapeutics for human diseases.

### Adding missing residues to protein structure files.

Almost all protein and nucleic acid crystallographic structure files have missing residues. The reason for it is that the most flexible parts of biopolymers are disordered in crystals and therefore positions of their atoms cannot be accurately determined. These atoms, however may be crucial for MD simulations (e.g. loops connecting functional domains, nucleic acid chains, incomplete aminoacid side chains ... etc). For realistic simulation we need to build a model contating all atoms.

#### Adding missing residues using SWISS MODEL

1. Download structure and sequence files from PDB database:

~~~
$ wget https://files.rcsb.org/download/6n4o.pdb
$ wget https://www.rcsb.org/fasta/entry/6N4O/download -O 6n4o.fasta
~~~
{: .bash}

Extract full sequence of chain A:

~~~
$ grep -A1 "Chain A" 6n4o.fasta > 6n4o_chain_A.fasta
~~~
{: .bash}

Extract chain A atoms from 6n4o.pdb using VMD.
Load VMD module and launch VMD:
~~~
$ module load vmd
$ vmd
~~~
{: .bash}

In VMD prompt execute the following commands:

~~~
vmd > mol new 6n4o.pdb
vmd > set sel [atomselect top "chain A"]
vmd > $sel writepdb 6n4o_chain_A.pdb
vmd > quit
~~~
{: .bash}

2. Navigate to [SWISS MODEL](https://swissmodel.expasy.org)
3. Click Start Modelling,
4. Click User Template,
5. Paste a full sequence of your protein or upload 6n4o_chain_A.fasta,
6. Upload your structure file 6n4o_chain_A.pdb missing residues,
7. Click Build Model.
8. Download the homology model, and check it with your original structure.

Were all missing residues added?
- SWISS MODEL does not add terminal aminoacids

#### Other homology modeling servers
Another homology modeling server [i_TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/) (Iterative Threading ASSEmbly Refinement uses advanved protocol and is capable for threading terminal fragments. The downside is that i_TASSER process is much longer (about 60 hours for protein like 6n4o) and and positions of all atoms will be optimized. The resuls of i_TASSER model is in file 6N4O_i-TASSER_model_chainA.pdb.

### Aligning protein models.
i-TASSER procedure changes orientation of the protein and slightly optimizes positions of all atoms. We will keep the original atom positions and take only the terminal end from the i-TASSER model. To be able to combine the i-TASSER model with the original 6n4o coordinates we need to align these two structures.

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

Now we can combine i-TASSER and SWISS models.
~~~
grep -h ATOM 6n4o_resid_1-21.pdb 6N4O_SWISS_PROT_model_chainA.pdb > 6n4o_chain_A_complete.pdb
~~~
{: .bash}

### Mutating residues
PDB entry 6N4O is structure of the catalytically inactive hAgo2 mutant D669A. To construct the active form we need to revert this mutation.

To do this we need to delete from alanine 669 all atoms that are not present in aspartate. Then we change residue name of alanine 669 to ASP. Before we do this let's check what atoms are in residue 669:
~~~
$ grep 'A 669' 6n4o_chain_A_complete.pdb
~~~
{: .bash}

If you are familiar with aminoacid structures you remember that alanine sidechain is made of only one beta carbon atom (CB). All aminoacids except glycine have beta carbon as well. So there is nothing to delete. All we need to do is to change resName of all five alanine 669 atoms to ASP.
This can be done using stream editor:
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
Catalytic site of Ago2 is located in PIWI domain and is comprized of the acidic residues D597, E637 and D669. It is known that Ago2 requires divalent metal ion near its catalytic site to inactivate RNA. 6N4O PDB file does not have this ion, but another Ago2 structure, 4W5O does. We can align these two files and them copy Mg2+ ion from 4W5O to out model.


### Adding missing residues to RNA structure files.

We first need to create a PDB file containing all RNA atoms placed at a reasonable positions. At this initial step we are not particularly concerned with the quality of the 3D structure because we will refine it afterwards.

Insertion of the missing residues could be done using freely available [ModeRNA server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone ModeRNA software. However, the automatic process used by ModeRNA server moves residues adjacent to the inserted fragment. In addition modeRNA server offers a limited set of options and the output PDB files will need more processing steps afterwards. This is not desirable in our case (we want to keep all experimental positions untouched). For these reasons we will use the standalone modeRNA package.

You can install ModeRNA on your own computer or use it on CC systems.

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

ModeRNA can only model a single strand of RNA, so we will model chains C and D separately. For this job we will need the following files:

#### Preparing structural templates for chains C and D.

Download 6n4o.pdb:
~~~
$ wget https://files.rcsb.org/download/6n4o.pdb
~~~
{: .bash}

Residue 6 in chain D has only phosphate and thus can not be used as a template. To prevent modeRNA from using it we need to delete residue 6. Leaving it in PDB file will lead to unwanted results.

Execute the following VMD commands to save chain C and chain D without residue 6:
~~~
mol new 6n4o.pdb
set sel [atomselect top "chain C or (chain D and not resid 6)"]
$sel writepdb 6n4o_chains_CD.pdb
~~~
{: .bash}

We created the file 6n4o_chains_CD.pdb suitable for use as a structural template.


#### Sequence alignment files for chains C and D.
Prepare two sequence alignment files for chains C and D. Each file should contain two sequences, the sequence of the model to be built and the template sequence.

FASTA file for chain C, 6n4o_C.fasta:
~~~
>Model
UGGAGUGUGACAAUGGUGUUU
>Template
UGGAGUGUG-CAAUGGUG-UU
~~~
{: .source}

FASTA file for chain D, 6n4o_D.fasta:
~~~
>Model
CAUUGUCACACUCCAAAA
>Template
CAUUG---CACUCCAA--
~~~
{: .source}

#### Python script to make ModeRNA models, make_models.py:
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

The modeRNA program will run and create two model files, chain_C_model_A.pdb and chain_D_model_B.pdb. examine them and make sure that all missing residues have been added and all residues present in the original pdb file have not been moved.

> ## How will ModeRNA server do the same task?
> Try using ModeRNA server. Compare automatically generated ModeRNA models with the original chains C and D. Did server only added missing residues without moving any other atoms?
>
{: .challenge}

#### Pitfalls with structure generated automatically by ModeRNA.

ModeRNA webserver inserts missing residies with the same residue sequence number (ResSeq) as the ResSeq of the residue preceding gap in the template. It uses insertion code (iCode) to differentiate inserted res. For example if residues 6-8 are missing they all will be assigned ResSeq 5 and iCode A, B, C. It is not possible to renumber residues in ModeRNA web server automatically, and it is not possible to change model chainID (it is always A). So you will need to take care of residue numbering and renaming chain D manually.

**References:**
1. [ModeRNA: a tool for comparative modeling of RNA 3D structure](https://doi.org/10.1093/nar/gkq1320)
2. [PDB file format](
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)

### Refining double stranded RNA.
As RNA strands were modeled independently of each other, the side chains of the complementary strands may clash. To obtain a good duplex structure we need to minimize energy of both chains together.

Refining of the double stranded RNA can be done using freely available [SimRNAweb server](http://iimcb.genesilico.pl/modernaserver/submit/model/) or standalone  [SimRNA software] (https://genesilico.pl/SimRNAweb). This program allows for RNA 3D structure modeling with optional restraints.

SimRNA features:
- Coarse-grained representation (5 atoms per residue)
- Monte Carlo method for sampling the conformational space
- The energy function is composed of statistical potential terms, derived from the observed frequencies of occurrence of various proximate structural patterns.

#### Installation of the SimRNA binary package:
Simulation on public SimRNAweb server may wait up to a few days in the queue, while it can be done on a local compouter in a couple of minutes. SimRNA is available as a binary distribution, so no installation is required. You only need to download and unpack the package:
~~~
$ wget https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz --no-check-certificate
$ tar -xf SimRNA_64bitIntel_Linux.tgz
~~~
{: .bash}

#### Input required for refinement of the missing residues in double stranded RNA.

SimRNAweb server requires user to submit RNA sequence, RNA structure in PDB format and a list of residues to be frozen in simulation. Command line SimRNA program requires only PDB structure file specially prepared to flag frozen atoms in the occupancy field.

1. Sequence:
```
UGGAGUGUGACAAUGGUGUUU CCAUUGUCACACUCCAAA
```
The convertion is that the first sequence corresponds to chain A and the second to chain B in the PDB strucure file.

2. List of residues not allowed to move (we don't want the program to move atoms resolved in the experimental structure).
```
 A:1-9,11-18,20-21;B:1-5,9-16
```
3. PDB file matching the sequence. All atoms must be present and chains must be named A (matching the first sequence) and B (matching the second). To prepare this file we need to rename chain A in the model of chain D to chain B and then combine models of chains C and D into one file:

Combine chains A and B:
~~~
cat chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

#### Preparing pdb file for simulation with standalone SimRNA program.
SimRNA expects all residues to have P atom. SimRNA web will do it automatically, but for standalone SimRNA program 5' terminal of chain D must be phosphorylated manually. There are several options to add phosphate.

The simplest fix is to rename O5' atom to P. You can do this and skip the next step.

**Adding 5' monophosphate with AmberTools/20.**

First the phosphorylated 5' terminal nucleotide should be renamed according to AMBER convention. Phosporylated terminals in AMBER have names A5,C5,G5,U5,DA5,DC5,DG5,DT5.
AmberTools/20 libraries of phosphorylated 5' terminal nucleotides are in the file 'terminal_monophosphate.lib', so we nned to load RNA force field and this file.
~~~
$ module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 ambertools/20
tleap -f leaprc.RNA.OL3 -I $EBROOTAMBERTOOLS/dat/leap/lib/
~~~
{: .bash}
The in Leap execute the commands:
~~~
loadoff terminal_monophosphate.lib
chainD = loadpdb chain_D_model_B.pdb
savepdb chainD chain_D5P.pdb
~~~
{: .bash}

We don't want to use PDB file prepared with leap for SimRNA because AMBER has different aminoacid naming conventions. So we just copy phosphate atoms from chain_D5P.pdb, paste them into chain_D_model_B.pdb and edit chain and residue IDs.

Now we can combine chain_C_model_A.pdb and chain_D_model_B.pdb files. VMD adds the title line to each saved pdb file. This line must be removed. This can be done with the command:
~~~
$ grep -vh CRYST1 chain_C_model_A.pdb chain_D_model_B.pdb > chains_CD_model_AB.pdb
~~~
{: .bash}

**Adding 5' monophosphate with [CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader).**
Charmm gui changes residue names to 3-letter code and changes chain ID to "R".

**Changing occupancy values in PDB files with VMD**
Standalone SimRNA program accepts PDB file where frozen atoms have occupancy 0 and completely free have occupancy 1 as an input.

~~~
mol new chains_CD_model_AB.pdb
set sel [atomselect top all]
$sel set occupancy 0
set sel [atomselect top "chain A and resid 10 19"]
$sel set occupancy 1
set sel [atomselect top "chain B and resid 6 7 8 17 18"]
$sel set occupancy 1
set sel [atomselect top all]
$sel writepdb chains_CD_model_AB_frozen.pdb
~~~
{: .bash}

#### Running simulation

You should have two file is the working directory:
chains_CD_model_AB_frozen.pdb and the file with SimRNA configuration settings, config.

Contents of the file 'config':
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

In the working directory make a symbolic link do the 'data' directory in SimRNA distribution. Assuming that you installed SimRNA in $HOME:
~~~
$ ln -s data ~/SimRNA_64bitIntel_Linux/data
~~~
{: .bash}
Then run the simulation:
~~~
$ ~/SimRNA_64bitIntel_Linux/SimRNA -P chains_CD_simRNA.pdb -c config -E 10
~~~
{: .bash}

-E \<number of replicas>. Turns on replica exchange mode.
Replica exchange mode is parallelized with OMP.
This command will run for about a minute and produce trajectory for each replica.


**Extracting a structure from the simulation trajectory:**

Extract the lowest energy frame from the trajectory of the first replica
~~~
$ ~/SimRNA_64bitIntel_Linux/trafl_extract_lowestE_frame.py chains_CD_simRNA.pdb_01.trafl
~~~
{: .bash}
Convert the lowest energy frame to PDB
~~~
$ ~/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs chains_CD_simRNA.pdb chains_CD_simRNA.pdb_01_minE.trafl 1 AA
~~~
{: .bash}

**References:**
1. [SimRNA: a coarse-grained method for RNA folding simulations and 3D structure prediction](https://doi.org/10.1093/nar/gkv1479)
2. [SimRNA manual](https://ftp.users.genesilico.pl/software/simrna/version_3.20/SimRNA_UserManual_v3_20_20141002.pdf)
3. [VMD TCL commands](https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node121.html)
