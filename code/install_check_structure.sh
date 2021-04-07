#!/bin/bash
module --force purge
module load StdEnv/2020 python
virtualenv ~/env-biobb
source ~/env-biobb/bin/activate
pip install biobb-structure-checking 

BIOBB_HOME=$HOME/env-biobb/lib/python3.8/site-packages
mkdir ~/bin
ln -s $BIOBB_HOME/biobb_structure_checking/check_structure.py ~/bin/check_structure 
chmod +x $BIOBB_HOME/biobb_structure_checking/check_structure.py
