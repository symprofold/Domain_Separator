# domain_separator

The tool takes a multidomain pdb file as input, identifies the structural domains and creates a fasta file with domain sequences separated by line breaks.

## Requirements
ChimeraX 1.6  
python modules: bs4, copy, glob, math, matplotlib, networkx, numpy, os, sys, traceback

## Running domain_separator
1.  copy pdb files into the folder input/
2.  start ChimeraX
3.  execute this command:
    run /path_to_project/Domain_Separator/domain_separator.py
    
## Output files
Output fasta files with domain sequences separated by line breaks are written into the folder input/.
