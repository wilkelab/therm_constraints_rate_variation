This is a directory containing python scripts that are needed to perform the evolutionary rate analysis used for the paper.

Contents:

rate_helper.py - A helper file that contains functions needed for the other python scripts

make_alignments.py - A file that convertes the fasta alignments to phylip format to by used in RAxML

get_evol_rates.py - This is a python script that extracts the rates from Rate4Site, maps them to the sequence alignments, and outputs a readable format.