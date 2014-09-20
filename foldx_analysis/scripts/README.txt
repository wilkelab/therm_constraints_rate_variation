This directory contains all of the python scripts used to generate the ddG values with Foldx.

Pipeline:
In order to perform the analysis start to finish, run the scripts in the following order:

1. download_structures.py 
	This will download all of the pdbs from the Huang dataset.
	After running this script, there should be a directory called pdbs with all of the structures
	In addition there should be another directory called relevant pdb files, with only the necessary chains

2. create_repair_list.py
	This is a simple script that just creates a list of all of the pdb names. 
	This pdb_list.txt is one of the inputs needed to repair the structures in Foldx.
	
3. create_runfiles.py
	***This is run after you repair all of pdbs.***
	This is a script that creates all of the pdb lists and runfiles needed to run Position Scan 
	in Foldx for each protein in the dataset. It also creates a paramlist that was used running
	FoldX on a HPC.
	
4. combine_raw_files.py
	***This is run after you have run Position Scan on all of the pdbs using the runfiles generated with the above script.***
	In order to parallelize the Position Scan run, each sequence of each protein was split among several Position
	Scan runs. This means that sequence of each protein was split up and run using a separate Position Scan run.
	This script copies all of the files for a give protein from its separate ouput files and places them into
	one directory. At the end of running this script, each pdb should have one folder containing all of the 
	"energies" files for each residue
	
5. extract_foldx_data.py
	This file creates a file for each protein with the ddG info for each protein. Each row is a residue in the 
	pdb and each column is the ddG that occurs when you mutant the wildtype amino acid to another amino acid.
	The wildtype to itself was represented as zero.
	
	
The file ddG_var_helper.py is just a helper file used by some of the others. It just needs to be in the directory.