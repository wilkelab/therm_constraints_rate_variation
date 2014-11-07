This is a directory that contains the files that were used to calculate the ddG values using ddG Monomer. 

Contents:

calculate_ddGs_rosetta.sh - This is the shell script written to run the protocol on the HPC (we used TACC)

convert_to_cst_files.sh - Used for the minimization 

ddG_flags_files/
	Contains the Rosetta flag files used in all to calculate the ddG values using ddG Monomer. There is a flag file for each Rosetta ddG Monomer run. Each "run" corresponds to calculating all possible amino acid substitutions at one site.

minimize_flag_files/
	Contains the flag files for the minimization protocol prior to running the ddG calculation protocol.

pdb_lists/
	Contains the pdb_list.txt files used for ddG Monomer

renumbered_pdbs/
	Contains the five pdbs used. There had renumbered residues.

resfiles/
	Contains the resfiles used in ddG Monomer.

sp2_paper_talaris2013_scaled.wts - This scoring weights file used 
