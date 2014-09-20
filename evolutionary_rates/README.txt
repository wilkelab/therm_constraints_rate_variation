This is a directory that contains the evolutionary rate data that was generated using rate4site.
The evolutionary rates were calculated using the JTT model.

evol_rates_JTT/
	Contains the raw output files from rate4site using the JC model	
	true_rates/
		Raw output files for with the "unnormalized rates"
	z_rates/
		Raw output files with the z-rates		

evol_rates_JC/
	Contains the raw output files from rate4site aaJC model
	true_rates/
		Raw output files for with the "unnormalized rates"
	z_rates/
		Raw output files with the z-rates
	
phylip_alignments/
	The alignments in phylip alignments. These sequences were aligned using the software MAFFT.

raxml_trees/
	Contains the evolutionary tree corresponding to each alignment. Generated using RAxML.

scripts/
	Contains the scripts used the generate the all_evolutionary_rates.csv file

all_evolutionary_rates.csv - Contains the evolutionary rates for all proteins in the dataset