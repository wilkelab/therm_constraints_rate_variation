import rate_helper as rh
import numpy as np


def get_rates(pdb, chain, rate_type, z_norm):
	'''
	This is a function that extracts and maps the evolutionary rates to the pdb sequence
	
	Args:
		pdb: The pdb name
		map_file: The file that maps the pdb sequence to the alignment	
		rate_type: JTT or LG (the amino acid model using in the rate4site calculations)
	Returns:
		sites: The sites in the pdb that have been mapped to the alignment	
		aligned_rates: The rates mapped back to sites in the pdb	
	'''
		
	align_sites = []
	rates = []
	amino_acids = []
	pdbs = ["1COY", "1CZ1", "1E6E", "1H3I", "1VQ1"]
	if z_norm == True:
		rate_file = "../evol_rates_" + rate_type + "/"  + "z_rates/" + pdb + "_" + chain + "_" + rate_type + "_rates.txt"	
	else:
		rate_file = "../evol_rates_" + rate_type + "/" + "true_rates/" + pdb + "_" + chain + "_" + rate_type + "_rates_unscaled.txt"	
		
	map_name = "../../pdb_to_alignment_maps/" + pdb + "_" + chain + "_map.txt"		#Name of the file that maps the pdb sequence
	align_name = "../../mafft_alignments/" + pdb + "_" + chain + "_mafft.fasta"
	
	#Get the ref sequences (should be the pdb sequence)
	seqs, headers = rh.get_sequences(align_name)
	ref_index = headers.index(">" + pdb+chain)
	ref_seq = seqs[ref_index]	 
	ref_seq.strip()
	
	input = open(rate_file)
	data = input.readlines()
	data = data[11:]
	data.pop()
	data.pop()
	
	if(rate_type == "JTT" and pdb == "1L7Q"):
		
		for line in data: #Get the evolutionary rates 
			site = line[2:8].strip()
			amino = line[9:12].strip()
			#rate = line[12:20].strip()
			#rate = line[11:19].strip()
			rate = line[10:19].strip()
			align_sites.append(site)
			amino_acids.append(amino)
			rates.append(rate)
			
	elif(rate_type == "JC" and pdb in pdbs):
		for line in data: #Get the evolutionary rates 
			site = line[2:8].strip()
			amino = line[9:12].strip()
			rate = line[10:19].strip()
			align_sites.append(site)
			amino_acids.append(amino)
			rates.append(rate)	
	
	else:
		for line in data: #Get the evolutionary rates 
			site = line[2:8].strip()
			amino = line[9:12].strip()
			rate = line[12:20].strip()
			align_sites.append(site)
			amino_acids.append(amino)
			rates.append(rate)
	
	sites = extract_sites(map_name)
	aligned_rates = []
	counter = 0
	ref_seq_array = rh.get_sequence_array(ref_seq)

	
	#Align the calculated rates to the proper amnio acid in the pdb using the pdb_to_alignment_map created earlier (in match_frequencies.py)
	for j in xrange(0, len(ref_seq_array)):
		if(ref_seq[j] == '-'):
			#print str(j+1), ref_seq_array[j], sites[j], "NA"
			aligned_rates.append("NA")
		else:
			aligned_rates.append(rates[counter])
			counter = counter + 1
	return sites, aligned_rates

def extract_sites(map_file):
	'''
	This is a function extracts the site information from the mapped file
	
	Args:
		pdb: The pdb name
		map_file: The file that maps the pdb sequence to the alignment
	
	Returns:
		sites: The sites in the pdb that have been mapped to the alignment	
		
	'''
	sites = []
	site_data = np.genfromtxt(map_file, delimiter = "\t", dtype = "string", skip_header = 1, usecols = (1)) #Read in the data
	sites = []
	for s in site_data:
		sites.append(s.strip()) #String off the new lines and append the site to the sites list
	return sites

def main():
	pdbs, chains = rh.get_pdb_info("../../Huang_et_al_data/TableS1.csv") #Get the pdbs and their corresponding chains
	outfile = "../all_evolutionary_rates.csv"
	out = open(outfile, "w")
	out.write('"pdb","chain","site","zr4s_JTT","r4s_JTT,"zr4s_JC","r4s_JC"\n')
	
		
	for i in xrange(0, len(pdbs)):

		pdb = pdbs[i]
		chain = chains[i]
		
		#if (pdb == '1A4L' or pdb =='1AGM' or pdb == '1JS4' or pdb == '1L1L' or pdb == '1L7Q' or '1OG1'):
		#	continue		
		
		#if (pdb == '1A4L' or pdb =='1AGM' or pdb == '1JS4' or pdb == '1L1L' or pdb == '1L7Q' or pdb == '1OG1' or pdb == '1ONR' or pdb == '1VQ1'):
		#	continue
		
		#if (pdb == '1A65' or pdb == '1AGM' or pdb == '1CB8' or  pdb == '1CDG' or  pdb == '1COY' or pdb == '1CZ1' or  pdb == '1E6E' or  pdb == '1EH5' or pdb == '1EUG' or pdb == '1EUU' or pdb == '1EX1' or pdb == '1GOG' or pdb == '1H19' or pdb == '1H3I' or pdb == '1HPL' or pdb == '1OG1' or  pdb == '1OGO' or  pdb == '1PA9' or  pdb == '1VQ1'):
		#	continue
		
		print pdb
		
		sites, rates_JTT_z_norm = get_rates(pdb, chain, "JTT", z_norm = True)
		sites, rates_JTT = get_rates(pdb, chain, "JTT", z_norm = False)
		
		sites, rates_JC_z_norm = get_rates(pdb, chain, "JC", z_norm = True)
		sites, rates_JC = get_rates(pdb, chain, "JC", z_norm = False)
		
		#print sites
		#print rates_JTT
		#print rates_LG
		
		if (len(sites) != len(rates_JTT_z_norm) or len(sites) != len(rates_JTT) or len(rates_JTT) != len(rates_JTT_z_norm)):
			print "Sites not mapped correctly to Rates!"
			print "Length of Sites: ", len(sites)
			print "Length of True Rates: ", len(rates_JTT)
			print "Length of Z Normed Rates: ", len(rates_JTT_z_norm)
			break
		
		j = 0
		while ( j < len(sites)): #For each site
			if (sites[j] != "NA"): #If the site is NOT NA (Meaning the the site is mapped to a site in the pdb...
				out.write('"' + pdb + '"' + ',' + '"' + chain + '"' + ',' + str(sites[j]) + ',' + str(rates_JTT_z_norm[j]) + ',' + str(rates_JTT[j]) +  ',' + str(rates_JC_z_norm[j]) + ',' + str(rates_JC[j]) + "\n") #Right out the frequency information for that site
				j = j + 1
			else: #else do not it
				j = j + 1
		
	out.close()		
if __name__ == "__main__":
	main()