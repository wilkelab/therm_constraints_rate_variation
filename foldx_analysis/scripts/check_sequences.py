import sys
import re
import os, subprocess
from Bio import PDB
from Bio import AlignIO
from Bio.PDB.Polypeptide import PPBuilder 
import numpy as np
import ddG_var_helper

def main(sys = sys.argv):	
	chains = []
	pdbs = []
	pdb_file = open("../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list

	#res_lists = []
	#aa_list = ""
	huang_res_lengths =  np.genfromtxt("../../Huang_et_al_data/TableS1.csv", delimiter = ",", dtype = "float", skip_header = 1, usecols = 7)	
	for i in xrange(0, len(pdbs)):
		pdb = pdbs[i]
		chain = chains[i]
		print pdb		
		try:
			ppb = PPBuilder()	
			p=PDB.PDBParser(QUIET = True)
			s=p.get_structure('X', "../tacc_files/repaired_pdbs/RepairPDB_" + pdb + "_" + chain + ".pdb")
			pp = ppb.build_peptides(s)[0]
			seq = str(pp.get_sequence())	
		
			#Get the structures
			ref_struct = s[0]
			ref_chain = ref_struct[chain]
			ref_residues = []
			ref_res_nums = []
			for res in ref_chain:
				print res
				ref_residues.append( res.resname )
				ref_res_nums.append(res.id[1])				
		except KeyError:
			print "Something is wrong with the mapped residues"
			#continue
		seq_array = ddG_var_helper.convert_to_one_letter_code(ref_residues)

		#print pdb + " Length: " + str(len(seq_array))

		print "My Length Ref Lengths: " + str(len(ref_res_nums))
		print "Huang Number Residues: "  + str(huang_res_lengths[i])
		print ref_res_nums
		
	#huang_res_names = np.genfromtxt("../../Huang_et_al_data/TableS2.csv", dtype = "str", delimiter = ",", skip_header = 1, usecols = (0,1,2, 3))
	#huang_chains =  np.genfromtxt("../../Huang_et_al_data/TableS2.csv", dtype = "str", delimiter = ",", skip_header = 1, usecols = (0,1,2, 3))

	
	
	
	
	#huang_data = np.loadtxt("../../Huang_et_al_data/TableS1.csv", delimiter = ",", skip_header = 1, usecols = (0,2))
	#print huang_data
	
if __name__ == "__main__":
	main(sys.argv)
	