import sys
import re
import os, subprocess
import rosetta_helper as rh
import numpy as np
from Bio import PDB
from Bio import AlignIO

aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] 
ordered_aa_list = ['G', 'A', 'L', 'V', 'I', 'P', 'R', 'T', 'S', 'C', 'M', 'K', 'E', 'Q', 'D', 'N', 'W', 'Y', 'F', 'H'] #List of amino acids

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

def main(argv = sys.argv):
	headers = resdict.keys()
	subprocess.call("mkdir -p ../../rosetta_ddG/", shell = True)
	pdbs = ["2ACY", "1LBA", "1LJL", "1BP2", "1PYL"]
	chains = ["A", "A", "A", "A", "A"] 
	pdb_lengths = [98, 146, 130, 123, 95]
	
	for i in xrange(0, len(pdbs)):
		ref_residues = [] #Will hold residue name info (ex. Valine, Proline etc) 
		ref_res_nums = [] #Will hold the position numbers corresponding to each amino acid
		chain = chains[i]
		pdb = pdbs[i]
		print "PDB: ", pdb
		outfile = "../../rosetta_ddG/" + pdbs[i] + "_" + chains[i] + "_rosetta_ddG.txt"	
		out = open(outfile, "w") #If directory exists create an outfile a#Open the output file for the ddG data
		out.write("SITE\tGLY\tALA\tLEU\tVAL\tILE\tPRO\tARG\tTHR\tSER\tCYS\tMET\tLYS\tGLU\tGLN\tASP\tASN\tTRP\tTYR\tPHE\tHIS\n")
		try:
			p=PDB.PDBParser(QUIET = True)  #Get the structure
			s=p.get_structure('X', "../../foldx_analysis/tacc_files/repaired_pdbs/RepairPDB_" + pdb + "_" + chain + ".pdb")
			ref_struct = s[0]
			ref_chain = ref_struct[chain]
			
			for res in ref_chain: #Get the residue names and the residue positions for the protein
				ref_residues.append( res.resname )
				ref_res_nums.append(res.id[1])				
		except KeyError:
			print "Something is wrong with the mapped residues"

		print "i: ", i
		for j in xrange(1, (pdb_lengths[i]+1)):
			res = ref_res_nums[j-1]
			#print "j, res_num: , Residue: ", j, res, ref_residues[j-1]
			ddG_dict = {}
			pdb_dir = "../../raw_rosetta_output/" + pdbs[i] + "_" + chains[i] + "_" + str(j)
			data_file = pdb_dir + "/ddg_predictions.out"
			if(os.path.exists(data_file)): #Check if that PDB has a directory full of rosetta output		
				print "Directory: " + pdb_dir + " exists and has data" 
	
				res_counter = 0
				reference = ref_residues[res_counter]
				resfile = pdb_dir + "/ddg_predictions.out"
				ddG_array = np.genfromtxt(resfile, usecols = (2) , skip_header = 1, dtype=None) #Get the ddGs from the rosetta file
				
				k = 0
				for aa in aa_list: #Create a dictionary that maps the amino acids to the ddG of mutation to that amino acid
					ddG_dict[aa] = ddG_array[k]	
					k = k + 1			
				ddG_string = str(res) + "\t" #Create a string to print to the file
				for aa in ordered_aa_list: 
					ddG_value = '%.6f' % ddG_dict[aa]
					ddG_string = ddG_string +  "\t" + ddG_value
				out.write(ddG_string + "\n")
				print "ddG Result File Created!!!"
			elif(os.path.exists(data_file)!= True and ref_residues[j-1] == 'CYS'): #Cysteines in disulfide bonds do not have ddG data (were not mutated)
				print "Cysteine Residue"	
				ddG_string = str(res) + "\t"
				for z in xrange(0, len(ordered_aa_list)):
					ddG_string = ddG_string +  "\t" + "NA"
				out.write(ddG_string + "\n")
			else:
				print "DATA FILE NOT THERE!!!"
			
if __name__ == "__main__":
	main()