import sys
import re
import os, subprocess
import numpy as np
from Bio import PDB
from Bio import AlignIO

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

def main(argv = sys.argv):
	headers = resdict.keys()
	subprocess.call("mkdir -p ../../foldx_ddG/", shell = True)
	chains = []
	pdbs = []
	pdb_file = open("../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Huang et al data file
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list
	
	#For testing
	#pdbs = ["132L", "1A2T"]
	#chains = ["A", "A"]
	
	for i in xrange(0, len(pdbs)):
		pdb_dir = "../../raw_foldx_output/" + pdbs[i] + "_" + chains[i]
		if(os.path.exists(pdb_dir)): #Check if that PDB has a directory full of foldx output
			outfile = "../../foldx_ddG/" + pdbs[i] + "_" + chains[i] + "_foldx_ddG.txt"	 
			out = open(outfile, "w") #If directory exists create an outfile#Open the output file for the ddG data
			out.write("SITE\tGLY\tALA\tLEU\tVAL\tILE\tPRO\tARG\tTHR\tSER\tCYS\tMET\tLYS\tGLU\tGLN\tASP\tASN\tTRP\tTYR\tPHE\tHIS\n")
			print "Directory: " + pdb_dir + " exists" 
			ref_residues = [] #Will hold residue name info (ex. Valine, Proline etc) 
			ref_res_nums = [] #Will hold the position numbers corresponding to each amino acid
			chain = chains[i]
			pdb = pdbs[i]
			try:
				#ppb = PPBuilder()	
				p=PDB.PDBParser(QUIET = True)  #Get the structure
				s=p.get_structure('X', "../tacc_files/repaired_pdbs/RepairPDB_" + pdb + "_" + chain + ".pdb")
				ref_struct = s[0]
				ref_chain = ref_struct[chain]
				for res in ref_chain: #Get the residue names and the residue positions for the protein
					ref_residues.append( res.resname )
					ref_res_nums.append(res.id[1])				
			except KeyError:
				print "Something is wrong with the mapped residues"
		
			res_counter = 0
			for res in ref_res_nums: #For position in the protein chain
				reference = ref_residues[res_counter]
				res_file = pdb_dir + "/energies_" + str(res) + "_RepairPDB_" + pdb + "_" + chain + ".txt" #Name of the file with the calculated dGs
				
				#These two lines dGs (in dG_array) and the corresponding mutants (in aminos)
				aminos = np.genfromtxt(res_file, delimiter = "\t", usecols = (0) , skip_header = 1, dtype=None)
				dG_array = np.genfromtxt(res_file, delimiter = "\t", usecols = (1),  dtype=None )
				aminos = aminos[:20] #Only need the mutants not the WTref dG
				dG_array = dG_array[:21]
				m = len(dG_array)
				ddGs = dG_array - (dG_array[0]*(np.ones(m))) #Calculate the change in dG
				ddGs = ddGs[1:]
				headers = []
								
				for aa in aminos: #Only need the three letter code
					headers.append(aa[0:3])
				
				#These lines format the ddG values and write them to a file
				if headers[0] == reference:
					ddG_value = '%.6f' % 0.0 #Format the ddG values to all have six decimal places
				else:
					ddG_value = '%.6f' % ddGs[0]
				ddG_string = str(res) + "\t" + ddG_value 								
				for k in xrange(1, len(headers)):
					if headers[k] == reference:
						ddG_value = '%.6f' % 0.0
					else:
						ddG_value = '%.6f' % ddGs[k]
					ddG_string = ddG_string +  "\t" + ddG_value
				if(res_counter < (len(ref_res_nums) - 1)):
					ddG_string = ddG_string + "\n"	
				out.write(ddG_string)
				res_counter = res_counter + 1
			out.close()		
			print "ddG Result File Created!!!"
	
if __name__ == "__main__":
	main()