import sys
import re
import back_var_helper
import os, subprocess
from Bio import PDB
from Bio import AlignIO
from Bio.PDB.Polypeptide import PPBuilder 

#Create 
def create_residue_list(pdb, chain):
	print pdb, chain
	aa_list = ""
	try:
		ppb = PPBuilder()	
		p=PDB.PDBParser(QUIET = True)
		s=p.get_structure('X', "../../relevant_pdbs_chains/" + pdb + "_" + chain + ".pdb")
		pp = ppb.build_peptides(s)[0]
		seq = str(pp.get_sequence())
		print seq
		
		#Get the structures
		ref_struct = s[0]
		ref_chain = ref_struct[chain]
		ref_residues = []
		ref_res_nums = []
		for res in ref_chain:
			ref_residues.append( res.resname )
			#print res.resname #Prints the residue name
			#print res.id #Prints the residue id
			ref_res_nums.append(res.id[1])				
	except KeyError:
		print "Something is wrong with the mapped residues"
		#continue
	seq_array = back_var_helper.convert_to_one_letter_code(ref_residues)
	#print seq_array
	for i in xrange(0, len(seq_array)):
		aa_list = aa_list + ","
		res_num = ref_res_nums[i]
		aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
	aa_list = aa_list + ";"	
	print aa_list
	return aa_list

def create_runfile(pdb, chain):
	filename = "runfiles/runfile_" + pdb + "_" + chain + ".txt"
	listname = "pdb_list_" + pdb + "_" + chain + ".txt"
	
	out = open(filename,"w") 
	res_string = create_residue_list(pdb, chain)
	out.write("<TITLE>FOLDX_runscript;" + "\n")
	out.write("<JOBSTART>#;" + "\n")
	out.write("<PDBS>#;" + "\n")
	out.write("<BATCH>" + listname + ";" + "\n")
	out.write("<COMMANDS>FOLDX_commandfile;" + "\n")
	out.write("<PositionScan>#" + res_string + "\n")
	out.write("<END>#;" + "\n")
	out.write("<OPTIONS>FOLDX_optionfile;" + "\n")
	out.write("<Temperature>298;" + "\n")
	out.write("<R>#;" + "\n")
	out.write("<pH>7;" + "\n")
	out.write("<IonStrength>0.050;" + "\n")
	out.write("<water>-CRYSTAL;" + "\n")
	out.write("<metal>-CRYSTAL;" + "\n")
	out.write("<VdWDesign>2;" + "\n")
	out.write("<OutPDB>false;" + "\n")
	out.write("<pdb_hydrogens>false;" + "\n")
	out.write("<complex_with_DNA> true;" + "\n")
	out.write("<END>#;" + "\n")
	out.write("<JOBEND>#;" + "\n")
	out.write("<ENDFILE>#;")
	out.close()
	
def create_paramlist(pdbs, chains):
	out = open("paramlist_mutate_pdbs", "w")
	for i in xrange(len(pdbs)):
		pdb_file = "RepairPDB_" + pdbs[i] + "_" + chains[i] + ".pdb"
		filename = "runfile_" + pdbs[i] + "_" + chains[i] + ".txt"
		if (i < len(pdbs) - 1):
			out.write("./foldx_mutate.sh " + pdbs[i] + "_" + chains[i] + " " + pdb_file + " " + filename + "\n")
		else:
			out.write("./foldx_mutate.sh " + pdbs[i] + "_" + chains[i] + " " + pdb_file + " " + filename)
	out.close()

def main(argv = sys.argv):
	chains = []
	pdbs = []
	pdb_file = open("../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list

	subprocess.call("mkdir -p pdb_lists", shell = True)
	subprocess.call("mkdir -p runfiles", shell = True)
	for i in xrange(len(pdbs)):	
		print pdbs[i]
		filename = "pdb_lists/pdb_list_" + pdbs[i] + "_" + chains[i] + ".txt"
		#print filename
		out = open(filename, "w")
		out.write("RepairPDB_" + pdbs[i] + "_" + chains[i] + ".pdb")
		out.close()
		create_runfile(pdbs[i], chains[i])
		create_paramlist(pdbs, chains)
	
if __name__ == "__main__":
	main(sys.argv)