import re, os, math, string, subprocess
import ddG_var_helper
from Bio import PDB

#Description: This script that downloads all of the pdbs

def main():
	process = subprocess.call("mkdir -p ../../pdbs", shell = True) #Creates the directory that pdbs will be downloaded to
	process = subprocess.call("mkdir -p ../../relevant_pdbs_chains", shell = True) #Creates the directory that pdbs will be downloaded to

	pdbs = []
	chains = []
	i = 0
	output = open("../nondownloaded_pdbs.txt" , "w") #Open a file to print to pdbs that couldn't download
	pdb_file = open("../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list

	pdbl = PDB.PDBList() #Make a PDB.List() object
	for pdb in pdbs:
		chain = chains[i]
		try: #Try to download the pdb
			pdbl.retrieve_pdb_file(pdb, pdir = '../../pdbs/.') #This line downloads the PDB
			rename_string = "mv " + "../../pdbs/pdb"+ pdb.lower() + ".ent " + "../../pdbs/" + pdb + ".pdb"
			process = subprocess.call(rename_string, shell=True) #Rename the file to the .pdb extension
		except IOError as e:
			output.write(pdb + "\n")
			print "I/O error({0}): {1}".format(e.errno, e.strerror)		
		except:
			print "Unexpected error: ", sys.exc_info()[0]
			output.write(pdb + "\n")	
			raise
		output.close()
		
		p=PDB.PDBParser()
		try: #Try to get the structure 
			s=p.get_structure('X', "../../pdbs/"+ pdb + '.pdb')
			extracted_name = pdb + "_" + chain + ".pdb"
			back_var_helper.extract(s, chain, -2000, 10000, extracted_name)
		except IOError as e: #Probably was originally downloaded (probably an obsolete structure)
			print "PDB was not downloaded!"			
		i = i + 1

	cp_string = "cp *.pdb  ../../relevant_pdbs_chains" #Copy the extracted chains to a new directory
	process = subprocess.call(cp_string, shell = True)

	rm_string = "rm *.pdb" 
	process = subprocess.call(rm_string, shell = True) #Delete the extracted chains in this directory


if __name__ == "__main__":
	main()