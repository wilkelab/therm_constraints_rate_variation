import sys, os, subprocess

#Description: This is a script that copies all of the FoldX output from the result directories in a common directory for each pdb.
  
def main(argv = sys.argv):	
	chains = []
	pdbs = []
	pdb_file = open("../../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Huang et al data file
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list

	#For each pdb
	for i in xrange(len(pdbs)):
		pdb = pdbs[i]
		chain = chains[i]
		raw_out_dir = "../raw_foldx_output/" + pdb + "_" + chain
		subprocess.call("mkdir -p " + raw_out_dir, shell = True) #Create a result directory 
		for j in xrange(1, 15):
			result_dir = "../FoldX_results/" + pdb + "_" + chain + "_" + str(j)
			if(os.path.exists(result_dir)):
				print "Result Directory: " + result_dir + " exists!"
				subprocess.call("cp " + result_dir +  "/energies_* " +  raw_out_dir, shell = True) #Copy results into the directory
				print "Copied files!"
	
if __name__ == "__main__":
	main()