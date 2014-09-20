#This is a script that creates the list used for repairing the pdbs before the FoldX analysis
def main():
	output = open("pdb_list.txt","w")
	pdbs = []
	chains = []
	i = 0
	pdb_file = open("../../Huang_et_al_data/TableS1.csv", "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines extract the info from the Huang et al data file
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list
	
	i = 0
	for i in xrange(0, len(pdbs)): # Writes the pdbs to a list so they can be repaired
		filename = pdbs[i] + "_" + chains[i] + ".pdb"
		if (i <len(pdbs) - 1):
			output.write("pdbs/" + filename + "\n")
		else:
			output.write("pdbs/" + filename)
		i = i+1
	output.close()
	
if __name__ == "__main__":
	main()