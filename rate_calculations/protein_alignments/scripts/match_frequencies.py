#!/usr/bin/python
import cor_analysis_helper as ch
import subprocess, os

'''
Last Edited By: Eleisha Jackson (June 23, 2014)
Description: This script create a map called *_map.txt. This map maps a position in the alignment to each position in the pdb sequence.
It also creates a MAFFT alignment for each muscle alignment. This alignment contains all sequences from the muscle alignment along with the pdb sequences (the first sequence in the alignment).
Must have MAFFT installed to use this script.
'''

def create_pdb_alignment_map(seqs, pdb_res_nums, map_filename):
	'''
	Creates a file that maps a pdb seqeunce onto a given alignments. 
	Args: 
		seqs: A list of aligned sequences 
		pdb_res_nums: A list of residue numbers from the pdb file
		map_filename: The file name to write the map to
	
	Returns:
		Does not return anything but creates a file with a map
	'''
	
	alignment_length = len(seqs[0]) 
	aligned_pdb_seq = seqs[0]
	#Make map list
	map_outfile = open(map_filename, "w")
	map_outfile.write("align_pos\tpdb_pos\tpdb_aa\n")
	res_counter = 0
	for j in xrange(0, alignment_length):
		map_outfile.write(str(j+1) + "\t")
		if(aligned_pdb_seq[j]!= "-"):
			map_outfile.write(str(pdb_res_nums[res_counter]) + "\t" + aligned_pdb_seq[j])
			res_counter = res_counter + 1
		else:
			map_outfile.write("NA" + "\t" + "-")
		map_outfile.write("\n")

def main():
	pdbs, chains = ch.get_pdb_info("../../../Huang_et_al_data/TableS1.csv") #Get the pdbs and their corresponding chains
	for i in xrange(0, len(pdbs)):
		subprocess.call("mkdir -p  ../mapped_mafft_alignments", shell = True)
		pdb = pdbs[i]
		chain = chains[i]
		pdb_file = "../../../ddG_calculations/foldX/repaired_pdbs/" + "RepairPDB_" + pdb + "_" + chain + ".pdb"
		align_file = "../../../Huang_et_al_data/alignments/" + pdb + chain + ".muscle_fasta"
		align_outfile = "../mapped_mafft_alignments/" + pdb + "_" + chain + "_mafft.fasta"
		map_filename = "../pdb_to_alignment_maps/"+ pdb + "_" + chain + "_map.txt"
		
		#Open the pdb, Get the sequence from the pdb
		#Get the residue list (make a function)
		(pdb_seq, pdb_res_nums) = ch.get_info_from_pdb(pdb_file, chain)	
			
		#Get unaligned sequences from the alignment
		(all_seqs, all_headers) = ch.get_unaligned_seqs(align_file)
			
		#Then add this to the top of the file
		all_seqs.insert(0, pdb_seq)
		all_headers.insert(0, ">" + pdb + "_" + chain)
			
		#Then write all of the sequences to a file, Align the sequences in a file
		ch.align_seqs_mafft(all_seqs, all_headers, align_outfile)	
		
		#Extract the sequences from the aligned file	
		(seqs, headers) = ch.get_sequences(align_outfile)
			
		#Count the number of non-gap residues and then compare to the residue number length
		if (len(pdb_res_nums) != ( len(seqs[0]) - seqs[0].count("-"))):
			print "Number of residues extracted does not match number of residues in pdb seq"
			print "Number of residues: ",len(pdb_res_nums)
			print "Length of sequence: ", seqs[0].count("-") 
			break
		
		create_pdb_alignment_map(seqs, pdb_res_nums, map_filename)

if __name__ == "__main__":
	main()