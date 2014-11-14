import rate_helper as rh
import subprocess, os

'''
Description: This is a script that creates the phylip alignments used in the evolutionary rate analysis
'''

def main():
	pdbs, chains = rh.get_pdb_info("../../../Huang_et_al_data/TableS1.csv")
	process = subprocess.call("mkdir -p ../phylip_alignments", shell = True)
	for i in xrange(0, len(pdbs)):
		pdb = pdbs[i]
		chain = chains[i]
		temp_file = "temp.fasta"
		out_file = "../phylip_alignments/" + pdb + "_" + chain + "_aligned.phy"
		alignment = "../../protein_alignments/mapped_mafft_alignments/" + pdb + "_" + chain + "_mafft.fasta"
		print pdb, chain
		seqs, headers = rh.get_sequences(alignment)
		seq1 = seqs.pop(0)
		head1 = headers.pop(0)
		rh.write_seq_to_file(seqs, headers, temp_file)
		rh.convert_alignment_format(temp_file, "fasta", "phylip-relaxed", out_file, all_caps = True)
		os.unlink(temp_file)
		

if __name__ == '__main__':
	main()