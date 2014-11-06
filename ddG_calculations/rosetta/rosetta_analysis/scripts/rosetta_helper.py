#!/usr/bin/python
import re, os, math, string, subprocess
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
from Bio import PDB

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

_hydrogen=re.compile("[123 ]*H.*") 

class ChainSelector(object): 
	"""
	Adapted from the extract() function of Biopython 
	"""

	def __init__(self, chain_id, start, end, model_id=0): 
		self.chain_id=chain_id 
		self.start=start 
		self.end=end 
		self.model_id=0 

	def accept_model(self, model): 
		# model - only keep model 0 
		if model.get_id()==self.model_id: 
			return 1 
		return 0 

	def accept_chain(self, chain): 
		if chain.get_id()==self.chain_id: 
			return 1 
		return 0 
		
	def accept_residue(self, residue): 
		# residue - between start and end 
		hetatm_flag, resseq, icode=residue.get_id() 
		if hetatm_flag!=" ": 
			# skip HETATMS 
			return 0 
		if self.start<=resseq<=self.end: 
			return 1 
			return 0 

	def accept_atom(self, atom): 
		# atoms - get rid of hydrogens 
		name=atom.get_id() 
		if _hydrogen.match(name): 
			return 0 
		else: 
			return 1 

def extract(structure, chain_id, start, end, filename): 
	""" 
	Write out selected residues from a particular given pdb structure to filename. 
	
	Args: 
		structure: The four letter pdb code for the protein (this pdb must currently exist!)
		chain_id: The chain in the pdb you are trying to extract
		start: The beginning residue that you want to start your extraction from
		end:
		filename: The filename (ex. my_protein.pdb) of the pdb that you will extract the selected residues to  	
	Returns:
		A pdb file that contains only the selected residues from the original pdb structure given as an argument
		
	Example of usage
		p=PDBParser()
		s=p.get_structure('X', '1RUZ.pdb')
		extract(s, 'H', 40, 60, 'extracted.pdb')
 
	""" 
	sel=ChainSelector(chain_id, start, end) 
	io=PDBIO() 
	io.set_structure(structure) 
	io.save(filename, sel) 


def get_pdb_info(data_file):
	pdbs = []
	chains = []
	pdb_file = open(data_file, "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list		
	return pdbs, chains
	
def get_sequence_array(seq):
	"""
	Turns a string representing a string into a list of characters.
	
	Args:
		seq: A string representing a sequence (ex. nucleotide, protein)
		
	Returns:
		seq_array: A list of single letter elements representing the sequence (ex. each amino acid in a protein is an element in the list)
	
	"""

	seq_array = []
	seq_length = len(seq)
	count = 0
	while count < seq_length:
		aa = seq[count]
		seq_array.append(aa)
		count = count + 1
	return seq_array

def array_to_seq(seq_array):
	"""
	Takes an array that represents a sequence and returns a string representation
	
	Args:
		seq_array: An array where each element is an character in a sequence
	
	Returns:
		seq: A string representation of the sequence
		
	"""

	seq = "" #Create an empty string
	for i in xrange(0, len(seq_array)): #For each element in the seq array
		seq = seq + seq_array[i]	 #Append the character to the empty string
	return seq

def get_info_from_pdb(pdb_file, chain):
	"""
	Extracts the sequences from of Protein Database File (PDB) file
	
	Args:
		pdb_file = A string that is the path to the file (just the name if in the same directory)
		chain = The chain in the pdb for the sequence you want to extract
			
	Returns:
		seq: A string representation of the sequence in one letter code
			
	"""
	p=PDB.PDBParser()
	s=p.get_structure('X', pdb_file)				
	ref_struct = s[0]
	ref_chain = ref_struct[chain]
	ref_residues = []
	ref_res_nums = []
	for res in ref_chain:
		#print res
		ref_residues.append( res.resname )
		ref_res_nums.append(res.id[1])
	seq_array = convert_to_one_letter_code(ref_residues)
	seq = array_to_seq(seq_array)
	return seq, ref_res_nums

def convert_to_one_letter_code(residue_list):
	"""
	Takes a list of amino acids in letter-code and changes it to one-letter code.
	
	Args:
		residue_list: A list with an amino acid sequence represented in three-letter code
	
	Returns:
		one_letter_residue_list: A list whose elements are the amino acids in one-letter code
	
	"""
	one_letter_residue_list = []
	
	for res in residue_list:
		abbrev = resdict[res]
		one_letter_residue_list.append(abbrev)
	return one_letter_residue_list	
