import molecular_weights
import sys

def FASTA_iterator(fasta_filename):

	fd = open(fasta_filename,"r")
	sequence = ""
	for line in fd:
		if line[0]==">":
			if len(sequence)>0:
				try:
					yield Protein(identifier, sequence)
				except:
					pass
			identifier = line[1:].strip()
			sequence = ""
		else:
			sequence+=line.strip()
	fd.close()

	if len(sequence)>0:
		try:
			yield Protein(identifier, sequence)
		except:
			pass


class Protein(object):
	"""
	Class to define a Protein
	"""

	def __init__(self, identifier, sequence):
		"""
		Initializes the protein instance with a given identifier and sequence
		"""
		self.identifier = identifier
		self.sequence = sequence

	def get_identifier(self):
		"""
		Getter method to obtain the identifier of the protein
		"""
		return self.identifier

	def get_sequence(self):
		"""
		Getter method to get the sequence string
		"""
		return self.sequence

	def get_mw(self):
		"""
		Calculate the molecular weight of the protein instance
		as the sum of the molecular weight of all the monomers
		Return a float number
		"""
		return sum( molecular_weights.aminoacid_mw.setdefault(aa, 0) for aa in self.sequence )

	def has_subsequence(self, subsequence):
		"""
		Check if the sequence of sequence_obj is contained in
		the Sequence object
		"""
		return subsequence.get_sequence() in self.sequence

	def get_length(self):
		"""
		Return and integer corresponding to the length of the sequence.
		"""
		return len(self.sequence)

if __name__ == "__main__":

   for protein in FASTA_iterator(sys.argv[1]):
       print(protein.get_length(protein))
       print("%s (%d): %.4f" %(protein.get_identifier(), protein.get_length(), protein.get_mw()))


