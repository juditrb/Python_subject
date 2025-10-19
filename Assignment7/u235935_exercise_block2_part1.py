#Exercise 1

class Protein(): #defining a new class object
    def __init__(self, identifier, sequence): #to create a new instance of a class it is used __init__
        self.identifier = identifier
        self.sequence = sequence
    def get_identifier(self): #representing identifier instance of the protein class by means of self. It is also represented the sequence, molecular weight, the subsequence and the length
        return self.identifier
    def get_sequence(self):
        return self.sequence
    def get_mw(self):
        mw = 0
        amino_acid_weights = {
        "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
        "Q": 146.2, "E": 147.1, "G": 75.1, "H": 155.2, "I": 131.2,
        "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
        "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1
        }
        for aa in self.sequence:
            mw += amino_acid_weights[aa]
            return mw
    def has_subsequence(self, Protein):
        return Protein.get_sequence() in self.sequence
    def get_length(self):
        return len(self.sequence)
    
a = Protein('Protein 1', 'ADEFWGHICPFVMWN')

print(a.get_identifier())
print(a.get_sequence())
print(a.get_mw())
print(a.has_subsequence(Protein('b', 'WGHI')))
print(a.get_length())

#Exercise 2

def FASTA_iterator(filename):
    
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")

            yield Protein(identifier, sequence) #it is using the yield keyword, which indicates that it's inside a generator function

#for each protein of the fasta_file, it is printed the identifier, the sequence, the molecular weight, if the subsequence defined is inside this protein and the length of the protein
            #the instances of class protein are printed for each protein of the fasta_file
for protein_obj in FASTA_iterator('example_fasta_file.fa'):
    print(f'{protein_obj.get_identifier()},{protein_obj.get_sequence()},{protein_obj.get_mw()},{protein_obj.has_subsequence(Protein("Protein 2", "WGHI"))},{protein_obj.get_length()}')