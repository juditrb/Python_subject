
from sequence_dictionaries import protein_weights
from sequence_dictionaries import dna_weights, rna_weights
#print(protein_weights)
from sequence_dictionaries import protein_letters, dna_letters, rna_letters
from sequence_dictionaries import rna_table, dna_table
from sequence_dictionaries import rna_table_back
from sequence_dictionaries import rna_stop_codons, rna_start_codons, dna_stop_codons, dna_start_codons, dna_complement, rna_complement


class Sequence():
    
    alphabet = set() #creating an alphabet 
    
    def __init__(self, identifier, sequence, alphabet):
        self.__identifier = identifier
        self.__sequence = sequence
    
        if not set(self.__sequence).issubset(self.alphabet):
            invalid = set(self.__sequence) - self.alphabet
            raise ValueError(f"Impossible to create instance: {invalid} not possible")
        
    def get_identifier(self):
        return self.__identifier
    
    def get_sequence(self):
        return self.__sequence
    
    def get_mw(self): 
        mw = 0
        for aa in self.__sequence:
            mw += protein_weights[aa]
        return mw
        
        if isinstance(self, ProteinSequence):
            weigths_dict = protein_weights
        elif isinstance(self, DNASequence):
            weigths_dict = dna_weights
        elif isinstance(self, RNASequence):
            weigths_dict = rna_weights
            
        return sum(weigths_dict[monomer] for monomer in self.get_sequence())
    
    def has_subsequence(self, sequence_obj):
        return sequence_obj.get_sequence() in self.__sequence
    
    #defining behaviour for the classes defined
    
    def __len__(self):
        return len(self.get_sequence())
    
    def __eq__(self, other):
        return self.get_sequence() == other.get_sequence()
    
    def __ne__(self, other):
        return self.get_sequence() != other.get_sequence()
    
    def __add__ (self, other):
        if type(self) != type(other): #if the sequences are different type, error message is printed
            raise ValueError(f"The sequences must be the same type")
        
        idents = self.get_identifier() + "+" + other.get_identifier()
        seqs = self.get_sequence() + other.get_sequence()

        return self.__class__(idents, seqs, alphabet=self.alphabet)
    
    #Sequence[i], returning the sequence element at position i. Position 0 is the first one. 
    def __getitem__(self,res):
        return self.get_sequence()[res]
    
    #it is returned a boolean if the sequence has a specific item
    def __contains__(self, item):
        return item in self.get_sequence()
    
    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))
    
    #sorting sequences according to their molecular weight 
    def __lt__(self, other):
        return self.get_mw() < other.get_mw() 
    def __le__(self, other):
        return self.get_mw() <= other.get_mw() 
    def __gt__(self, other):
        return self.get_mw() > other.get_mw() 
    def __ge__(self, other):
        return self.get_mw() >= other.get_mw() 
    
    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))

class ProteinSequence(Sequence):
    alphabet = set(protein_letters)
    
class NucleotideSequence(Sequence):
    
    def translate(self):

        if isinstance(self, DNASequence):
            codon_dict = dna_table
            start_codon = dna_start_codons
            stop_codon = dna_stop_codons
        elif isinstance(self, RNASequence):
            codon_dict = rna_table
            start_codon = rna_start_codons
            stop_codon = rna_stop_codons
        else:
            raise ValueError('Invalid sequence type')
    
        sequence = self.get_sequence()
        aa_list = []

        has_start_index = [sequence.find(codon) for codon in start_codon \
            if codon in sequence]
        
        if has_start_index:
            start_i = min([sequence.find(codon) for codon in start_codon \
                if codon in sequence])
            for i in range(start_i, len(sequence), 3):
                codon = sequence[i: i+3]
                if codon in stop_codon or len(codon) != 3 :
                    break
                aa_list.append(codon_dict[codon])
            
            return ''.join(aa_list)

        else:
            raise ValueError('Object has no start codon')

#defining DNASequence class
class DNASequence(NucleotideSequence):
    alphabet = set(dna_letters)
    
    def transcribe(self):
        return ''.join([dna_complement[DNA] for DNA in self.get_sequence()])

#defining RNASequence class

class RNASequence(NucleotideSequence):
    alphabet = set(rna_letters)

    def reverse_transcribe(self):
        return ''.join([rna_complement[RNA] for RNA in self.get_sequence()])

#Example usage:

#Creating instances of sequences
protein_seq1 = ProteinSequence(identifier="P1", sequence="MVKLH", alphabet=set(protein_letters))
protein_seq2 = ProteinSequence(identifier="P2", sequence="MVKLH", alphabet=set(protein_letters))

protein_seq1 == protein_seq2
    #super(Nucleotide, self).__init__(ProteinSequence) #returning the parent class