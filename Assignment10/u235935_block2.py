import sys
import os
from sequence_dictionaries import protein_weights
from sequence_dictionaries import dna_weights, rna_weights
from sequence_dictionaries import protein_letters, dna_letters, rna_letters
from sequence_dictionaries import rna_table, dna_table
from sequence_dictionaries import rna_table_back
from sequence_dictionaries import rna_stop_codons, rna_start_codons, dna_stop_codons, dna_start_codons, dna_complement, rna_complement

    
class IncorrectSequenceLetter(Exception):
    def __init__(self, letter, class_name):
        self.letter = letter
        self.class_name = class_name
    
    def __str__(self):
        return "The sequence item %s is not found in the alphabet of class %s" % (self.letter, self.class_name)

class Sequence():
    
    alphabet = set() #creating an alphabet 
    
    def __init__(self, identifier, sequence, alphabet):
        self.__identifier = identifier
        self.__sequence = sequence

        if not all(letter in alphabet for letter in sequence):
            invalid_letter = next(letter for letter in sequence if letter not in alphabet)
            raise IncorrectSequenceLetter(invalid_letter, self.__class__.__name__)
                
        
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
        sys.stderr.write("Calculating sequence length")
        return len(self.get_sequence())
    
    def __eq__(self, other):
        sys.stderr.write("Comparing the length")
        return self.get_sequence() == other.get_sequence()
    
    def __ne__(self, other):
        sys.stderr.write("Comparing the instances of two sequences")
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
    
    
def FASTA_iterator(filename, class_name):
    
    class_map = {
        'ProteinSequence': ProteinSequence,
        'DNASequence': DNASequence,
        'RNASequence': RNASequence
    }
    
    if class_name not in class_map:
        raise ValueError(f"Invalid class_name: {class_name}. It must be one of {', '.join(class_map.keys())}")
    
    sequence_class = class_map[class_name]
    
    with open(filename, 'r') as fd:
        file = fd.read()
        chains = file.split(">")[1:]
    
        for chain in chains:
            try: 
                identifier = chain.split("\n", 1)[0]
                sequence = chain.split("\n", 1)[1].replace("\n", "")
            
                seq_instance = sequence_class(identifier, sequence, alphabet=sequence_class.alphabet)

                yield seq_instance
            
            except IncorrectSequenceLetter as e:
                sys.stderr.write("ValueError: {}\n".format(e))
                continue
                
def process_files(filepath, outputpath=None):
    output_sequences = []

    if not os.path.exists(filepath):
        print(f"Error: Input file '{filepath}' does not exist.")
        return
    
    #using the FASTA_iterator to read sequences from the file
    sequence_iterator = FASTA_iterator(filepath, 'DNASequence')

            #processing each sequence
    for sequence in sequence_iterator:
            #example of usage: printing the identifier and sequence
        identifier = sequence.get_identifier()
        length = sequence.get_sequence()
        molecular_weight = sequence.get_mw()

        print(f"{identifier}\t{length}\t{molecular_weight}")

            #adding the processed sequence to the list
        output_sequences.append(sequence)

    #optionally, write the processed sequences to an output file
    if outputpath:
        with open(outputpath, 'w') as output_file:
            for sequence in output_sequences:
                output_file.write(f"{sequence.get_identifier()} {sequence.get_sequence()} {sequence.get_mw()}\n")
    
def main():
    
    inputpath = None
    outputpath = None
    
    if len(sys.argv) == 1:
        inputpath = '.'
    elif len(sys.argv) == 2:
        inputpath = sys.argv[1]
    elif len(sys.argv) == 3:
        inputpath = sys.argv[1]
        outputpath = sys.argv[2]
    else:
        print("Usage: python script.py [inputpath] [outputpath]")
        sys.exit(1)
        
    if not os.path.exists(inputpath):
        print(f"Error: Input path '{inputpath}' does not exist.")
        sys.exit(1)
        
    process_files(inputpath, outputpath)

if __name__ == "__main__":
    main()
    
#terminal: python 
#Example usage:

#Creating instances of sequences
#protein_seq1 = ProteinSequence(identifier="P1", sequence="MVKLHB", alphabet=set(protein_letters))
#protein_seq2 = ProteinSequence(identifier="P2", sequence="MVKLH", alphabet=set(protein_letters))

#protein_seq1 == protein_seq2

#print(protein_seq1.get_sequence())
    #super(Nucleotide, self).__init__(ProteinSequence) #returning the parent class
    
# Example usage for Protein sequences
#protein_iterator = FASTA_iterator('example_fasta_file.fa', 'ProteinSequence')
#for protein_sequence in protein_iterator:
#    print(protein_sequence.get_identifier())
#    print(protein_sequence.get_sequence())
#    print(protein_sequence.get_mw())

# Example usage for DNA sequences
#dna_iterator = FASTA_iterator('example_fasta_file.fa', 'DNASequence')
#for dna_sequence in dna_iterator:
#    print(dna_sequence.get_identifier())
#    print(dna_sequence.get_sequence())
#    print(dna_sequence.transcribe())