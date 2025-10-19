import sys
import os
import os.path
import sequence_data
import urllib.request
import re
import argparse
import gzip
import random

class IncorrectSequenceLetter(ValueError):

    def __init__(self, letter, class_name):
        self.letter = letter
        self.class_name = class_name

    def __str__(self):
        return "The sequence item %s is not found in the alphabet of class %s" %( self.letter,
                                              self.class_name )



def FASTA_iterator(fasta_filename, sequence_class):

    fd = open(fasta_filename,"r")
    sequence = ""
    for line in fd:
        if line[0]==">":
            if len(sequence)>0:
                try:
                    yield sequence_class(identifier, sequence)
                except IncorrectSequenceLetter as e:
                    sys.stderr.write("%s\n" %e)

            identifier = line[1:].strip()
            sequence = ""
        else:
            sequence+=line.strip()
    fd.close()

    if len(sequence)>0:
        try:
            yield sequence_class(identifier, sequence)
        except IncorrectSequenceLetter as e:
            sys.stderr.write("%s\n" %e)


class Sequence(object):
    """
    Generic class to define a BioPolymer.
    A biopolymer is composed by an identifier and a sequence of monomers
    """

    alphabet = set()    # Class attribute to define the possibles monomers
    mw = {}         # Class attribute to define the monomer molecular weights

    def __init__(self, identifier, sequence):

        self.__identifier = identifier

        for letter in sequence:
            if letter not in self.alphabet:
                raise IncorrectSequenceLetter(letter, type(self).__name__)

        self.__sequence = sequence
        self.__mw = None

    def get_identifier(self):
        """
        Getter method to obtain the identifier of the sequence
        """
        return self.__identifier

    def get_sequence(self):
        """
        Getter method to get the sequence string
        """
        return self.__sequence

    def get_mw(self):
        """
        Calculate the molecular weight of the Sequence instance
        as the sum of the molecular weight of all the monomers
        Return a float number
        """
        if self.__mw is None:
            self.__mw = sum( self.mw[letter] for letter in self.__sequence )
        return self.__mw

    def has_subsequence(self, sequence_obj ):
        """
        Check if the sequence of sequence_obj is contained in 
        the Sequence object
        """
        return sequence_obj.get_sequence() in self.__sequence


    def __len__(self):
        return len(self.__sequence)


    def __eq__(self, other):
        return type(self) == type(other) and self.get_sequence() == other.get_sequence()


    def __add__(self, other):
        if self.__class__ == other.__class__:
            return self.__class__( identifier = "%s+%s" %(  self.get_identifier(),
                                    other.get_identifier() ),
                           sequence = self.get_sequence() + other.get_sequence() )
        else:
            raise TypeError("Not possible to concatenate different types of sequences")

    def __getitem__(self, key):
        return self.__sequence[key]

    def __lt__(self, other):
        return len(self)<len(other)
        
    def __hash__(self):
        return (self.get_identifier(), self.get_sequence()).__hash__()

    def __contains__(self, seq):
        return seq in self.get_sequence()


class ProteinSequence(Sequence):
    """
    Protein Sequence object. The monomers of ProteinSequence are
    the aminoacids
    """
    
    #Override all specific Class attributes for ProteinSequence
    alphabet = set(sequence_data.protein_letters)
    mw = sequence_data.protein_weights

    reverse_RNA = {'A': 'GCU', 'C': 'UGU', None: 'UAA', 'E': 'GAG', 'D': 'GAU', 'G': 'GGU', 'F': 'UUU', 'I': 'AUU', 'H': 'CAU', 'K': 'AAG', 'M': 'AUG', 'L': 'UUG', 'N': 'AAU', 'Q': 'CAG', 'P': 'CCU', 'S': 'UCU', 'R': 'CGU', 'T': 'ACU', 'W': 'UGG', 'V': 'GUU', 'Y': 'UAU'}
    reverse_DNA = {'A': 'GCT', 'C': 'TGT', None: 'TAA', 'E': 'GAG', 'D': 'GAT', 'G': 'GGT', 'F': 'TTT', 'I': 'ATT', 'H': 'CAT', 'K': 'AAG', 'M': 'ATG', 'L': 'TTG', 'N': 'AAT', 'Q': 'CAG', 'P': 'CCT', 'S': 'TCT', 'R': 'CGT', 'T': 'ACT', 'W': 'TGG', 'V': 'GTT', 'Y': 'TAT'}

    def reverse_translate_to_RNA(self):
        return RNASequence( identifier = self.get_identifier()+"_reverse_translated",
                    sequence = "".join( [ self.reverse_RNA[aa] for aa in self.get_sequence() ] ) )

    def reverse_translate_to_DNA(self):
        return DNASequence( identifier = self.get_identifier(),
                                    sequence = "".join( [ self.reverse_DNA[aa] for aa in self.get_sequence() ] ) )


    def get_uniprot_attribute(self, attribute_name ):
        """Retrives a specific attribute value from uniprot webpage"""

        if re.match("\w{6}",self.get_identifier()) is None:
            raise ValueError("%s does not correspond to a uniprot identifier" %self.get_identifier())
            
        url_fd = urllib.request.urlopen("http://www.uniprot.org/uniprot/%s.txt" %self.get_identifier())
    
        values = []

        for line in url_fd:
            line = line.decode("utf-8")
            if line.startswith("DR"):
                fields = line[2:].strip().split(";")
                if fields[0].lower() == attribute_name.lower():
                    values.append(fields[1].strip())
        return values


class NucleotideSequence(Sequence):
    """
    Nucleotide Sequence Object. The monomers of NucleotideSequence
    are not defined, as they can be ribonucleotides or deoxyribonucleotides
    """

    complement = {}
    start_codons = set()
    end_codons = set()

    translation_dict =  {}
    
    def translate(self):
        started = False
        translated_sequence = ""
        sequence = self.get_sequence()
        for i in range(0,len(sequence),3):
            codon = sequence[i:i+3]
            if started is False:
                if codon in self.start_codons:
                    started = True
                    translated_sequence = self.translation_dict[codon]
            elif codon in self.stop_codons:
                break
            else:
                translated_sequence+=self.translation_dict[codon]
        
        return ProteinSequence( identifier = self.get_identifier()+"_translated",
                    sequence = translated_sequence )


class DNASequence(NucleotideSequence):

    alphabet = set(sequence_data.dna_letters)
    mw = sequence_data.dna_weights

    complement = sequence_data.dna_complement
    stop_codons = set(sequence_data.dna_stop_codons)
    start_codons = set(sequence_data.dna_start_codons)

    translation_dict = sequence_data.dna_table
    
    def transcribe(self):
        return RNASequence( identifier = self.get_identifier()+"_transcribed", sequence = self.get_sequence().replace("T","U") )


class RNASequence(NucleotideSequence):
    alphabet = set(sequence_data.rna_letters)
    mw = sequence_data.rna_weights
    complement = sequence_data.rna_complement

    stop_codons = set(['UAA', 'UAG', 'UGA'])
    start_codons = set(['UUG', 'CUG', 'AUG'])

    translation_dict = sequence_data.rna_table

    def reverse_transcribe(self):
                return DNASequence( identifier = self.get_identifier()+"_transcribed", sequence = self.get_sequence().replace("U","T") )


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Given DNA FASTA file(s), calculate the sequence length and molecular weight of their corresponding translated proteins. Outputs the output sorted by sequence length")

    parser.add_argument(    '-i', '--input',
                dest= "infile",
                action= "store",
                default= "./",
                help="Input FASTA file or directory containing fasta file")

    parser.add_argument(    '-o', '--output',
                dest="outfile",
                action="store",
                default= None,
                help="Ouput file. If not defined, it prints output to standard output.")

    parser.add_argument(    '-v', '--verbose',
                dest="verbose",
                action="store_true",
                default= False,
                help= "Print progression log to standard error")

    parser.add_argument(    '-p', '--pattern',
                dest="pattern",
                action="store",
                default=None,
                help="Regular expression pattern to search in the translated sequences")

    parser.add_argument(    '-r', '--random',
                            dest="random_output_num",
                action="store",
                default= None,
                type= int,
                help="Random how many sequence to print")

    options = parser.parse_args()


    # CAPTURING THE INPUT FILE(s)
    input_path = options.infile

    if os.path.isdir(input_path):
        list_files = [ os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith(".fa") or f.endswith(".fasta") ]
    elif os.path.isfile(input_path):
        list_files = [ input_path ]
    else:
        list_files = []

    if options.verbose:
        sys.stderr.write("%d FASTA files found.\n" %len(list_files))


    # READING THE SEQUENCE FILES AND STORE THEM IN A LIST
    output_list = []

    if options.pattern is not None:
        pattern_re = re.compile(options.pattern)


    for filename in list_files:
        for dna in FASTA_iterator(filename, DNASequence):
            prot = dna.translate()
            if options.pattern is None:
                output_list.append( ( prot.get_identifier(), len(prot), prot.get_mw()) )
            elif pattern_re.search(prot.get_sequence()):        
                output_list.append( ( prot.get_identifier(), len(prot), prot.get_mw()) )
            
        if options.verbose:
            sys.stderr.write("%s finished.\n" %filename)

    if options.verbose:
        sys.stderr.write("%s sequences found.\n" %len(output_list))

    if options.outfile is None:
        out_fd = sys.stdout
    else:
        if options.outfile.endswith(".gz"):
            out_fd = gzip.open(options.outfile, "wt")
        else:
            out_fd = open(options.outfile, "w")

    # WRITE OUTPUT
    if options.random_output_num is not None:
        if options.random_output_num > len(output_list):
            sys.stderr.write("Number of required sequences to be printed is greater than the number of available sequences\n")
            sys.exit(1)

        output_to_print = random.sample(output_list, options.random_output_num)
        
    else:
        output_to_print = output_list

    # SORT THE SEQUENCES BY THE LENGTH
    if options.verbose:
        sys.stderr.write("Sorting the sequences...\n")
    output_to_print.sort(key=lambda x: x[1], reverse=True)
    if options.verbose:
        sys.stderr.write("Sort process finished.\n")

    for identifier, length, mw in output_to_print:
        out_fd.write("%s\t%s\t%s\n" %(  identifier,
                        length,
                        mw))

    out_fd.close()

    if options.verbose:
        sys.stderr.write("Program finished correctly.")


    
