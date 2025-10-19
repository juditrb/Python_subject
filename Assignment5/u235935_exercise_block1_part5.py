#Exercise 1
    #Exercise 1.1
def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield (identifier, sequence)
def get_proteins_ratio_by_residue_threshold(filename, residue, relative_threshold=0.03, absolute_threshold=10):
    def frequencies_calcul(protein, residue):
        absolute = protein.count(residue)
        relative = absolute / len(protein) if len(protein) > 0 else 0
        return absolute, relative

    total_proteins = 0
    protein_match = 0

    for identifier, sequence in FASTA_iterator(filename):
        total_proteins += 1
        if len(sequence) == 0:
            continue

        absolute, relative = frequencies_calcul(sequence, residue)
        if absolute >= absolute_threshold and relative >= relative_threshold:
            protein_match += 1

    return protein_match / total_proteins if total_proteins > 0 else 0

ratio = get_proteins_ratio_by_residue_threshold('example_fasta_file.fa', "A")
print("Ratio:", ratio)

    #Exercise 1.2
def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield (identifier, sequence)
def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
    
        with open(output_filename, "w") as out_file:
            
            for identifier, sequence in FASTA_iterator(filename):
                    
                    first_n_aa = sequence[:first_n] # The first/last residues of the protein sequence shown in the output depend on the value of the variable defined in the function
                    last_m_aa = sequence[-last_m:]
                    
                    aa_freq = {} #creating a dictionary to store amino acid frequencies
                    for aa in sequence:
                        if aa in aa_freq: # If the current amino acid is present in amino acid frequency dictionary, add one to the counter for that amino acid, if the amino acid is not in the dictionary, add it and set it to 1
                            aa_freq[aa] += 1
                        else: 
                            aa_freq[aa] = 1
                    aa_freq_string = ','.join([f"{aa}:{aa_freq[aa]}" for aa in aa_freq]) # For each amino acid in aa frequency dictionary, put together the amino acids with their frequency using a tabulator (the aa: frequency)
                
                    out_file.write(f"{identifier}\t{first_n_aa}\t{last_m_aa}\t{aa_freq_string}\n") # The output file will contain the protein_id, first_n, last_m and aa_frequency separated by a tabulator
                    
print_sequence_summary('example_fasta_file.fa', 'output_summary.txt', first_n=10, last_m=10)

#Exercise 2
def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence

def get_max_sequence_length_from_FASTA_file(filename):
    max_length = 0
    for _, sequence in FASTA_iterator(filename):
        sequence_length = len(sequence)
        max_length = max(max_length, sequence_length)
    
    print(max_length)

get_max_sequence_length_from_FASTA_file('example_fasta_file.fa')

#Exercise 3
def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
            
def get_min_sequence_length_from_FASTA_file(filename):
    min_length = float('inf') #the first sequence length has to be the new minimum
    for _, sequence in FASTA_iterator(filename): #for each sequence, it is calculated the length and the minimum length is updated
        sequence_length = len(sequence)
        min_length = min(min_length, sequence_length)
            
    print(min_length)

get_min_sequence_length_from_FASTA_file('example_fasta_file.fa')

#Exercise 4

def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
def get_longest_sequences_from_FASTA_file(filename):
    max_length = 0
    sequences_with_max_length = []

    for identifier, sequence in FASTA_iterator(filename):
        seq_length = len(sequence)
        #it is calculated the sequence length and then an if conditional is used to keep track of the maximum length along with the corresponding identifier and sequence
        if seq_length > max_length:
            max_length = seq_length
            sequences_with_max_length = [(identifier, sequence)]
        elif seq_length == max_length:
            sequences_with_max_length.append((identifier, sequence))
    #the sequences are sorted from longest to shortest by means of sorted function
    sorted_sequences = sorted(sequences_with_max_length, key=lambda x: x[0].lower())

    return sorted_sequences

result = get_longest_sequences_from_FASTA_file('example_fasta_file.fa')
print(result)

#Exercise 5

def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
def get_shortest_sequences_from_FASTA_file(filename):
    min_length = float('inf')
    sequences_with_min_length = []

    for identifier, sequence in FASTA_iterator(filename):
        seq_length = len(sequence)
        #it is calculated the sequence length and then an if conditional is used to keep track of the minimum length along with the corresponding identifier and sequence
        if seq_length == min_length:
            sequences_with_min_length.append((identifier, sequence))
        elif seq_length < min_length:
            min_length = seq_length
            sequences_with_min_length = [(identifier, sequence)]
    #the sequences are sorted from shortest to longest by means of sorted function
    return sorted(sequences_with_min_length, key=lambda x: x[0].lower())

result = get_shortest_sequences_from_FASTA_file('example_fasta_file.fa')
print(result)

#Exercise 6

#defining a function that iterates over a FASTA file, yield identifiers and sequences.
def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
            
def get_molecular_weights(filename):
    #initializing a molecular weight dictionary to store molecular weights for each protein
    molecular_weights = {}
    #dictionary containing molecular weights for each amino acids.
    amino_acid_weights = {
        "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
        "Q": 146.2, "E": 147.1, "G": 75.1, "H": 155.2, "I": 131.2,
        "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
        "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1
    }
    #iterate over identifiers and sequences from the FASTA_file.
    for identifier, sequence in FASTA_iterator(filename):
        #calculating the protein weight summing individual amino acid weights in the sequence
        protein_weight = round(sum(amino_acid_weights.get(aa, 0) for aa in sequence), 1)
        #it is stored the molecular weight in the dictionary
        molecular_weights[identifier] = protein_weight
        
    return molecular_weights

result = get_molecular_weights('example_fasta_file.fa')
print(result)


#Exercise 7

def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
            
def get_lowest_molecular_weights(filename):
    
    first_min_protein = None
    molecular_weights = get_molecular_weights(filename)
    min_molecular_weights = min(molecular_weights.values())
    
    for identifier, sequence in FASTA_iterator(filename):
        if molecular_weights[identifier] == min_molecular_weights:
            if first_min_protein is None:
                first_min_protein = (identifier, sequence)
    
    return first_min_protein

result = get_lowest_molecular_weights('example_fasta_file.fa')
print(result)

#Exercise 8


def FASTA_iterator(filename):
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield identifier, sequence
            
def get_mean_molecular_weight(filename):

    molecular_weights = get_molecular_weights(filename)
    mean_molecular_weights = round((sum(molecular_weights.values())/len(molecular_weights)), 1)
    
    for identifier, sequence in FASTA_iterator(filename):
    
        return mean_molecular_weights

result = get_mean_molecular_weight('example_fasta_file.fa')
print(result)