
#1. Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
#float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
#or equal than a given threshold provided as an argument named “relative_threshold” and
#having an absolute frequency of the same residue higher or equal than a given threshold
#provided as an argument named “absolute_threshold” for a given residue. The function
#should be named as follows, with the same arguments definition:

def get_proteins_ratio_by_residue_threshold(filename, residue, relative_threshold=0.03, absolute_threshold=10):
    def frequencies_calcul(protein, residue): # 'Frequencies_calcul' function is defined to calculate the relative/absolute frequency depending on protein and residue variables
        absolute = protein.count(residue) # The absolute frequency is calculated through counting the quantity of a specific residue in a protein 
        relative = absolute / len(protein) if len(protein) > 0 else 0 # The relative frequency is calculated by the division of the absolute frequency and the protein length
        return absolute, relative
    
    with open(filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:] # This command is spliting the entire file using ">" as a delimiter, to diferenciate one protein to another (next to the ">" symbol the protein identifier and the protein sequence are found)
            #each "proteins" element is a protein sequence entry
        proteins = [protein.split("\n", 1)[1].replace("\n", "") for protein in proteins if protein] # Each protein entry is split into two parts: the header line and the sequence, [1] is selecting the second element of the split
        #The "replace" part concatenates the sequence into a continuous string and the "for" loop is iterating over each protein, only if it is not empty
        #print(count_proteins)
        #print(proteins)
        
        total_proteins = len(proteins)
        if total_proteins == 0:
            return 0
        
        protein_match = 0 #Initialization of protein_match variable
        
        for protein in proteins: # In this loop it is calculated the absolute and relative frequencies and they are compared to the absolute/relative threshold previously defined
            absolute, relative = frequencies_calcul(protein, residue)
            if absolute >= absolute_threshold and relative >= relative_threshold: #If both premises are met, the "protein match" variable increments the counter plus 1
                protein_match += 1
                
        return protein_match / total_proteins if total_proteins > 0 else 0
        
ratio = get_proteins_ratio_by_residue_threshold('example_fasta_file.fa', "M") #It is returned the ratio value depending on the residue we are interested in and the fa/fasta file used
print("Ratio:", ratio)

#2. Given a protein FASTA file (filename), save on a output file named output_filename the
#protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency
#in the protein of all the aminoacids found in the protein (the aminoacids that do not appear
#in the protein should not be shown). The fields must be separated by a tabulator, and one
#protein by line.
def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
    
    with open(filename, 'r') as fd:
        proteins = fd.read().split(">")[1:]
        #print(proteins)
    
        with open(output_filename, "w") as out_file:
            for protein in proteins:
                if protein:
                    lines = protein.split("\n") # The 'protein' string is including the content of a protein entry including the identifier and the sequence lines
                    protein_id = lines[0]  # The protein identifier is in the first line 'lines[0]'
                    #print(protein_id)
                    sequence = ''.join(lines[1:]) # The sequence variable stores the protein sequences without the identifier, it goes from lines[1] to the end of the corresponding protein sequence
                    #print(sequence)
                    
                    first_n_aa = sequence[:first_n] # The first/last residues of the protein sequence shown in the output depend on the value of the variable defined in the function
                    last_m_aa = sequence[-last_m:]
                    
                    aa_freq = {} #creating a dictionary to store amino acid frequencies
                    for aa in sequence:
                        if aa in aa_freq: # If the current amino acid is present in amino acid frequency dictionary, add one to the counter for that amino acid, if the amino acid is not in the dictionary, add it and set it to 1
                            aa_freq[aa] += 1
                        else: 
                            aa_freq[aa] = 1
                    aa_freq_string = ','.join([f"{aa}:{aa_freq[aa]}" for aa in aa_freq]) # For each amino acid in aa frequency dictionary, put together the amino acids with their frequency using a tabulator (the aa: frequency)
                
                    out_file.write(f"{protein_id}\t{first_n_aa}\t{last_m_aa}\t{aa_freq_string}\n") # The output file will contain the protein_id, first_n, last_m and aa_frequency separated by a tabulator
                    
print_sequence_summary('example_fasta_file.fa', 'output_summary.txt', first_n=10, last_m=10)