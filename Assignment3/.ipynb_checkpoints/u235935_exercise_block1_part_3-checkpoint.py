# Create a function that, given a multi-line protein FASTA file (fasta_filename) and a
# “sub-sequences” file (subsequences_filename) (one sequence in each line),
# calculates the proportion of proteins in the FASTA file containing at least N-times
#(number_of_repetitions) each of the sub-sequences (exactly equal). Save it in an
#output file with the specified format, ordered by the proportion value
#(descending order)

def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename, output_filename, number_of_repetitions=3):

    with open('sequence_fragments.txt', 'r') as file:
        subseqdic = {}
        for line in file:
            subseqdic[line.strip()] = 0

    with open('example_fasta_file.fa', 'r') as fd:
        fastaf = fd.read().split('>')[1:]
        proteins = ["".join(protein.split("\n")[1:]) for protein in fastaf]
    
    for protein in proteins:
        for subseq in subseqdic.keys():
            if protein.count(subseq) >= 3:
                subseqdic[subseq] += 1
    

    total_proteins = len(proteins)

    sorted_subseqdic = dict(sorted(subseqdic.items(), key=lambda x:x[1], reverse=True))

    with open('output.txt', 'w') as out_file:
        out_file.write(f"#Number of proteins: {total_proteins:>{43}}\n")
        out_file.write(f"#Number of subsequences: {len(subseqdic.keys()):>{39}}\n")
        out_file.write("#Subsequence proportions:\n")
        for subseq, count in sorted_subseqdic.items():
            proportion = round(count / total_proteins, 4)
            out_file.write(f"{subseq} {count:>{20}} {proportion:>{40}}\n")

calculate_aminoacid_frequencies('example_fasta_file.fa', 'sequence_fragments.txt', 'output.txt', 3)