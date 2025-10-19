from data import molecular_weights

def FASTA_iterator(fasta_filename):
    """
    Generator function to read multiline fasta files

    Yield a tuple (id, sequence) for each sequence in the FASTA file
    """

    fd = open(fasta_filename)

    sequence = ""

    for line in fd:

        if line[0] == ">":

            if len(sequence)>0:
                yield((identifier, sequence))

            sequence = ""
            identifier = line[1:].strip()

        else:
            sequence += line.strip()

    if len(sequence)>0:
        yield((identifier, sequence))

    fd.close()


def get_max_sequence_length_from_FASTA_file(fasta_filename):
    """
    Return the maximum  length of a protein in a protein FASTA File """

    return max( map(lambda x: len(x[1]), FASTA_iterator(fasta_filename)) )


def get_min_sequence_length_from_FASTA_file(fasta_filename):

    return min(map(lambda x: len(x[1]), FASTA_iterator(fasta_filename)))


def get_longest_sequences_from_FASTA_file(fasta_filename):
    """

    """

    current_max_length = 0
    current_max_seqs = []

    for seq_id, sequence in FASTA_iterator(fasta_filename):
        if len(sequence) > current_max_length:
            current_max_seqs = [(seq_id, sequence)]
            current_max_length = len(sequence)
        elif len(sequence) == current_max_length:
            current_max_seqs.append((seq_id, sequence))

    return sorted(current_max_seqs, key=lambda x: x[0].upper())


def get_shortest_sequences_from_FASTA_file(fasta_filename):
    """

    """

    current_min_length = 1000000000
    current_min_seqs = []

    for seq_id, sequence in FASTA_iterator(fasta_filename):
        if len(sequence) < current_min_length:
            current_min_seqs = [(seq_id, sequence)]
            current_min_length = len(sequence)
        elif len(sequence) == current_min_length:
            current_min_seqs.append((seq_id, sequence))

    return sorted(current_min_seqs, key=lambda x: x[0].upper())


def get_molecular_weight(sequence):
    """
    """

    return sum( [molecular_weights.aminoacid_mw.setdefault(aa,0) for aa in sequence] )


    aminoacids =  molecular_weights.aminoacid_mw.keys()

    for aa in aminoacids:
        mw += molecular_weights.aminoacid_mw[aa] * sequence.count(aa)

    return mw

def get_molecular_weights(fasta_filename):
    """
    """

    proteins_mw = {}

    for seq_id, sequence in FASTA_iterator(fasta_filename):
        proteins_mw[seq_id] = get_molecular_weight(sequence)

    return proteins_mw



def get_sequence_with_min_molecular_weight( fasta_filename ):

    sorted_mws = sorted(FASTA_iterator(fasta_filename), key=lambda x: get_molecular_weight(x[1]))
    return sorted_mws[0]



def get_mean_molecular_weight( fasta_filename ):
    """
    Returns an integer corresponding to the mean of the molecular weights of all the proteins
    """
    mws = get_molecular_weights(fasta_filename)

    return sum(mws.values()) / len(mws)
