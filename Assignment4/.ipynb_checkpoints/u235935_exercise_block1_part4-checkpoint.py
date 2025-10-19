
#Exercise1

def FASTA_iterator(fasta_filename):
    with open('example_fasta_file.fa', 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]

        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            sequence = protein.split("\n", 1)[1].replace("\n", "")
            yield (identifier, sequence)
        
for identifier, sequence in FASTA_iterator('example_fasta_file.fa'):
    print(f"({identifier}, {sequence})")

#Exercise2
    
def FASTA_iterator(fasta_filename):
    with open(fasta_filename, 'r') as fd:
        file = fd.read()
        proteins = file.split(">")[1:]
    
        for protein in proteins:
            identifier = protein.split("\n", 1)[0]
            yield (identifier)

def compare_fasta_file_identifiers(fasta_filenames_list):
    result = {
        "intersection": set(),
        "union": set(),
        "frequency": {},
        "specific": {}
    }
    all_identifiers = []
    file_identifiers = {}

    for filename in fasta_filenames_list:
        identifiers = set()

        for identifier in FASTA_iterator(filename):
            identifier_lower = identifier.lower()  #all identifiers are attached to the pre-initialized set, kept in lowercase in order to perform insensitive comparison between identifiers
            identifiers.add(identifier_lower)

            result["frequency"][identifier_lower] = result["frequency"].get(identifier_lower, 0) + 1 #it is updating the frequency dictionary, if the identifier is not present in the dictionary, it is returned 0 by default. Then it is summed 1 to the value obtained

        file_identifiers[filename] = identifiers
        all_identifiers.extend(identifiers) #it is adding all unique identifiers from the current file to the list of all identifiers

    #finding intersection and union of identifiers
    #the following line code is determining the set of identifiers that are common
    result["intersection"] = set.intersection(*map(set, file_identifiers.values())) #the map function is used to apply the set function to each element of file_identifiers.values()
    result["union"] = set(all_identifiers)

    #finding specific identifiers for each file
    for filename, identifiers in file_identifiers.items():
        specific_identifiers = identifiers.difference(result["intersection"])
        
        #the identifiers with frequency 1 are included in the "specific" dictionary
        specific_identifiers = {id for id in specific_identifiers if result["frequency"][id] == 1}
        
        result["specific"][filename] = specific_identifiers

    return result

fasta_files = ['uniprot_sprot_sample.fasta', 'uniprot_sprot_sample2.fasta', 'uniprot_sprot_sample3.fasta']
result_dict = compare_fasta_file_identifiers(fasta_files)

for key, value in result_dict.items(): #printing the contents of result_dict in a formatted way
    if key == "specific":
        print(f"{key}:")
        for filename, specific_identifiers in value.items(): #in 'specific' dictionary, the filenames are keys and specific identifiers are values
            print(f"  {filename}: {specific_identifiers}")
    else:
        print(f"{key}: {value}") #this part is executed for keys other than 'specific'