import sys

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path):
    residues_positions = {}

    with open(pdb_file_path, 'r') as fd:
        file = fd.readlines()

        for line in file:
            if line.startswith('ATOM'):
                chain_id = line[21]
                residue_num = int(line[22:26])
                residue_name = line[17:20].strip()
                position_x = float(line[30:38].strip())
                position_y = float(line[38:46].strip())
                position_z = float(line[46:54].strip())

                residue_key = f"{residue_num}_{residue_name}"

                if chain_id not in residues_positions:
                    residues_positions[chain_id] = {}

                if residue_key not in residues_positions[chain_id]:
                    residues_positions[chain_id][residue_key] = []

                residues_positions[chain_id][residue_key].append((position_x, position_y, position_z))

    #calculating min distances for each pair of residues within each chain
    min_distances = {}

    for chain_id, residues in residues_positions.items():
        min_distances[chain_id] = {}

        residue_keys = list(residues.keys())
        for i in range(len(residue_keys) - 1):
            for j in range(i + 1, len(residue_keys)):
                min_distance = float('inf')
                #iterating over positions of residues i and j
                for position_i in residues[residue_keys[i]]:
                    for position_j in residues[residue_keys[j]]:
                        #calculating the Euclidean distance between position_i and position_j
                        distance = sum((a - b) ** 2 for a, b in zip(position_i, position_j))
                        min_distance = min(min_distance, distance)
                #storing the min distance for the pair residues
                min_distances[chain_id][(residue_keys[i], residue_keys[j])] = min_distance

    #calculating mean distance for each chain
    mean_distances = {}

    for chain_id, distances in min_distances.items():
        mean_distance = sum(distances.values()) / len(distances)
        mean_distances[chain_id] = mean_distance

    #the result is printed following the specified format, the mean distance for each chain with 4 decimal positions
    for chain_id, mean_distance in mean_distances.items():
        print(f"{chain_id}: {mean_distance:.4f}")

    return mean_distances

if __name__ == "__main__":
    import sys

    #checking if a PDB file path is provided as a command line argument
    if len(sys.argv) > 1:
        pdb_file_path = sys.argv[1]
    else:
        #reading pdb file from the standard input in case it is not provided as a command line argument
        pdb_file_path = None

    #calling the function and printing the mean distances
    calculate_pdb_chain_mean_minimum_distances(pdb_file_path)

#terminal implementation: python u235935_exercise_block1_part6.py /path/pdb/file