import sys
import math

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path=None):
  """
  """

  if pdb_file_path is None:
    input_fd = sys.stdin
  else:
    input_fd = open(pdb_file_path)

  pdb_data = {}
  result = {}

  for line in input_fd:
    if line.startswith("ATOM"):
      chain = line[21]
      residue = line[22:26]
      coordinates = (float(line[28:38]), float(line[38:46]), float(line[46:54]))
      pdb_data.setdefault(chain, {}).setdefault(residue, []).append(coordinates)

  for chain, residues_dict in pdb_data.items():
    distances =[]
    residues = list(residues_dict.keys())
    for residue1_index in range(len(residues)):
      for residue2_index in range(residue1_index):
        residue_pair_distances = []
        for atom1 in residues_dict[residues[residue1_index]]:
          for atom2 in residues_dict[residues[residue2_index]]:
            x1, y1, z1 = atom1
            x2, y2, z2 = atom2
            distance = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
            residue_pair_distances.append(distance)
        
        distances.append(min(residue_pair_distances))
        

    #print(chain, len(residues), len(distances))

    if len(distances)>0:
      result[chain] = sum(distances)/len(distances)

  return result


if __name__ == "__main__":

  if (len(sys.argv)>1):
    input_file = sys.argv[1]
  else:
    input_file = None

  result = calculate_pdb_chain_mean_minimum_distances(input_file)

  for chain in result:
      print("%s: %.4f" %(chain, result[chain]))


