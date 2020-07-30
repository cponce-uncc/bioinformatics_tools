from Bio import SeqIO
from Bio import Seq
import csv


# Variables
num_chromosomes = 10
output_handle = "ChineseAmber_draft.pseudomol"

# Initialize list of sequence objects based on number of chromosomes in genome
sequences = []
for x in range(1, num_chromosomes):
    sequences.append(Seq)

# Extract data from sorted_best_hits (no threshold) file
with open('sorted_best_hit_coords_no_t.csv', newline='') as file:
    reader = csv.reader(file)
    coord_list = list(reader)
    coord_list.pop(0)  # Remove header information

# Reformat list extracted from .csv into dictionary
coord_dict = {}

for entry in coord_list:
    key = entry[0]
    value = [entry[1], entry[2]]
    coord_dict.update({key: value})

# Import .contigs.fasta data

# # Write to .fasta file
# with open("example.fasta", "w") as output_handle:
#     SeqIO.write(sequences, output_handle, "fasta")

print(coord_dict)