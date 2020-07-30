from Bio import SeqIO
from Bio import Seq
import csv

# Variables
num_chromosomes = 10
output_handle = "ChineseAmber_draft.pseudomol"
contig_file = "ChineseAmber.contigs.fasta"

# Initialize list of sequence objects based on number of chromosomes in genome
sequences = []
for x in range(1, num_chromosomes):
    sequences.append('')

# Extract data from sorted_best_hits (no threshold) file
with open('sorted_best_hit_coords_no_t.csv', newline='') as file:
    reader = csv.reader(file)
    coord_list = list(reader)
    coord_list.pop(0)  # Remove header information

for entry in coord_list:
    contig_name, starting_position, chromosome_num = entry

    # Get content of match and write it to the corresponding chromosome (HIGHLY INEFFICIENT TIME COMPLEXITY)
    for contig in SeqIO.parse(contig_file, "fasta"):
        if contig.id == contig_name:
            sequences[int(chromosome_num)] += contig.seq  # Append data to chromosome

# # Write to .fasta file
# with open("example.fasta", "w") as output_handle:
#     SeqIO.write(sequences, output_handle, "fasta")

print(sequences)
