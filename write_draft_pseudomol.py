from Bio import SeqIO
from Bio.Seq import Seq
import csv

# Variables
num_chromosomes = 10
output_handle = "ChineseAmber_draft.pseudomol"
contig_file = "ChineseAmber.contigs.fasta"

# Initialize list of sequence objects based on number of chromosomes in genome
sequences = []
for x in range(0, num_chromosomes):
    sequences.append('')

# Extract data from sorted_best_hits (no threshold) file
with open('unsorted_test.csv', newline='') as file:
    reader = csv.reader(file)
    coord_data = list(reader)
    coord_data.pop(0)  # Remove header information

# Translate coords data to dict format
contig_data = {}
for contig in SeqIO.parse(contig_file, "fasta"):
    contig_data.update({contig.id: contig.seq})

x = 0  # Variable to track current contig

for entry in coord_data:

    contig_name, starting_position, chromosome_num = entry  # Unpack contig data
    contig_sequence = str(contig_data.get(contig_name))  # Get contig sequence data
    chromosome_sequence = str(sequences[(int(chromosome_num) - 1)])  # Obtain and
    sequences[int(chromosome_num) - 1] = chromosome_sequence + contig_sequence # Add new contig data to chromosome sequence
    x += 1  # Append tracker

    print(f"CONTIG NO. {x} on CHR. {chromosome_num}, POS. {starting_position} OUT OF {len(coord_data)} CONTIGS")

for entry in sequences:
    print("ENTRY NO. X LEN: ", len(entry))

with open("example.fasta", "w") as ofile:
    chromosome_num = 1
    for entry in sequences:

        ofile.write(f">chr {chromosome_num} len = {len(entry)}\n")  # Write sequence tag and length
        ofile.write(entry + '\n')  # Write sequence of chromosome with a newline
        chromosome_num += 1  # Iterate

    ofile.close()