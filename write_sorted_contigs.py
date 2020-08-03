from Bio import SeqIO
import csv

# Rearrange contig records themselves instead of sorting into chromosomes

# Variables
num_chromosomes = 10
output_handle = "ChineseAmber_draft.pseudomol"
contig_file = "ChineseAmber.contigs.fasta"

# Initialize list of sequence objects based on number of chromosomes in genome
sequences = []
for x in range(0, num_chromosomes):
    sequences.append('')

# Extract data from sorted_best_hits (no threshold) file
with open('sorted_test.csv', newline='') as file:
    reader = csv.reader(file)
    coord_data = list(reader)
    coord_data.pop(0)  # Remove header information

# Translate coords data to dict format
contig_data = {}
for contig in SeqIO.parse(contig_file, "fasta"):
    contig_data.update({contig.id: contig.seq})

x = 1  # Counter variable for printing progress

with open("sorted_contigs.fasta", "w") as file:
    for entry in coord_data:
        contig_id, starting_position, chromosome_num = entry  # Unpack contig data
        contig_sequence = str(contig_data.get(contig_id))  # Get contig sequence data
        file.write(f">{contig_id} len={len(contig_sequence)} chr={chromosome_num}\n")  # Write sequence tag and length
        file.write(contig_sequence + '\n')

        print(f"CONTIG NO. {x} on CHR. {chromosome_num}, POS. {starting_position} OUT OF {len(coord_data)} CONTIGS")
        x += 1  # Append counter

    file.close()