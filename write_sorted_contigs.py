from Bio import SeqIO
import csv

# Rearrange contig records themselves instead of sorting into chromosomes

# Variables
num_chromosomes = 10
coord_file = "ChineseAmber.coords"  # Change to the local path of the .coords file of the line
contig_assembly_file = "ChineseAmber.contigs.fasta"  # Edit concatenated string if using assembler other than Canu
sorted_file = "ChineseAmber_updated_contigs.csv"  # Change to local path of output of contig_chunking.py
output_file = "ChineseAmber_sorted_contigs" + ".fasta"  # Change to desired name of output file

# Initialize list of sequence objects based on number of chromosomes in genome
sequences = []
for x in range(0, num_chromosomes):
    sequences.append('')

# Extract data from sorted_best_hits (no threshold) file
with open(sorted_file, newline='') as file:
    reader = csv.reader(file)
    coord_data = list(reader)

# Translate coords data to dict format
contig_data = {}
for contig in SeqIO.parse(contig_assembly_file, "fasta"):
    contig_data.update({contig.id: contig.seq})

x = 1  # Counter variable for printing progress

with open(output_file, "w") as file:
    for entry in coord_data:
        contig_id, starting_position, chromosome_num = entry  # Unpack contig data
        contig_sequence = str(contig_data.get(contig_id))  # Get contig sequence data
        file.write(f">{contig_id} len={len(contig_sequence)} chr={chromosome_num}\n")  # Write sequence tag and length
        file.write(contig_sequence + '\n')

        print(f"CONTIG NO. {x} on CHR. {chromosome_num}, POS. {starting_position} OUT OF {len(coord_data)} CONTIGS")
        x += 1  # Append counter

    file.close()
