import csv
from Bio import SeqIO
import matplotlib.pyplot as plt

contig_file = "ChineseAmber.contigs.fasta"
best_hit_file = "ChineseAmber_sorted.csv"
linkage_map_file = "linkage_map_s1_data.csv"


def extract_data(file_name):
    # Extract data from sorted_best_hits (no threshold) file
    with open(file_name, newline='') as file:
        reader = csv.reader(file)
        extracted_data = list(reader)
        extracted_data.pop(0)  # Remove header information
        return extracted_data


def plot_data(data_set_one, c='blue', line_name='insert line name here'):
    x = []
    y = []
    z = []

    for entry in data_set_one:
        contig_size = int(entry[2]) / 1000000
        query_coverage = float(entry[1]) / 100
        x.append(contig_size)
        y.append(query_coverage)
        z.append(query_coverage * contig_size)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    csfont = {'fontname': 'Times New Roman'}

    ax.scatter(x, y, z, marker='o', c=c)
    plt.title(f"Size, query coverage, and nucleotide coverage in {line_name} contigs", **csfont)
    plt.xlabel('Contig size (Mbp)', **csfont)
    plt.ylabel('Query coverage % (0-1)', **csfont)
    ax.set_zlabel('Nucleotides covered (Mbp)', **csfont)
    plt.show()


best_hit_data = extract_data(best_hit_file)
linkage_data = extract_data(linkage_map_file)

# Get linkage contigs
linkage_contigs = []
for entry in linkage_data:
    linkage_contigs.append(entry[0])

# Translate coords data to dict format
contig_data = {}
for contig in SeqIO.parse(contig_file, "fasta"):
    contig_data.update({contig.id: len(contig)})

filtered_best_hit_data = []
for entry in best_hit_data:
    if entry[0] not in linkage_contigs:
        filtered_best_hit_data.append([entry[0], entry[1], contig_data.get(entry[0])])

filtered_linkage_data = []
for entry in linkage_data:
    filtered_linkage_data.append([entry[0], entry[1], contig_data.get(entry[0])])

plot_data(filtered_best_hit_data, line_name='C. Amber')
plot_data(filtered_linkage_data, c='red', line_name='C. Amber')
