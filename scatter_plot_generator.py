import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
import statistics

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


def plot_data(data_set_one, c='blue'):
    x = []
    y = []
    z = []

    for entry in contig_data:
        mbp_covered = int(entry[2]) / 1000000
        query_coverage = float(entry[1]) / 100
        y.append(query_coverage)
        x.append(mbp_covered)
        z.append(query_coverage * mbp_covered)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    csfont = {'fontname': 'Times New Roman'}

    ax.scatter(z, y, x, marker='o', c=c)
    plt.title("'Best hit' contig statistics plot", **csfont)
    plt.ylabel('Query coverage percentage (0-1)', **csfont)
    plt.xlabel('Covered nucleotides in Mbp', **csfont)
    ax.set_zlabel('Contig size in Mbp', **csfont)
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

plot_data(filtered_best_hit_data)
plot_data(filtered_linkage_data, c='red')

