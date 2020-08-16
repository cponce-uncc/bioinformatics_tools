"""
This program analyzes the quality of .fasta files and prints the results out. The analysis includes the following metrics:
    - Number of contigs
    - BP count
    - Largest contig (BP)
    - N50, N90, L50, and L90 statistics
    - Individual base counts and ratios

To use the program
Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.
"""
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import glob
import csv


# Function not created by me
def file_browser(ext=" "):
    return [f for f in glob.glob(f"*{ext}")]


def get_targeted_files():
    targeted_extensions = {}

    with open('targeted_extensions.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        for row in csv_reader:
            print(row, row[2])

            if row[2] == '1':  # Add file type to targeted dictionary if enabled in supp
                targeted_extensions.update({row[0]: row[1]})

    return targeted_extensions


def analyze_quality(file_path):
    contig_sizes = []
    num_contigs = 0  # Keeps track of number of contigs in file
    num_bp = 0  # bp = base pairs, keeps track of total base pairs
    largest_bp = 0  # Keeps track of contig with most base pairs
    a_count, t_count, c_count, g_count, n_count = 0, 0, 0, 0, 0

    for contig in SeqIO.parse(file_path, "fasta"):  # 3307 discovered records

        num_contigs += 1  # Iterate to account for a single continuous sequence in the assembly

        seq_length = len(contig.seq)  # Get length of each sequence
        contig_sizes.append(seq_length)
        num_bp += seq_length  # Add sequence length to total length

        # Updates largest_bp if necessary
        if seq_length > largest_bp:  # Check if contig length is greater than previous greatest length
            largest_bp = seq_length  # Reassigns largest_bp to updated value

        # Iterate count based on number of base-pairs found in each contig
        a_count += contig.seq.count('A')
        t_count += contig.seq.count('T')
        g_count += contig.seq.count('G')
        c_count += contig.seq.count('C')
        n_count += contig.seq.count('N')

    # Calculate  metrics
    contig_sizes.sort(reverse=True)
    contig_sum = sum(contig_sizes)
    contig_count = 0
    n50 = 0
    n90 = 0
    l50 = 0
    l90 = 0

    for contig in contig_sizes:
        contig_count += contig
        l50 += 1
        if contig_count >= contig_sum * .5:
            n50 = contig
            break

    contig_count = 0
    for contig in contig_sizes:
        contig_count += contig
        l90 += 1
        if contig_count >= contig_sum * .9:
            n90 = contig
            break

    return num_contigs, num_bp, largest_bp, n50, l50, n90, l90, a_count, t_count, g_count, c_count, n_count


def print_results(file_path):
    num_contigs, num_bp, largest_bp, n50, l50, n90, l90, a_count, t_count, \
        g_count, c_count, n_count = analyze_quality(file_path)

    print(f"Results for {file_path}... \n"
          f"Number of contigs: {num_contigs} \n"
          f"Bp count: {num_bp} ({round((num_bp / 732200000) * 100, 2)}% of reference) \n"
          f"Largest contig: {largest_bp} bp \n"
          f"N50: {n50} bp \n"
          f"N90: {n90} bp \n"
          f"L50: Contig no. {l50} (Located at {round((l50 / num_contigs) * 100, 2)}% of total contig count) \n"
          f"L90: Contig no. {l90} (Located at {round((l90 / num_contigs) * 100, 2)}% of total contig count)\n"
          f"'A' count: {a_count }\n"
          f"'T' count: {t_count }\n"
          f"'G' count: {g_count }\n"
          f"'C' count: {c_count }\n"
          f"'N' count: {(n_count)}\n"
          f"Total count: {a_count + t_count + g_count + c_count + n_count}\n")


def graph_contig_sizes(file_path, graph_type='bar'):
    contig_lengths = []
    contig_names = []
    line_name = file_path.split('.')[0]  # Get line name (text before period)

    for record in SeqIO.parse(file_path, "fasta"):  # Iterate through each contig in the file
        contig_lengths.append(len(record.seq))
        contig_names.append(len(record.id))

    if graph_type == 'histogram':  # Histogram plot
        plt.hist(contig_lengths, 25)
        plt.ylabel("Number of contigs in bin")  # Weird issue here, ask during call
        plt.xlabel("Size (Mbp)")
        plt.title(f'Contig length distribution across {file_path}')

    else:  # Bar graph default plot
        y_pos = np.arange(len(contig_lengths))
        plt.bar(y_pos, contig_lengths, align='center', alpha=0.5)
        plt.ylabel("Size (Mbp)")  # Weird issue here, ask during call
        plt.xlabel("Contig no.")
        plt.title(f'Contig lengths (Mbp) across {file_path}')

    plt.savefig(f'{line_name}_{graph_type}.png', dpi=500)  # Save figure based on line name and graph type
    plt.show()


if __name__ == "__main__":

    # Get targeted files through targeted_extensions.csv file
    targeted_files = {'.pseudomol': 'MUMmer'}

    # Analyze each file in the .csv file
    for extension in targeted_files.keys():
        files = file_browser(extension)

        for file in files:
            print_results(file)
