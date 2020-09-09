"""
This program processes MUMmer .coords files. It will take the best hit (based on query coverage %) for each contig,
and then sort the data by chromosome starting position (S1 column), and output the results as a .csv containing
the contig ID, chromosome starting position, and chromosome number of the hit.

Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.

TASK LIST:
    - Filter hits outside of wrong chromosome out
    - Find average point in contig hits
"""
from Bio import SeqIO
import csv

line_file = "ChineseAmber.coords"  # Change to local file path for .coords file of the line
output_file = "ChineseAmber_best_qc_hits.csv"  # Change to desired name for output file (Best query coverage % hits)


def get_assembly_contigs(file_path):
    contig_names = []  # Holds all contig names

    for record in SeqIO.parse(file_path, "fasta"):  # Iterate through each contig in the file
        contig_names.append(record.id)  # Add contig name to list

    return set(contig_names)  # Return set of contig names to remove redundant entries


def get_coords_info(file_path, threshold=0.0):
    file = open(file_path, 'rt')
    line_no = 0
    all_query_hits = []
    best_hits = {}

    for line in file:  # Iterate through every line in the file

        if line_no > 2:  # Disregard header lines

            separated_contents = line.split('|')  # Split contents using the | character as the delimiter

            # Extract contig name for each line
            contig_name = separated_contents[6].split('\t')[1] \
                .replace('\n', ' ') \
                .strip()

            # Get query coverage for each hit
            query_coverage = separated_contents[5].split('   ')[2].strip()

            # Query length
            seq_lengths = separated_contents[4].split(' ')
            print(seq_lengths)

            x = 0
            query_length = 0
            for entry in seq_lengths:
                if entry.strip().isnumeric():
                    x += 1
                if (x == 2) and entry.strip().isnumeric():
                    query_length = int(entry.strip())

            print(query_length)

            # Calculate nucleotide coverage
            nucleotide_coverage = query_length * (float(query_coverage) / 100)

            # Get S1data from .coords file
            chromosome_coordinate_data = separated_contents[0].split(' ')

            # Eliminate blank data in coordinate information
            s1_data = ''  # Create variable for s1 data to prevent reference before assignment warning
            for entry in chromosome_coordinate_data:
                if entry != '':
                    s1_data = entry  # S1 data, represents starting coordinate of match with respect to chr
                    break  # Stop search at first non-blank entry

            # Get chromosome tag (RESOLVE ISSUE DURING CALL)
            if '_' in separated_contents[6]:
                chromosome_num = separated_contents[6] \
                    .split('\t')[0] \
                    .split('_')[1]

            else:  # Work here
                chromosome_num = separated_contents[6] \
                    .split(' ')[1] \
                    .split('\t')[0] \
                    .split('r')[1]

                if chromosome_num.isnumeric() and (int(chromosome_num) != 10):
                    chromosome_num = chromosome_num.replace('0', '')

                print(f'num {chromosome_num} end')

            # coordinate_info.append([s1_data, chromosome_num])

            if chromosome_num.isnumeric():
                all_query_hits.append([contig_name, query_coverage, s1_data, chromosome_num])

        line_no += 1  # Iterate line counter

    all_query_hits.sort()
    print(all_query_hits)

    for entry in all_query_hits:

        # Unpack values
        contig_name = entry[0]
        query_coverage = float(entry[1])
        s1_data = int(entry[2])
        chromosome_num = int(entry[3])

        if query_coverage > 0:  # Only writes contig to best_hits dictionary if threshold is reached

            if contig_name in best_hits.keys():  # Checks if contig is already represented in best hits dictionary

                if best_hits.get(contig_name)[0] < query_coverage:  # Update dict if entry is greater
                    best_hits.update({contig_name: [query_coverage, s1_data, chromosome_num]})

            else:  # Put entry into dictionary if contig has no entry yet
                best_hits.update({contig_name: [query_coverage, s1_data, chromosome_num]})

    return best_hits



def sort_best_hit_data(best_hit_data):
    best_hit_contigs = list(best_hit_data.keys())
    best_hit_positions = list(best_hit_data.values())
    all_data = []
    index = 0

    while index < len(best_hit_contigs):
        all_data.append([best_hit_contigs[index], best_hit_positions[index][0], best_hit_positions[index][1], best_hit_positions[index][2]])
        index += 1

    sorted_data = sorted(all_data, key=lambda x: (float(x[3]), float(x[2])))

    return sorted_data


def create_csv(new_file_name, contig_data):
    print('your data', contig_data)
    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)

        for entry in contig_data:
            print(entry)
            writer.writerow(entry)

        file.close()


filtered_results = get_coords_info(line_file)
sorted_coords = sort_best_hit_data(filtered_results)
create_csv(output_file, sorted_coords)
