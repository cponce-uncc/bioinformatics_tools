"""
This program processes MUMmer .coords files. It will take the best hit (based on query coverage %) for each contig,
and then sort the data by chromosome starting position (S1 column), and output the results as a .csv containing
the contig ID, chromosome starting position, and chromosome number of the hit.

Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.
"""
from Bio import SeqIO
import csv
import statistics
import numpy as np


def reject_outliers(data, m=2):
    
    return data[abs(data - np.mean(data)) < m * np.std(data)]


def get_assembly_contigs(file_path):
    contig_names = []  # Holds all contig names

    for record in SeqIO.parse(file_path, "fasta"):  # Iterate through each contig in the file
        contig_names.append(record.id)  # Add contig name to list

    return set(contig_names)  # Return set of contig names to remove redundant entries


def eliminate_blank_data(coordinate_list, index=1):
    if index == 1:
        for item in coordinate_list:
            if item != '':
                return item
    if index == 2:
        counter = 0
        for item in coordinate_list:
            if item != '':
                counter += 1
                if counter == 2:
                    return item


def get_coords_info(file_path, threshold=0.0):
    file = open(file_path, 'rt')
    line_no = 0
    all_query_hits = []
    hits_by_contig = {}

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
            #print(seq_lengths)

            x = 0
            query_length = 0
            for entry in seq_lengths:
                if entry.strip().isnumeric():
                    x += 1
                if (x == 2) and entry.strip().isnumeric():
                    query_length = int(entry.strip())

            # print(query_length)

            # Calculate nucleotide coverage
            nucleotide_coverage = query_length * (float(query_coverage) / 100)

            # Get S1data from .coords file
            s1_data = eliminate_blank_data(separated_contents[0].split(' '))

            e1_data = eliminate_blank_data(separated_contents[0].split(' '), 2)

            s2_data = eliminate_blank_data(separated_contents[1].split(' '))

            e2_data = eliminate_blank_data(separated_contents[1].split(' '), 2)

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

                #print(f'num {chromosome_num} end')

            # coordinate_info.append([s1_data, chromosome_num])

            if chromosome_num.isnumeric():
                all_query_hits.append([contig_name, s1_data, e1_data, s2_data, e2_data, chromosome_num])

        line_no += 1  # Iterate line counter

    for entry in all_query_hits:

        hits_by_contig.update({entry[0]: []})

    for entry in all_query_hits:
            print('test', entry[1:6])
            current_data = hits_by_contig.get(str(entry[0]))
            current_data.append(entry[1:6])
            hits_by_contig.update({entry[0]: current_data})

    # Filter contigs by number of hits and chromosome concentration
    disqualified_contigs = []
    min_hits = 45  # Minimum number of hits contig must have to pass filter
    for contig in hits_by_contig.keys():
        all_contig_hits = hits_by_contig.get(contig)
        if len(all_contig_hits) < min_hits:
            disqualified_contigs.append(contig)

    # Remove each contig that did not meet the specifications
    for contig in disqualified_contigs:
        hits_by_contig.pop(contig)

    # Iterate through each contig in the dictionary
    for contig in hits_by_contig.keys():
        all_contig_hits = hits_by_contig.get(contig)

        # Count occurrence of hits across all chromosomes
        chromosome_nums = []

        for entry in all_contig_hits:
            chromosome_nums.append(int(entry[4]))

        # Find main chromosome
        main_chromosome = statistics.mode(chromosome_nums)

        # Eliminate any hit not on the main contig
        updated_contigs = []
        for entry in all_contig_hits:
            if int(entry[4]) == main_chromosome:
                updated_contigs.append(entry)

        # Eliminate any hit that has outlier with respect to the chromosome position data
        all_s1_data = []
        all_s2_data = []

        for entry in updated_contigs:
            print("EBTW", entry)
            all_s1_data.append(entry[0])
            all_s2_data.append(entry[2])

        print(all_s1_data, all_s2_data)

        # Convert to numpy array
        all_s1_data_array = np.asarray(all_s1_data)
        all_s2_data_array = np.asarray(all_s2_data)

        # Remove outlier from data
        reject_outliers(all_s1_data_array, 2)
        reject_outliers(all_s2_data_array, 2)

        print(all_s1_data_array)


        # Sort contigs by s1 data
        sorted_contigs = (sorted(updated_contigs, key=lambda x: (float(x[0]))))

        # Update dictionary data
        hits_by_contig.update({contig: sorted_contigs})


    print('final:', hits_by_contig)
    return list(hits_by_contig)


def update_best_hits(best_hits, linkage_file_path):
    # Import linkage map file
    linkage_coordinate_info = []
    with open(linkage_file_path, 'r') as file:
        for line in file:
            linkage_coordinate_info.append(line)

    # Override best hits file
    for entry in linkage_coordinate_info:
        contig_name, query_coverage, s1_data, chromosome_num = entry.split(',')
        best_hits.update({contig_name: [float(query_coverage), int(s1_data), int(chromosome_num)]})

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
        writer.writerow(["Contig Name", "% Query Coverage"])

        for entry in contig_data:
            print('entry: ', entry)
            writer.writerow(entry)

        file.close()


filtered_results = get_coords_info('ChineseAmber.coords')
# create_csv('second_approach_contigs.csv', filtered_results)

