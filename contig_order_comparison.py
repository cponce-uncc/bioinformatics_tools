"""
This program compares the order of contigs between two different files. To use the program, edit the variables that
hold the file paths and edit the arguments passed to the get_common_contigs() function on line 96.

Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.
"""

import csv
import difflib

# Files to compare
threshold_hits_file_one = "chunked_contigs_only.csv"
threshold_hits_file_two = "tx463_ChineseAmber_sorted.csv"
threshold_hits_file_three = 'Rio_ChineseAmber_sorted.csv'
threshold_hits_file_four = 'ChineseAmber_sorted_nucleotide_coverage.csv'
linkage_maps_file = 'ChineseAmber_Linkage_Map.csv'


def convert_to_csv(old_file_path, new_file_path):
    old_file = open(old_file_path, "r")

    with open(new_file_path, 'w', newline='') as file:
        writer = csv.writer(file)

        for line in old_file:
            old_content = line.split(' ')

            new_content = [old_content[1].replace('S', 'tig')
                               .split('_')[0], old_content[1].replace('S', 'tig') \
                               .split('_')[1], old_content[0]]

            writer.writerow(new_content)

    old_file.close()  # Close old .txt file
    file.close()  # Close new .csv file


def get_common_contigs(list_one, list_two):
    print('1', list_one)
    print('2', list_two)
    new_list_one = []
    new_list_two = []

    for entry in list_one:
        if entry in list_two:
            new_list_one.append(entry)

    for entry in list_two:
        if entry in list_one:
            new_list_two.append(entry)

    return new_list_one, new_list_two


def extract_linkage_contigs(linkage_maps_file):
    linkage_contigs = []
    linkage_map = list(open(linkage_maps_file, 'r'))
    # Extract all contigs represented in linkage map
    for line in linkage_map:
        line = line.replace("\n", '') \
                .split(",")[0]

        if line not in linkage_contigs:
            linkage_contigs.append(line)

    return linkage_contigs


def extract_coords_contigs(threshold_hits_file):
    threshold_hits = list(open(threshold_hits_file, 'r'))
    old_threshold_contigs = []
    # Extract all contigs represented in threshold .coords file
    for line in threshold_hits:
        line = line.replace("\n", '') \
            .split(",")[0] \

        if line not in old_threshold_contigs:  # Only add one of each contig to represented contig list
            old_threshold_contigs.append(line)

    old_threshold_contigs.pop(0)  # Remove header from data
    return old_threshold_contigs


def write_to_csv(list_one, list_two, new_file_name):
    position = 1
    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        for entry in list_one:
            writer.writerow([position, entry, list_two[position - 1]])
            position += 1


# Extract the data from each file and filter the data to only represent contigs common between both sources
linkage_contigs = extract_coords_contigs(threshold_hits_file_one)
threshold_contigs = extract_coords_contigs(linkage_maps_file)
new_threshold_contigs, linkage_contigs = get_common_contigs(threshold_contigs, linkage_contigs)

# Create SequenceMatcher object to compare two contig orders and print results
sequence_matcher = difflib.SequenceMatcher(None, linkage_contigs, new_threshold_contigs)
print(f"Similarity in order between updated sorted contigs\nand linkage map contigs: {sequence_matcher.ratio() * 100}%")

# Create a csv with the results of the order comparison
#write_to_csv(linkage_contigs, new_threshold_contigs, 'Rio_ChineseAmber_comparison.csv')

