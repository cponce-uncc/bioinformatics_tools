"""
This program compares Mummer .coords files to the .fasta files they originated from. It will check for contig
representation in the .coords files, and find the best match (Highest query coverage) for each contig represented in the
.coords file. The results can be print out or be made into a CSV table. To use this program, use create_csv() or
print_results() and pass the filenames as arguments.

Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.
"""
from Bio import SeqIO
import csv


def get_assembly_contigs(file_path):
    contig_names = []  # Holds all contig names

    for record in SeqIO.parse(file_path, "fasta"):  # Iterate through each contig in the file
        contig_names.append(record.id)  # Add contig name to list

    return set(contig_names)  # Return set of contig names to remove redundant entries


def get_coords_info(file_path, threshold=0.0, return_coords=True):
    file = open(file_path, 'rt')
    line_no = 0
    all_query_hits = []
    contig_coords = {}
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

            if float(query_coverage) > threshold:
                coordinate_info = []

                if return_coords:
                    # Get S1data from .coords file
                    chromosome_coordinate_data = separated_contents[0].split(' ')

                    # Eliminate blank data in coordinate information
                    for entry in chromosome_coordinate_data:
                        if entry != '':
                            s1_data = entry  # S1 data, represents starting coordinate of match with respect to chr
                            break  # Stop search at first non-blank entry

                    # Get chromosome tag
                    chromosome_num = separated_contents[6] \
                        .split('\t')[0] \
                        .split('_')[1]

                    coordinate_info.append([s1_data, chromosome_num])

                    # print(contig_name, contig_coord)
                    contig_coords.update(
                        {contig_name: coordinate_info})  # Will throw error if there is headers are not line 0

                else:
                    if float(query_coverage) >= threshold:
                        all_query_hits.append([contig_name, query_coverage])

        line_no += 1  # Iterate line counter

    for entry in all_query_hits:
        if entry[0] in best_hits.keys():

            if best_hits.get(entry[0]) < entry[1]:  # Update dictionary if new entry is greater
                best_hits.update({entry[0]: entry[1]})

        else:  # Put entry into dictionary if contig has no entry yet
            best_hits.update({entry[0]: entry[1]})

    # Return data based on optional arguments passed to function
    if return_coords is True:
        return contig_coords

    else:
        return best_hits


def create_filtered_csv(new_file_name, filtered_results):
    contig_names = list(filtered_results.keys())
    coordinate_data = list(filtered_results.values())
    s1_data = []
    chromosome_nums = []

    for entry in coordinate_data:
        for contig_data in entry:
            s1_data.append(contig_data[0])
            chromosome_nums.append(contig_data[1])

    i = 0

    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Contig Name", "% Query Coverage"])

        for entry in filtered_results:
            writer.writerow([contig_names[i], s1_data[i], chromosome_nums[i]])
            i += 1

        file.close()


def create_sorted_csv(new_file_name, sorted_results):
    contig_names = []
    s1_data = []
    chromosome_nums = []

    for entry in sorted_results:
        contig_names.append(entry[0])
        s1_data.append(entry[1])
        chromosome_nums.append(entry[2])

    i = 0

    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Contig Name", "S1 Data", "Chromosome Tag"])

        for entry in filtered_results:
            writer.writerow([contig_names[i], s1_data[i], chromosome_nums[i]])
            i += 1

        file.close()


def create_quality_csv(assembly_file_path, coord_file_path, new_file_name):
    assembly_names = get_assembly_contigs(assembly_file_path)
    best_hits = get_coords_info(coord_file_path)

    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Contig Name", "Represented in .coords", "Best Query Hit"])

        for entry in assembly_names:
            if entry in best_hits.keys():
                writer.writerow([entry, "Yes", best_hits.get(entry)])
            else:
                writer.writerow([entry, "No", "N/A"])

        file.close()


def sort_best_hit_data(file_path):
    best_hit_data = []
    best_hit_csv = list(open(file_path, 'r', newline=''))
    best_hit_csv.pop(0)  # Remove headers
    for entry in best_hit_csv:
        hit_data = entry.split(',')
        hit_data[2] = hit_data[2].replace('\r\n', '')
        best_hit_data.append(hit_data)

    sorted_data = sorted(best_hit_data, key=lambda x: (float(x[2]), float(x[1])))
    return sorted_data


filtered_results = get_coords_info('ChineseAmber.coords', 0.0, return_coords=True)
create_filtered_csv('no_threshold_best_hits.csv', filtered_results)
sorted_results = sort_best_hit_data('no_threshold_best_hits.csv')
print(sorted_results)
create_sorted_csv('sorted_best_hit_coords_no_t.csv', sorted_results)
