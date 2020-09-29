"""
This program processes MUMmer .coords files. It will attempt to find clusters of hits for each contig to obtain
positional data (S1, E1). If no cluster is available, the program will find the hit with the best query coverage to
obtain the positional data. This is the second approach to accurately extracting positional data from the results of
pairwise alignment (.coords file).

Written by Cristian Ponce during Summer 2020 for the Cooper Lab (University of North Carolina, Charlotte)
at the North Carolina Research Campus.
"""
import csv
import statistics
import matplotlib.pyplot as plt

# EDIT THESE VARIABLE
line_file = "ChineseAmber.coords"  # Change to local file path of the .coords file of the line
best_hit_file = "ChineseAmber_best_qc_hits.csv"  # Change to local file path of coordinate_data_processing.py output
output_file = "ChineseAmber_updated_contigs.csv"  # Change to desired name of updated contigs file


def reject_outliers(data, percentage=10):
    data.sort()  # Sort dataset (ascending)

    # Determine how many entries to remove from dataset
    data_size = len(data)
    shave_amount_left = (data_size / 100) * percentage
    shave_amount_right = data_size - shave_amount_left

    # Remove disqualifying data
    new_data = []
    for entry in data:
        if ((data.index(entry) + 1) > shave_amount_left) and ((data.index(entry) + 1) < shave_amount_right):
            new_data.append(entry)

    return new_data


def eliminate_blank_data(coordinate_list, index=1):
    if index == 1:  # Extracts first entry in list of data
        for item in coordinate_list:
            if item != '':
                return item
    if index == 2:  # Extracts second entry in list of data
        counter = 0
        for item in coordinate_list:
            if item != '':
                counter += 1
                if counter == 2:
                    return item


def get_chunked_coords_info(file_path):
    file = open(file_path, 'rt')
    line_no = 0
    all_query_hits = []
    hits_by_contig = {}
    filtered_data = {}

    for line in file:  # Iterate through every line in the file

        if line_no > 2:  # Disregard header lines

            separated_contents = line.split('|')  # Split contents using the | character as the delimiter

            # Extract contig name for each line
            contig_name = separated_contents[6].split('\t')[1] \
                .replace('\n', ' ') \
                .strip()

            # Get S1data from .coords file
            s1_data = eliminate_blank_data(separated_contents[0].split(' '))

            e1_data = eliminate_blank_data(separated_contents[0].split(' '), 2)

            s2_data = eliminate_blank_data(separated_contents[1].split(' '))

            e2_data = eliminate_blank_data(separated_contents[1].split(' '), 2)

            query_coverage = separated_contents[5].split('   ')[2].strip()

            # Get chromosome tag
            if '_' in separated_contents[6]:
                chromosome_num = separated_contents[6] \
                    .split('\t')[0] \
                    .split('_')[1]

            else:
                chromosome_num = separated_contents[6] \
                    .split(' ')[1] \
                    .split('\t')[0] \
                    .split('r')[1]

                if chromosome_num.isnumeric() and (int(chromosome_num) != 10):
                    chromosome_num = chromosome_num.replace('0', '')

            if chromosome_num.isnumeric():
                all_query_hits.append([contig_name, s1_data, e1_data, s2_data, e2_data, chromosome_num, query_coverage])

        line_no += 1  # Iterate line counter

    for entry in all_query_hits:
        hits_by_contig.update({entry[0]: []})

    for entry in all_query_hits:
        current_data = hits_by_contig.get(str(entry[0]))
        current_data.append(entry[1:7])
        hits_by_contig.update({entry[0]: current_data})

    # Filter contigs by number of hits and chromosome concentration
    disqualified_contigs = []
    min_hits = 45  # Minimum number of hits contig must have to pass filter
    for contig in hits_by_contig.keys():
        all_contig_hits = hits_by_contig.get(contig)
        if len(all_contig_hits) < min_hits:
            disqualified_contigs.append(contig)

    # Create list of each contig before removal of disqualifying contigs
    all_contig_ids = hits_by_contig.keys()

    # Remove each contig that did not meet the specifications
    for contig in disqualified_contigs:
        hits_by_contig.pop(contig)

    # Iterate through each contig in the dictionary
    for contig in all_contig_ids:
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
            all_s1_data.append(int(entry[0]))
            all_s2_data.append(int(entry[2]))

        # Shows plots for a select contig before data filtering
        if contig == "tig00000009":
            plt.hist(all_s1_data)
            plt.title("Starting position distribution before filtering for tig00000009")
            plt.show()

        # Remove outlier from data
        filtered_s1_data = reject_outliers(all_s1_data)
        filtered_s2_data = reject_outliers(all_s2_data)

        # Get important information from filtered data
        first_s1 = int(filtered_s1_data[0])
        last_s2 = int(filtered_s2_data[-1])
        total_size = abs(last_s2 - first_s1)
        filtered_data.update({contig: [first_s1, last_s2, total_size, main_chromosome]})

    return filtered_data


# Combine chunked contigs and non-chunked contigs
def combine_best_hit_data(chunked_hit_dict, highest_query_file_path):
    combined_data_dict = {}
    # Open highest query coverage percentage file and write data to list
    query_coverage_file = open(highest_query_file_path, "r")
    query_data = []

    for line in query_coverage_file:
        content = line.replace("\n", "")\
                        .split(",")
        # print(content)
        query_data.append([content[0], content[2], content[3]])

    # Get each contig id
    contig_ids = []
    for entry in query_data:
        contig_ids.append(entry[0])

    # Convert chunked hit data to a list
    chunked_hit_list = []
    for contig in chunked_hit_dict.keys():
        contig_data = chunked_hit_dict.get(contig)
        chunked_hit_list.append([contig, contig_data[0], contig_data[3]])

    for entry in query_data:
        combined_data_dict.update({entry[0]: [entry[1], entry[2]]})

    # Overwrite contigs represented in chunked hit list with chunked data
    for entry in chunked_hit_list:
        combined_data_dict.update({entry[0]: [entry[1], entry[2]]})

    # Convert combined data back into list
    combined_data_list = []
    for contig in combined_data_dict.keys():
        contig_data = combined_data_dict.get(contig)
        combined_data_list.append([contig, contig_data[0], contig_data[1]])

    sorted_data = sorted(combined_data_list, key=lambda x: (float(x[2]), float(x[1])))

    return sorted_data


def create_csv(new_file_name, contig_data):
    # Converts data to list and then writes, if variable passed to function is a dict
    if isinstance(contig_data, dict):
        contig_data_list = []
        for contig in contig_data.keys():
            data = contig_data.get(contig)
            contig_data_list.append([contig, data[0], data[1], data[2], data[3]])

        with open(new_file_name, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Contig Name", "% Query Coverage"])

            for entry in contig_data_list:
                writer.writerow(entry)

            file.close()

    # Writes data normally if variable is a list
    else:
        with open(new_file_name, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Contig Name", "% Query Coverage"])

            for entry in contig_data:
                writer.writerow(entry)

            file.close()


chunked_results = get_chunked_coords_info(line_file)
print(chunked_results)
combined_data = combine_best_hit_data(chunked_results, best_hit_file)
print(combined_data)
create_csv(output_file, combined_data)
