import csv

linkage_maps_file = "ChineseAmber_Linkage_Map.csv"
coords_file = "ChineseAmber.coords"
import numpy as np
import matplotlib.pyplot as plt


def coord_to_dict(coords_file):
    file = open(coords_file, 'rt')
    line_no = 0
    coords_dict = {}

    for line in file:  # Iterate through every line in the file

        if line_no > 2:  # Disregard header lines

            separated_contents = line.split('|')  # Split contents using the | character as the delimiter

            # Extract contig name for each line
            contig_name = separated_contents[6].split('\t')[1] \
                .replace('\n', ' ') \
                .strip()

            # Get query coverage for each hit
            query_coverage = separated_contents[5].split('   ')[2].strip()

            # Get S1 data
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

            # Update if contig is new or if the new contig is of a higher coverage
            if (contig_name not in coords_dict.keys()) or (float(query_coverage) > float(coords_dict.get(contig_name)[0])):
                coords_dict.update({contig_name: [query_coverage, s1_data, chromosome_num]})

        line_no += 1

    return coords_dict


def analyze_represented_contigs(linkage_map_file, coords_dict):
    linkage_contig_names = []
    query_coverage_data = []

    # Open linkage map file and reformat
    with open(linkage_map_file, newline='') as file:
        reader = csv.reader(file)
        linkage_map = list(reader)

    # Fix formatting issue with "_" in linkage map data
    for entry in linkage_map:
        contig_name = entry[1].split('_')[0]
        linkage_contig_names.append(contig_name)

    for name in linkage_contig_names:
        query_coverage, starting_position, chromosome_num = coords_dict.get(name)

        # Update query coverage data
        query_coverage_data.append(float(query_coverage))

    plt.hist(query_coverage_data, 100)
    plt.xlim(xmin=0, xmax=20)
    # plt.ylim(ymin=0, ymax=50)
    plt.ylabel("Number of contigs in bin")  # Weird issue here, ask during call
    plt.xlabel("% Query match to reference (Rio))")
    plt.title(f'Best hit query coverage data for contigs represented in linkage map')
    plt.show()


coords_dict = coord_to_dict(coords_file)
print(coords_dict.get('tig00001270'))
analyze_represented_contigs(linkage_maps_file, coords_dict)
