"""Translate linkage map data into list format (contig id, coordinate info, chromosome)
Initialize empty dictionary to hold data for linkage map contigs (contig_id, query coverage, s1_data, chromosome_num)

For entry in linkage map data:
- Iterate through each hit on the correct chromosome in the .coords file

     if the coordinate data falls in the correct range of the hit:

               if the hit that matches the coordinate data is of the highest Q coverage or not in the dict:

                              Add hit to dictionary

Overwrite the entries for contigs that are represented in the linkage map in the current sorted contig csv based on the
new data."""
# Eliminate blank data in coordinate information
import csv


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


# Convert .coords to csv
# Translate linkage map data into list format
linkage_map_filepath = "ChineseAmber_Linkage_Map.csv"
sorted_contigs_filepath = "sorted_test.csv"
coord_filepath = "ChineseAmber.coords"

# Translate linkage map into list format
linkage_map_data = []
linkage_map_file = open(linkage_map_filepath, 'r')

for line in linkage_map_file:
    line_data = line.split(',')
    line_data[2] = line_data[2].replace('\n', '')
    linkage_map_data.append(line_data)

# Translate .coord file into dict format (key=contig_id, values=[query_coverage,s1,s2,e2,chromosome_num])
coord_data = {}
coord_file = open(coord_filepath, 'r')
line_counter = 0

for line in coord_file:
    if line_counter > 2:  # Avoid headers
        separated_contents = line.split('|')

        # Extract contig name for each line
        contig_name = separated_contents[6].split('\t')[1] \
            .replace('\n', ' ') \
            .strip()

        # Get query coverage for each hit
        query_coverage = separated_contents[5].split('   ')[2].strip()

        # Get S1data from .coords file
        s1_data = eliminate_blank_data(separated_contents[0].split(' '))

        s2_data = eliminate_blank_data(separated_contents[1].split(' '))

        e2_data = eliminate_blank_data(separated_contents[1].split(' '), 2)

        # Get chromosome tag
        chromosome_num = separated_contents[6] \
            .split('\t')[0] \
            .split('_')[1]

        # Create entry if contig is new
        if contig_name not in coord_data.keys():
            coord_data.update({contig_name: ['']})
            current_data = coord_data.get(contig_name)
            current_data.append([query_coverage, s1_data, s2_data, e2_data, chromosome_num])
            coord_data.update({contig_name: current_data})
        else:
            if contig_name in coord_data.keys():
                current_data = coord_data.get(contig_name)
                current_data.append([query_coverage, s1_data, s2_data, e2_data, chromosome_num])
                coord_data.update({contig_name: current_data})

    line_counter += 1

z = 0
best_linkage_contigs = {}
for entry in linkage_map_data:
    contig_id, coordinate_data, chromosome_num = entry  # Unpack values

    coordinate_hits = coord_data.get(contig_id)

    for hit in coordinate_hits:
        if hit != '':
            if int(hit[4]) == int(chromosome_num):
                if int(coordinate_data) in range(int(hit[3]), int(hit[2])):
                    if contig_id not in best_linkage_contigs.keys():
                        best_linkage_contigs.update({contig_id: [hit[0], hit[1], hit[4]]})
                        z += 1
                    else:
                        if float(hit[0]) > float(best_linkage_contigs.get(contig_id)[0]):
                            best_linkage_contigs.update({contig_id: [hit[0], hit[1], hit[4]]})
                            z += 1

# Convert dictionary to list
best_linkage_contigs_list = []
for entry in best_linkage_contigs:
    query_coverage, s1_data, chromosome_num = best_linkage_contigs.get(entry)
    best_linkage_contigs_list.append([entry, query_coverage, s1_data, chromosome_num])

# Write data to csv file
with open('linkage_map_s1_data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for entry in best_linkage_contigs_list:
        writer.writerow(entry)
