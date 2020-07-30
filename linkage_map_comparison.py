# NOTE: I only found a single contig (tig00002499) at the 70% threshold level. 15% yielded better results.
import csv

"""I KNOW IT SAID 78"""

threshold_hits_file = 'fifteen_percent_threshold.csv'
linkage_maps_file = 'ChineseAmber_Linkage_Map.csv'


def convert_to_csv(old_file_path, new_file_path):
    old_file = open(old_file_path, "r")

    with open(new_file_path, 'w', newline='') as file:
        writer = csv.writer(file)

        for line in old_file:
            old_content = line.split(' ')
            old_content[1] = old_content[1].replace('S', 'tig')
            old_content[2] = old_content[2].replace('\n', '')
            writer.writerow(old_content)

    old_file.close()  # Close old .txt file
    file.close()  # Close new .csv file


threshold_hits = list(open(threshold_hits_file, 'r'))
linkage_map = list(open(linkage_maps_file, 'r'))

threshold_contigs = []
linkage_contigs = []

# Extract all contigs represented in linkage map
for line in linkage_map:
    line = line.replace("\n", '') \
            .split(",")[1] \
            .split("_")[0]
    linkage_contigs.append(line)
linkage_contigs = set(linkage_contigs)  # Cast to set to allow for use of intersection method

# Extract all contigs represented in threshold .coords file
for line in threshold_hits:
    line = line.replace("\n", '') \
        .split(",")[0] \

    if line not in threshold_contigs:  # Only add one of each contig to represented contig list
        threshold_contigs.append(line)
threshold_contigs.pop(0)  # Remove header from data
threshold_contigs = set(threshold_contigs)  # Cast to set to allow for use of intersection method

common_contigs = threshold_contigs.intersection(linkage_contigs)  # Find shared contigs

shared_threshold_contigs = threshold_contigs.intersection(common_contigs)  # Remove elements specific to set
shared_linkage_contigs = linkage_contigs.intersection(common_contigs)  # Remove elements specific to set

# Print results
print(f'There are {len(common_contigs)} contigs shared between the linkage map and the threshold file.')
print(shared_linkage_contigs)
print(shared_threshold_contigs)


# convert_to_csv('Chinese_Amber_Results_onemap.txt', 'ChineseAmber_Linkage_Map.csv')