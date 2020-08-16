# NOTE: I only found a single contig (tig00002499) at the 70% threshold level. 15% yielded better results.
import csv
import difflib

threshold_hits_file = 'sorted_test_plus_linkage.csv'
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
    list_one = set(list_one)
    list_two = set(list_two)
    return list_one.intersection(list_two)


threshold_hits = list(open(threshold_hits_file, 'r'))
linkage_map = list(open(linkage_maps_file, 'r'))

old_threshold_contigs = []
linkage_contigs = []

# Extract all contigs represented in linkage map
for line in linkage_map:
    line = line.replace("\n", '') \
            .split(",")[0]

    if line not in linkage_contigs:
        linkage_contigs.append(line)

# Extract all contigs represented in threshold .coords file
for line in threshold_hits:
    line = line.replace("\n", '') \
        .split(",")[0] \

    if line not in old_threshold_contigs:  # Only add one of each contig to represented contig list
        old_threshold_contigs.append(line)

old_threshold_contigs.pop(0)  # Remove header from data


new_threshold_contigs = []

for entry in old_threshold_contigs:
    if entry in linkage_contigs:
        new_threshold_contigs.append(entry)

with open('new_common_contig_order.csv', 'w', newline='') as file:
    x = 0
    y = 0
    writer = csv.writer(file)
    writer.writerow(['position', '.coords contigs', 'linkage contigs'])

    while x < len(new_threshold_contigs):
        writer.writerow([x + 1, new_threshold_contigs[x], linkage_contigs[x]])
        if new_threshold_contigs[x] == linkage_contigs[x]:
            y += 1
        x += 1

print(linkage_contigs)
print(new_threshold_contigs)
sequence_matcher = difflib.SequenceMatcher(None, linkage_contigs, new_threshold_contigs)
print(f"Similarity in order between updated sorted contigs\nand linkage map contigs: {sequence_matcher.ratio() * 100}%")


convert_to_csv('Chinese_Amber_results_onemap.txt', 'ChineseAmber_Linkage_Map.csv')
