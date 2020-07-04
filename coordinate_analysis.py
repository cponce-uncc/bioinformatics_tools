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


def get_coords_info(file_path, return_names=True, return_hits=True):
    file = open(file_path, 'rt')
    line_no = 0
    contig_names = []
    all_query_hits = []
    best_hits = {}

    for line in file:  # Iterate through every line in the file

        if line_no > 2:
            separated_contents = line.split('|')  # Split contents using the | character as the delimiter

            # Get query coverage for each hit
            query_coverage = separated_contents[5].split('   ')[2].strip()

            # Get contig name by taking seventh column, splitting it at the tab character, remove newline, and strip
            contig_name = separated_contents[6].split('\t')[1] \
                .replace('\n', ' ') \
                .strip()

            contig_names.append(str(contig_name))  # Add contig name to list
            all_query_hits.append([contig_name, query_coverage])

        line_no += 1  # Iterate line counter

    for entry in all_query_hits:

        if entry[0] in best_hits.keys():

            if best_hits.get(entry[0]) < entry[1]:  # Update dictionary if new entry is greater
                best_hits.update({entry[0]: entry[1]})

        else:  # Put entry into dictionary if contig has no entry yet
            best_hits.update({entry[0]: entry[1]})

    # Return data based on optional arguments passed to function
    if (return_names is True) and (return_hits is False):
        return contig_names

    elif (return_names is False) and (return_hits is True):
        return best_hits

    elif (return_names is True) and (return_hits is True):
        return contig_names, best_hits


def create_csv(assembly_file_path, coord_file_path, new_file_name):
    assembly_names = get_assembly_contigs(assembly_file_path)
    best_hits = get_coords_info(coord_file_path, return_names=False)

    with open(new_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Contig Name", "Represented in .coords", "Best Query Hit"])

        for entry in assembly_names:
            if entry in best_hits.keys():
                writer.writerow([entry, "Yes", best_hits.get(entry)])
            else:
                writer.writerow([entry, "No", "N/A"])

        file.close()


create_csv('ChineseAmber.contigs.fasta', 'ChineseAmber.coords', 'ChineseAmber_contig_analysis.csv')
