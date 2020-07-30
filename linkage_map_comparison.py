# Get two lists
# Filter out
import csv


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


convert_to_csv('Chinese_Amber_Results_onemap.txt', 'ChineseAmber_Linkage_Map.csv')