# bioinformatics_tools
About Repository
This repository contains scripts written to assist in de novo assembly during Summer 2020 at the Cooper Lab of the North Carolina Research Campus. The scripts are separated into two sections: one focusing on troubleshooting and analysis of output files (MUMmer pairwise alignment results, linkage maps, contig assemblies), and the second section focusing on sorting contigs by their position on a chromosome.

Repository link: https://github.com/cponce-uncc/bioinformatics_tools

*Scripts Contained in Repository:*
Assembly_quality.py: Calculates basic statistics to analyze the quality of the contig assembly file, including n50/l50, n90/l90, total bp count, bp percentages, number of contigs, and the contig with the largest bp.
    
    Input: Contig assembly file, targeted_extensions.csv (Used by the assembly_quality.py script to check for extensions to run analysis on automatically. Add a targeted extension in the following format (extension, assembler name). Examples are provided in the file for Canu, SmartDenovo, and WTBDG.)

Contig_order_comparison.py: Compares the order of contigs and the total order similarity (using SequenceMatcher). This can be used to compare the order of linkage maps or ordered contig files (In the format used in this repository).
    
    Input: Linkage map file (.csv) or sorted contig file (.csv) that you wish to compare

Scatter_plot_generator.py: Generates a 3D scatter plot of nucleotide coverage, size in bp, and query coverage percentage.
    
    Input: .coords file from MUMmer pairwise alignment

Linkage_map_analysis.py: Checks best-hit query coverage data for contigs represented in the linkage map file and outputs histogram with query coverage percentage data for represented contigs.
    
    Input: Linkage map file (.txt)

Extract_s1_data.py: Performs a reverse search through the .coords file using the mutation coordinate provided by the linkage map file in order to obtain the positional data of contigs represented in the linkage map.
    
    Input: Linkage map file (.csv), .coords file from MUMmer pairwise alignment

Coordinate_data_processing.py: Finds the hit with the highest query coverage percentage and extracts the positional data and then sorts the contigs based on this data. 
    
    Input: .coords file from MUMmer pairwise alignment

Contig_chunking.py: Iterates through the .coords file and chunks sequential hits together to obtain more accurate positional data (Contigs that qualified for this method covered roughly ⅗ of the genome) and combines this new positional data with the sorted contigs from coordinate_data_processsing.py.
    
    Input: Best hit (QC%) .csv file from cooridnate_data_processing.py
    
coordinate_data_processing.py (.csv) and .coords file from MUMmer pairwise alignment
    
    Input: .coords file from MUMmer pairwise alignment
    
Write_sorted_contigs.py: Takes final sorted contig file from contig_chunking.py and iterates through the .contigs.fasta file to grab the data for each contig and write it to a new .fasta file.

    Input: Updated, sorted contig list from contig_chunking.py (.csv), .coords file, contig assembly file 



*Unsorted to Sorted Contigs Guide*

    NOTE: Steps must be done in sequential order.

1.) Place .coord file for line in the local repository and edit the variable “line_file” in coordinate_data_processing.py to get the “best hits” .csv file.

2.) Edit the variable “line_file” in contig_chunking.py and execute to script to obtain a new .csv file updated with the contig chunk positions.

3.) Edit the variable “line_file” in write_sorted contigs and the variable “contigs_file” and execute the script to get the .fasta file containing the written, sorted contigs


