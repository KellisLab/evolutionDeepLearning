import csv
from collections import defaultdict
from pathlib import Path

# Input and output paths
input_file = "../subsampledBiosamples.txt"
file_handles = {}

try:
    # Open the input file
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')  # Use tab delimiter
        
        for row in reader:
            group = row['Group']  # Access the 'Group' column
            
            # Open a file for this group if not already open
            if group not in file_handles:
                output_file = f"{group}.txt"
                file_handles[group] = open(output_file, 'w', newline='')
                
                # Write the header row to the group file
                writer = csv.DictWriter(file_handles[group], fieldnames=reader.fieldnames, delimiter='\t')
                writer.writeheader()
            
            # Write the current row to the appropriate file
            writer = csv.DictWriter(file_handles[group], fieldnames=reader.fieldnames, delimiter='\t')
            writer.writerow(row)

finally:
    # Close all open file handles
    for handle in file_handles.values():
        handle.close()
