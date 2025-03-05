# load_bed_file.py
import os

def load_bed_file(file_path):
    """
    Load a bed file and return the data as a list of tuples.

    Parameters:
        file_path (str): Path to the bed file.

    Returns:
        list: List of tuples representing circles.
    """
    circle_data = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chr_name = fields[0].replace('chr', '')  # Remove 'chr' from the chromosome name
                start = int(fields[1])
                end = int(fields[2])
                circle_data.append((chr_name, start, end))
    return circle_data
