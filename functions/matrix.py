import os
import csv
import pandas as pd
from functions.load_bed_file import load_bed_file
from functions.compare_circular_data import are_circles_equal

def matrix_real(circular_dir, tools, sample, output_dir):
    """
    Create a matrix.csv file showing presence/absence of circles in different tools.

    Parameters:
        circular_dir: The base directory where tool folders are located.
        tools: List of tool names.
        sample: The sample name (e.g., "HCMS01").
        output_dir: The directory where the output matrix.csv will be saved.

    Returns:
        list: List of common circles.
    """
    # Create the header for the matrix file
    header = ['Circles'] + tools

    # Initialize a dictionary to store the presence/absence values with circles as keys
    matrix_dict = {}

    # Iterate through each tool to collect all unique circles
    for tool in tools:
        tool_dir = os.path.join(circular_dir, tool)
        compare_circles_file = os.path.join(tool_dir, f"{sample}.{tool}.bed")
        
        if os.path.exists(compare_circles_file):
            # Load circle data from the comparing file
            compare_circles_data = load_bed_file(compare_circles_file)
            
            if not compare_circles_data:
                print(f"No circles found in file: {compare_circles_file}")  # Debugging line

            for compare_circle in compare_circles_data:
                circle_key = f"{compare_circle[0]}:{compare_circle[1]}-{compare_circle[2]}"

                # Check if this circle (or an equivalent one) is already in the matrix_dict
                found = False
                for existing_circle_key in matrix_dict:
                    existing_circle = existing_circle_key.split(':')[1].split('-')
                    existing_circle_start = int(existing_circle[0])
                    existing_circle_end = int(existing_circle[1])

                    if are_circles_equal(compare_circle, [compare_circle[0], existing_circle_start, existing_circle_end]):
                        circle_key = existing_circle_key
                        found = True

                        # If circles are equal, keep the longer one
                        if existing_circle_end - existing_circle_start < compare_circle[2] - compare_circle[1]:
                            matrix_dict[circle_key] = matrix_dict.pop(existing_circle_key)
                        break

                # If this is a new circle, initialize its presence/absence list
                if not found:
                    matrix_dict[circle_key] = [0] * len(tools)

                # Mark this circle as present for the current tool
                tool_index = tools.index(tool)
                matrix_dict[circle_key][tool_index] = 1
        else:
            print(f"File does not exist: {compare_circles_file}")  # Debugging line

    # Convert matrix_dict to a list of rows for CSV writing
    if not matrix_dict:
        print("No data to write in the matrix.")  # Debugging line
    
    matrix = [[circle] + matrix_dict[circle] for circle in matrix_dict]
  
    # Create the directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write the matrix to a file inside the output directory
    matrix_file = os.path.join(output_dir, 'matrix.csv')
    with open(matrix_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)  # Write the header row
        writer.writerows(matrix)  # Write all matrix rows

def create_circle_presence_matrix(circular_bed, tools, circular_output, circular_type):
    """"
    Create a matrix indicating the presence of circles across different filtering types and tools.

    Parameters:
        circular_bed (str): Path to the simulated circles BED file.
        tools (list): List of tool names.
        circular_output (str): Output CSV file path.
        circular_type (str): Type of circle, used to build paths.

    Returns:
        None
    """
    def load_bed(file_path):
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return set()
        bed = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2])
        return set(map(tuple, bed.values))

    # Simulated circles
    simulated_file = os.path.join(circular_bed)
    simulated = load_bed(simulated_file)

    all_circles = set(simulated)

    # Dictionary to store results
    matrix = {}

    # Store simulated data
    for circle in simulated:
        matrix[circle] = {'Simulated': 1}

    # Filtering types and their corresponding subpaths
    filter_types = {
        'unfilter': 'unfilter',
        'filter-split': 'filter-split',
        'filter-duplicates': 'filter-duplicates',
        'filter': 'filter'
    }

    # Loop through each tool and each filter type
    for tool in tools:
        for label, folder in filter_types.items():
            bed_path = os.path.join(
                f'results/{circular_type}/insilico/{folder}/falsepositives/{tool}/cov30_{tool}.bed'
            )
            found = load_bed(bed_path)
            all_circles.update(found)
            for circle in found:
                if circle not in matrix:
                    matrix[circle] = {}
                matrix[circle][label] = 1

        # Fill missing columns with 0
    df = pd.DataFrame.from_dict(matrix, orient='index').fillna(0).astype(int)

    # Add Circle column with coordinates as string
    df['Circles'] = df.index.map(lambda x: f"{x[0]}:{x[1]}-{x[2]}")
    df.reset_index(drop=True, inplace=True)

    # Ensure all expected columns exist
    expected_cols = ['Simulated', 'unfilter', 'filter-split', 'filter-duplicates', 'filter']
    for col in expected_cols:
        if col not in df.columns:
            df[col] = 0

    # Reorder columns: Circle + expected
    matrix = df[['Circles'] + expected_cols]

    matrix.to_csv(circular_output, index=False)
