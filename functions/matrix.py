import os
import csv
from functions.load_bed_file import load_bed_file
from functions.compare_circular_data import are_circles_equal

def matrix(circular_dir, tools, sample, output_dir):
    """
    Create a matrix.csv file showing presence/absence of circRNAs in different tools.

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
        compare_circles_file = os.path.join(tool_dir, sample, f"{sample}.{tool}.bed")
        
        if os.path.exists(compare_circles_file):
            # Load circle data from the comparing file
            compare_circles_data = load_bed_file(compare_circles_file)
            
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

    # Convert matrix_dict to a list of rows for CSV writing
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
