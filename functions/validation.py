import pandas as pd
from functions.compare_circular_data import are_circles_equal

def load_bed(filepath):
    """
    Loads a BED file into a list of circles.

    Parameters:
        filepath (str): Path to the BED file.

    Returns:
        list of tuples: Each tuple represents a circle (chr, start, end).
    """
    df = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'start', 'end'])
    circles = list(df.itertuples(index=False, name=None))
    return circles

def format_circle(circle):
    """
    Formats a circle tuple into a string 'chr:start-end'.

    Parameters:
        circle (tuple): A tuple representing a circle (chr, start, end).

    Returns:
        str: Formatted circle string.
    """
    return f"{circle[0]}:{circle[1]}-{circle[2]}"

def build_tool_column(reference_circles, comparison_circles, threshold=20):
    """
    Creates a binary list indicating for each reference circle whether it matches
    any circle in the comparison set, within a given threshold.

    Parameters:
        reference_circles (list of tuples): Reference circles.
        comparison_circles (list of tuples): Circles to compare against.
        threshold (int, optional): Maximum distance for circles to be considered equal. Default is 20.

    Returns:
        list of int: 1 if matched, 0 if not, for each reference circle.
    """
    col = []
    for ref_circle in reference_circles:
        # Check if ANY circle in comparison_circles matches this ref_circle
        matched = any(are_circles_equal(ref_circle, comp_circle, threshold) for comp_circle in comparison_circles)
        col.append(1 if matched else 0)
    return col

def matrix_fp(circular_bed, tools, circular_output_path, circle_type, threshold=20):
    """
    Generates a matrix indicating the presence (1) or absence (0) of reference circles
    across multiple filtering folders and tools, then saves the results as CSV.

    Parameters:
        circular_bed (str): Path to the reference circles BED file.
        tools (list of str): List of tool names.
        circular_output_path (str): Output CSV file path.
        circle_type (str): Type of circle, used to build paths.
        threshold (int, optional): Threshold for circle matching. Default is 20.

    Returns:
        None
    """
    # Define the folders on disk vs. the column names you want
    filter_folders   = ['unfilter', 'filter-split',   'filter-duplicates', 'filter']

    # Load reference circles
    reference_circles = load_bed(circular_bed)
    circle_labels     = [format_circle(c) for c in reference_circles]

    # Start DataFrame with Circles + Simulated
    circular_output = pd.DataFrame({
        'Circles':  circle_labels,
        'Simulated': [1] * len(reference_circles)
    })

    # For each filtering folder, mark 1 if circle is in any tool under that folder
    for folder in filter_folders:
        presence = []
        for ref in reference_circles:
            found = False
            for tool in tools:
                fp = (
                    f'/data/benchmarking/results/'
                    f'{circle_type}/insilico/{folder}/truepositives/'
                    f'{tool}/cov30_{tool}.bed'
                )
                comp = load_bed(fp)
                if any(are_circles_equal(ref, c, threshold) for c in comp):
                    found = True
                    break
            presence.append(1 if found else 0)
        circular_output[folder] = presence

    # Save with exactly the header you asked for
    circular_output.to_csv(circular_output_path, index=False)

