# compare_circular_data.py
import os
import csv
from functions.load_bed_file import load_bed_file

def compare_bed_files(file1, file2, threshold=20):
    """
    Compare two bed files and find common circles.

    Parameters:
        file1 (str): Path to the first bed file.
        file2 (str): Path to the second bed file.
        threshold (int): Maximum allowed distance between start and end positions (default is 20).

    Returns:
        list: List of common circles.
    """
    circle_data1 = load_bed_file(file1)
    circle_data2 = load_bed_file(file2)
    common_circles = []
    matched_circles = set()
    for circle1 in circle_data1:
        for circle2 in circle_data2:
            if are_circles_equal(circle1, circle2, threshold) and circle2 not in matched_circles:
                common_circles.append(circle1)
                matched_circles.add(circle2)
                break
    return common_circles

def are_circles_equal(circle1, circle2, threshold=20):
    """
    Check if two circles overlap within a specified threshold.

    Parameters:
        circle1 (tuple): Tuple representing the first circle.
        circle2 (tuple): Tuple representing the second circle.
        threshold (int): Maximum allowed distance between start and end positions (default is 20).

    Returns:
        bool: True if the circles overlap within the threshold, False otherwise.
    """
    chr_name1, start1, end1 = circle1
    chr_name2, start2, end2 = circle2
    if chr_name1 != chr_name2:
        return False
    return abs(start1 - start2) <= threshold and abs(end1 - end2) <= threshold

def save_common_circles_to_file(common_circles, output_file):
    """
    Save common circles to a bed file.

    Parameters:
        common_circles (list): List of common circles.
        output_file (str): Path to the output bed file.
    """
    with open(output_file, 'w') as file:
        for circle in common_circles:
            chr_name, start, end = circle
            file.write(f"{chr_name}\t{start}\t{end}\n")

def save_circles_by_tool_and_coverage(circles, output_directory, tool, coverage):
    """
    Save circles to individual files based on tool and coverage.

    Parameters:
        circles (list): List of circles.
        output_directory (str): Path to the output directory.
        tool (str): Name of the tool.
        coverage (int): Coverage value.

    Returns:
        None
    """
    tool_dir = os.path.join(output_directory, tool)
    if not os.path.exists(tool_dir):
        os.makedirs(tool_dir)
    output_file = os.path.join(tool_dir, f"cov{coverage}_{tool}.bed")
    with open(output_file, 'w') as file:
        for circle in circles:
            chr_name, start, end = circle
            file.write(f"{chr_name}\t{start}\t{end}\n")

def organize_results(true_positives, false_negatives, false_positives, output_directory):
    """
    Organize results into folders and save circles.

    Parameters:
        true_positives (dict): Dictionary containing true positive circles.
        false_negatives (dict): Dictionary containing false negative circles.
        false_positives (dict): Dictionary containing false positive circles.
        output_directory (str): Path to the output directory.

    Returns:
        None
    """
    for tool in true_positives:
        tool_dir = os.path.join(output_directory)
        if not os.path.exists(tool_dir):
            os.makedirs(tool_dir)
        true_positives_tool_dir = os.path.join(tool_dir, "truepositives")
        if not os.path.exists(true_positives_tool_dir):
            os.makedirs(true_positives_tool_dir)
        for coverage in true_positives[tool]:
            save_circles_by_tool_and_coverage(true_positives[tool][coverage], true_positives_tool_dir, tool, coverage)
        false_negatives_tool_dir = os.path.join(tool_dir, "falsenegatives")
        if not os.path.exists(false_negatives_tool_dir):
            os.makedirs(false_negatives_tool_dir)
        for coverage in false_negatives[tool]:
            save_circles_by_tool_and_coverage(false_negatives[tool][coverage], false_negatives_tool_dir, tool, coverage)
        false_positives_tool_dir = os.path.join(tool_dir, "falsepositives")
        if not os.path.exists(false_positives_tool_dir):
            os.makedirs(false_positives_tool_dir)
        for coverage in false_positives[tool]:
            save_circles_by_tool_and_coverage(false_positives[tool][coverage], false_positives_tool_dir, tool, coverage)

def calculate_metrics(true_positives, false_positives, false_negatives):
    """
    Calculate precision, recall, and F-score.

    Parameters:
        true_positives (int): Number of true positives.
        false_positives (int): Number of false positives.
        false_negatives (int): Number of false negatives.

    Returns:
        float: Precision value.
        float: Recall value.
        float: F-score value.
    """
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) != 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) != 0 else 0
    fscore = (2 * precision * recall) / (precision + recall) if (precision + recall) != 0 else 0
    return precision, recall, fscore

def write_csv_file(data, header, output_file):
    """
    Write data to a CSV file.

    Parameters:
        data (dict): Dictionary containing data rows.
        header (list): List of column headers.
        output_file (str): Path to the output CSV file.

    Returns:
        None
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for key in data:
            writer.writerow([key] + data[key])

def write_metrics_to_csv(metrics, output_file, coverages):
    """
    Write precision, recall, and F-score metrics to a CSV file.

    Parameters:
        metrics (dict): Dictionary containing the metrics.
        output_file (str): Path to the output CSV file.

    Returns:
        None
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = ['Tool'] + coverages
        writer.writerow(header)
        for tool, values in metrics.items():
            writer.writerow([tool] + values)

def analyze_circular_data(circular_dir, circular_bed, output_dir, tools, coverages, threshold=20):
    """
    Analyze circular data across multiple tools and coverages.

    Parameters:
        circular_dir (str): Directory containing circular detection results for each tool.
        circular_bed (str): Path to the known circular BED file (ground truth).
        output_dir (str): Directory where results and statistics will be saved.
        tools (list of str): List of tool names to analyze.
        coverages (list of int): List of coverage levels to evaluate.
        threshold (int, optional): Maximum allowed distance between start and end positions (default is 20).

    Returns:
        None
    """
    common_circles_count = {}
    false_positives_count = {}
    false_negatives_count = {}
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Analyze circular data for each tool and coverage
    for tool in tools:
        tool_dir = os.path.join(circular_dir, tool)
        common_circles_count[tool] = {}
        false_positives_count[tool] = {}
        false_negatives_count[tool] = {}
        
        for coverage in coverages:
            filename = f"cov{coverage}_{tool}.bed"
            filepath = os.path.join(tool_dir, filename)
            common_circles = compare_bed_files(circular_bed, filepath, threshold)
            common_circles_count[tool][coverage] = len(common_circles)
            
            total_detected = len(load_bed_file(filepath))
            false_positives_count[tool][coverage] = total_detected - len(common_circles)
            
            total_known_circles = len(load_bed_file(circular_bed))
            false_negatives_count[tool][coverage] = total_known_circles - len(common_circles)
    
    # Create the "statistics" directory if it does not exist
    stats_directory = os.path.join(output_dir, 'statistics')
    if not os.path.exists(stats_directory):
        os.makedirs(stats_directory)

    # Write true positives CSV file
    true_positives_data = {}
    for tool in tools:
        true_positives_data[tool] = [str(common_circles_count[tool][cov]) for cov in coverages]
    write_csv_file(true_positives_data, ['Tool'] + [f"cov{cov}" for cov in coverages], os.path.join(stats_directory, 'truepositives.csv'))
    
    # Write false positives CSV file
    false_positives_data = {}
    for tool in tools:
        false_positives_data[tool] = [str(false_positives_count[tool][cov]) for cov in coverages]
    write_csv_file(false_positives_data, ['Tool'] + [f"cov{cov}" for cov in coverages], os.path.join(stats_directory, 'falsepositives.csv'))
    
    # Write false negatives CSV file
    false_negatives_data = {}
    for tool in tools:
        false_negatives_data[tool] = [str(false_negatives_count[tool][cov]) for cov in coverages]
    write_csv_file(false_negatives_data, ['Tool'] + [f"cov{cov}" for cov in coverages], os.path.join(stats_directory, 'falsenegatives.csv'))
    
    true_positives = {}
    false_negatives = {}
    false_positives = {}

    # Loop through each tool and coverage again
    for tool in tools:
        true_positives[tool] = {}
        false_negatives[tool] = {}
        false_positives[tool] = {}
        tool_dir = os.path.join(circular_dir, tool)
        for coverage in coverages:
            filename = f"cov{coverage}_{tool}.bed"
            filepath = os.path.join(tool_dir, filename)
            common_circles = compare_bed_files(filepath, circular_bed, threshold)
            true_positives[tool][coverage] = common_circles
            
            total_detected = load_bed_file(filepath)
            false_positives[tool][coverage] = [circle for circle in total_detected if circle not in common_circles]

            common_circles2 = compare_bed_files(circular_bed, filepath, threshold)
            total_known_circles = load_bed_file(circular_bed)
            false_negatives[tool][coverage] = [circle for circle in total_known_circles if circle not in common_circles2]
    
    organize_results(true_positives, false_negatives, false_positives, output_dir)

    # Create the header for the matrix file
    header = ['Circles']
    header.extend(tools)  # Add tool names as column headers

    # Load known circles from the provided bed file
    known_circles = load_bed_file(circular_bed)

    for coverage in coverages:
        matrix = []  # Reset the matrix for each coverage
        header = ['Circle'] + tools  # Create the header row for the matrix

        for known_circle in known_circles:
            row = [f"{known_circle[0]}:{known_circle[1]}-{known_circle[2]}"]  # Format circle information
            for tool in tools:
                tool_dir = os.path.join(circular_dir, tool)
                compare_circles_file = os.path.join(tool_dir, f"cov{coverage}_{tool}.bed")
                if os.path.exists(compare_circles_file):
                    # Load circle data from the comparing file
                    compare_circles_data = load_bed_file(compare_circles_file)
                    # Check if the known circle is in common with any circle in the comparing file for the current tool
                    found = False
                    for compare_circle in compare_circles_data:
                        if are_circles_equal(known_circle, compare_circle):
                            found = True
                            break
                    if found:
                        row.append(1)
                    else:
                        row.append(0)
                else:
                    row.append(0)  # If the comparing file doesn't exist, put 0
            matrix.append(row)

        # Write the matrix to a file inside the statistics directory for the current coverage
        matrix_file = os.path.join(stats_directory, f'matrix_cov{coverage}.csv')
        with open(matrix_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)  # Write the header row
            for row in matrix:
                writer.writerow(row)

    # Initialize dictionaries to store precision, recall, and F-score for each tool
    precision_by_tool = {tool: [] for tool in tools}
    recall_by_tool = {tool: [] for tool in tools}
    fscore_by_tool = {tool: [] for tool in tools}

    # Iterate over tools
    for tool in tools:
        # Initialize lists to store precision, recall, and F-score for each coverage
        precision_values = []
        recall_values = []
        fscore_values = []

        # Iterate over coverages
        for coverage in coverages:
            # Retrieve the counts for false positives and false negatives for the current tool and coverage
            false_positives = false_positives_count[tool][coverage]
            false_negatives = false_negatives_count[tool][coverage]

            # Calculate true positives (already available in common_circles_count)
            true_positives = common_circles_count[tool][coverage]

            # Calculate precision, recall, and F-score
            precision, recall, fscore = calculate_metrics(true_positives, false_positives, false_negatives)

            # Append metrics to lists
            precision_values.append(precision)
            recall_values.append(recall)
            fscore_values.append(fscore)

        # Store aggregated metrics for the current tool
        precision_by_tool[tool] = precision_values
        recall_by_tool[tool] = recall_values
        fscore_by_tool[tool] = fscore_values

    # Write aggregated precision, recall, and F-score metrics to CSV files
    precision_output_file = os.path.join(stats_directory, 'precision.csv')
    recall_output_file = os.path.join(stats_directory, 'recall.csv')
    fscore_output_file = os.path.join(stats_directory, 'fscore.csv')

    write_metrics_to_csv(precision_by_tool, precision_output_file, coverages)
    write_metrics_to_csv(recall_by_tool, recall_output_file, coverages)
    write_metrics_to_csv(fscore_by_tool, fscore_output_file, coverages)