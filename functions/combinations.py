import os
import itertools
import shutil
import csv
from math import isclose
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import median_abs_deviation, kruskal
import scikit_posthocs as sp
from functions.load_bed_file import load_bed_file
from functions.compare_circular_data import are_circles_equal, compare_bed_files

def find_union_circles(sets):
    """
    Merge multiple sets of circles into a single set.
    
    Parameters:
        sets (list of lists): List of sets of circles.

    Returns:
        list: Union set of circles, with each unique circle appearing only once.
    """
    union_set = []  # Use a list to maintain the order of insertion

    # Iterate over each set in the input sets
    for s in sets:
        # Iterate over each circle in the current set
        for circle in s:
            # Check if this circle (or an equivalent one) already exists in union_set
            if not any(are_circles_equal(circle, existing_circle) for existing_circle in union_set):
                # If no equivalent circle is found, add it to the union set
                union_set.append(circle)

    return union_set

def find_intersect_circles(sets):
    """
    Find circles that are common across multiple sets.

    Parameters:
        sets (list of lists): List of sets of circles.

    Returns:
        list: List of circles common across all sets, with each unique circle appearing only once.
    """
    # Start with the first set as the baseline
    common_circles = sets[0] if sets else []

    # Iterate through all other sets and find the intersection
    for s in sets[1:]:
        new_common = []
        for circle1 in common_circles:
            # Check if circle1 is equal to any circle in the current set `s`
            if any(are_circles_equal(circle1, circle2) for circle2 in s):
                # Add the circle to new_common only if it's not already present
                if not any(are_circles_equal(circle1, existing_circle) for existing_circle in new_common):
                    new_common.append(circle1)
        common_circles = new_common  # Update the common_circles with the new intersection

    return common_circles

def analyze_circular_combination(circular_dir, circular_bed, output_dir, tools, tool_abbreviations):
    """
    Analyze and generate combinations of circular RNA data from multiple tools with abbreviated names.

    Parameters:
        circular_dir (str): Directory containing circular detection results for each tool.
        circular_bed (str): Path to the known circular BED file (ground truth).
        output_dir (str): Path to the directory where output files will be stored.
        tools (list of str): List of tool names to analyze.
        tool_abbreviations (dict): Dictionary mapping tool names to abbreviations.

    Returns:
        None
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load known circles
    known_circles = load_bed_file(circular_bed)

    # Dictionary to store circles by tool
    circles_by_tool = {}

    # Load circles for each tool
    for tool in tools:
        tool_directory = os.path.join(circular_dir, tool)
        tool_circles = load_bed_file(os.path.join(tool_directory, f"cov30_{tool}.bed"))
        circles_by_tool[tool] = tool_circles

    # Process combinations of tools
    for r in range(1, len(tools) + 1):
        for combination in itertools.combinations(tools, r):
            # Use abbreviations for the file name
            key = '_'.join(tool_abbreviations[tool] for tool in combination)

            # Gather circles for the current combination
            current_sets = [circles_by_tool[tool] for tool in combination]

            # Create directories for intersect and union circles if they don't exist
            intersect_directory = os.path.join(output_dir, 'intersect')
            if not os.path.exists(intersect_directory):
                os.makedirs(intersect_directory)

            union_directory = os.path.join(output_dir, 'union')
            if not os.path.exists(union_directory):
                os.makedirs(union_directory)

            if len(combination) == 1:
                tool = combination[0]
                # Copy individual tool files for both intersect and union
                intersect_output_file = os.path.join(intersect_directory, f"{key}.bed")
                union_output_file = os.path.join(union_directory, f"{key}.bed")
                shutil.copyfile(os.path.join(circular_dir, tool, f"cov30_{tool}.bed"), intersect_output_file)
                shutil.copyfile(os.path.join(circular_dir, tool, f"cov30_{tool}.bed"), union_output_file)
            else:
                # Intersect circles
                intersect_circles = find_intersect_circles(current_sets)

                # Union circles
                union_circles = find_union_circles(current_sets)

                # Write intersect circles to file
                intersect_output_file = os.path.join(intersect_directory, f"{key}.bed")
                with open(intersect_output_file, 'w') as file:
                    for circle in intersect_circles:
                        file.write(f"{circle[0]}\t{circle[1]}\t{circle[2]}\n")

                # Write union circles to file
                union_output_file = os.path.join(union_directory, f"{key}.bed")
                with open(union_output_file, 'w') as file:
                    for circle in union_circles:
                        file.write(f"{circle[0]}\t{circle[1]}\t{circle[2]}\n")

def rosette(circular_bed, input_dir, output_dir, tool_abbreviations):
    """
    Generates rosette combinations of circular RNA data using provided tools.

    Parameters:
        circular_bed (str): Path to the known circular BED file (ground truth).
        input_dir (str): Directory containing the intersected results for each tool.
        output_dir (str): Directory where the rosette combinations will be stored.
        tool_abbreviations (dict): Mapping of tool combinations to their abbreviations.

    Returns:
        None
    """
    # Reverse the dictionary for quick lookup
    tools = list(tool_abbreviations.keys())

    # Dictionary to store circles by tool
    circles_by_tool = {}

    # Load circles for each tool
    for tool in tools:
        tool_file = os.path.join(input_dir, f"{tool}.bed")
        if not os.path.exists(tool_file):
            raise FileNotFoundError(f"File not found: {tool_file}")
        tool_circles = load_bed_file(tool_file)
        circles_by_tool[tool] = tool_circles

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process combinations of tools
    for r in range(1, len(tools) + 1):
        for combination in itertools.combinations(tools, r):
            # Use abbreviations for the file name
            key = '_'.join(tool_abbreviations[tool] for tool in combination)
            
            # Gather circles for the current combination
            current_sets = [circles_by_tool[tool] for tool in combination]

            if len(combination) == 1:
                tool = combination[0]
                # Copy individual tool files to the output directory
                intersect_output_file = os.path.join(output_dir, f"{key}.bed")
                shutil.copyfile(os.path.join(input_dir, f"{tool}.bed"), intersect_output_file)
            else:
                # Intersect circles
                intersect_circles = find_union_circles(current_sets)
                
                # Write intersect circles to the output directory
                intersect_output_file = os.path.join(output_dir, f"{key}.bed")
                with open(intersect_output_file, 'w') as file:
                    for circle in intersect_circles:
                        file.write(f"{circle[0]}\t{circle[1]}\t{circle[2]}\n")

# Combinaciones reales (usando guiones bajos en lugar de guiones medios)
eccDNA_real_combinations = [
    "CE", "CM", "CF", "EB", "EM", "CE_CM", "CE_CF", "CE_EB", "CE_EM", "CM_CF", 
    "CM_EB", "CM_EM", "CF_EB", "CF_EM", "EB_EM", "CE_CM_CF", "CE_CM_EB", 
    "CE_CM_EM", "CE_CF_EB", "CE_CF_EM", "CE_EB_EM", "CM_CF_EB", "CM_CF_EM", 
    "CM_EB_EM", "CF_EB_EM", "CE_CM_CF_EB", "CE_CM_CF_EM", "CE_CM_EB_EM", 
    "CE_CF_EB_EM", "CM_CF_EB_EM", "CE_CM_CF_EB_EM"
]

# Combinaciones adaptadas para Rosette
eccDNA_ros_combinations_mapping = {
    "CE_CM_CF": "CE_CM_CE_CF_CM_CF",
    "CE_CM_EB": "CE_CM_CE_EB_CM_EB",
    "CE_CM_EM": "CE_CM_CE_EM_CM_EM",
    "CE_CF_EB": "CE_CF_CE_EB_CF_EB",
    "CE_CF_EM": "CE_CF_CE_EM_CF_EM",
    "CE_EB_EM": "CE_EB_CE_EM_EB_EM",
    "CM_CF_EB": "CM_CF_CM_EB_CF_EB",
    "CM_CF_EM": "CM_CF_CM_EM_CF_EM",
    "CM_EB_EM": "CM_EB_CM_EM_EB_EM",
    "CF_EB_EM": "CF_EB_CF_EM_EB_EM",
    "CE_CM_CF_EB": "CE_CM_CE_CF_CE_EB_CM_CF_CM_EM_CF_EB",
    "CE_CM_CF_EM": "CE_CM_CE_CF_CE_EM_CM_CF_CM_EM_CF_EM",
    "CE_CM_EB_EM": "CE_CM_CE_EB_CE_EM_CM_EB_CM_EM_EB_EM",
    "CE_CF_EB_EM": "CE_CF_CE_EB_CE_EM_CF_EB_CF_EM_EB_EM",
    "CM_CF_EB_EM": "CM_CF_CM_EB_CM_EM_CF_EB_CF_EM_EB_EM",
    "CE_CM_CF_EB_EM": "CE_CM_CE_CF_CE_EB_CE_EM_CM_CF_CM_EB_CM_EM_CF_EB_CF_EM_EB_EM"
}

# Tool name to abbreviation mapping
eccDNA_tools_abbreviations = {
    'CIRCexplorer2': 'CE',
    'Circle-Map': 'CM',
    'Circle_finder': 'CF',
    'ecc_finder-bwa': 'EB',
    'ecc_finder-minimap2': 'EM',
    'segemehl': 'SE',
}

eccDNA_tool_abbreviations_r = {
    'CE_CM': 'CE_CM',
    'CE_CF': 'CE_CF',
    'CE_EB': 'CE_EB',
    'CE_EM': 'CE_EM',
    'CM_CF': 'CM_CF',
    'CM_EB': 'CM_EB',
    'CM_EM': 'CM_EM',
    'CF_EB': 'CF_EB',
    'CF_EM': 'CF_EM',
    'EB_EM': 'EB_EM',
}

combination = ['intersect', 'union', 'rosette']

def process_eccDNA_filtering(filtering, circular_bed, tools):
    # Define paths
    filtering_base_dir = f'data/insilico/eccDNA/{filtering}'
    output_base_dir = f'results/eccDNA/insilico/{filtering}'
    
    # Analyze circular combinations
    analyze_circular_combination(filtering_base_dir, circular_bed, output_base_dir, tools, eccDNA_tools_abbreviations)

    # Perform rosette operation
    circular_dir_ros = f'{output_base_dir}/intersect'
    output_dir_ros = f'{output_base_dir}/rosette'
    rosette(circular_bed, circular_dir_ros, output_dir_ros, eccDNA_tool_abbreviations_r)
    # Process combinations
    for combine in combination:
        combine_dir = f'{output_base_dir}/{combine}'
        calculate_statistics(combine_dir, circular_bed, output_base_dir, name=combine)

    # Paths for combination statistics
    union_path = f'{output_base_dir}/union_statistics.csv'
    rosette_path = f'{output_base_dir}/rosette_statistics.csv'
    intersect_path = f'{output_base_dir}/intersect_statistics.csv'
    output_path = f'{output_base_dir}/combination_stats.xlsx'

    # Generate combination stats
    generate_combination_stats(union_path, rosette_path, intersect_path, output_path, eccDNA_real_combinations, eccDNA_ros_combinations_mapping)

    # Move generated bed folders to results directory
    for folder_name in ['union', 'rosette', 'intersect']:
        src = os.path.join(filtering_base_dir, folder_name)
        dst = os.path.join(output_base_dir, folder_name)

        if os.path.exists(dst):
            shutil.rmtree(dst)  # Remove existing destination folder if it exists

        if os.path.exists(src):
            shutil.move(src, dst)

# Combinaciones reales (usando guiones bajos en lugar de guiones medios)
circRNA_real_combinations = [
    "CE", "CF", "CQ", "FC", "SE", "CE_CF", "CE_CQ", "CE_FC", "CE_SE", "CF_CQ", 
    "CF_FC", "CF_SE", "CQ_FC", "CQ_SE", "FC_SE", "CE_CF_CQ", "CE_CF_FC", 
    "CE_CF_SE", "CE_CQ_FC", "CE_CQ_SE", "CE_FC_SE", "CF_CQ_FC", "CF_CQ_SE", 
    "CF_FC_SE", "CQ_FC_SE", "CE_CF_CQ_FC", "CE_CF_CQ_SE", "CE_CF_FC_SE", 
    "CE_CQ_FC_SE", "CF_CQ_FC_SE", "CE_CF_CQ_FC_SE"
]

# Combinaciones adaptadas para Rosette
circRNA_ros_combinations_mapping = {
    "CE_CF_CQ": "CE_CF_CE_CQ_CF_CQ",
    "CE_CF_FC": "CE_CF_CE_FC_CF_FC",
    "CE_CF_SE": "CE_CF_CE_SE_CF_SE",
    "CE_CQ_FC": "CE_CQ_CE_FC_CQ_FC",
    "CE_CQ_SE": "CE_CQ_CE_SE_CQ_SE",
    "CE_FC_SE": "CE_FC_CE_SE_FC_SE",
    "CF_CQ_FC": "CF_CQ_CF_FC_CQ_FC",
    "CF_CQ_SE": "CF_CQ_CF_SE_CQ_SE",
    "CF_FC_SE": "CF_FC_CF_SE_FC_SE",
    "CQ_FC_SE": "CQ_FC_CQ_SE_FC_SE",
    "CE_CF_CQ_FC": "CE_CF_CE_CQ_CE_FC_CF_CQ_CF_SE_CQ_FC",
    "CE_CF_CQ_SE": "CE_CF_CE_CQ_CE_SE_CF_CQ_CF_SE_CQ_SE",
    "CE_CF_FC_SE": "CE_CF_CE_FC_CE_SE_CF_FC_CF_SE_FC_SE",
    "CE_CQ_FC_SE": "CE_CQ_CE_FC_CE_SE_CQ_FC_CQ_SE_FC_SE",
    "CF_CQ_FC_SE": "CF_CQ_CF_FC_CF_SE_CQ_FC_CQ_SE_FC_SE",
    "CE_CF_CQ_FC_SE": "CE_CF_CE_CQ_CE_FC_CE_SE_CF_CQ_CF_FC_CF_SE_CQ_FC_CQ_SE_FC_SE"
}
# Tool name to abbreviation mapping
circRNA_tools_abbreviations = {
    'CIRCexplorer2': 'CE',
    'circRNA_finder': 'CF',
    'CIRIquant': 'CQ',
    'find_circ': 'FC',
    'segemehl': 'SE',
} 

circRNA_tool_abbreviations_r = {
    'CE_CF': 'CE_CF',
    'CE_CQ': 'CE_CQ',
    'CE_FC': 'CE_FC',
    'CE_SE': 'CE_SE',
    'CF_CQ': 'CF_CQ',
    'CF_FC': 'CF_FC',
    'CF_SE': 'CF_SE',
    'CQ_FC': 'CQ_FC',
    'CQ_SE': 'CQ_SE',
    'FC_SE': 'FC_SE',
}

def process_circRNA_filtering(filtering, circular_bed, tools):
    # Define paths
    filtering_base_dir = f'data/insilico/circRNA/{filtering}'
    output_base_dir = f'results/circRNA/insilico/{filtering}'
    
    # Analyze circular combinations
    analyze_circular_combination(filtering_base_dir, circular_bed, output_base_dir, tools, circRNA_tools_abbreviations)

    # Perform rosette operation
    circular_dir_ros = f'{output_base_dir}/intersect'
    output_dir_ros = f'{output_base_dir}/rosette'
    rosette(circular_bed, circular_dir_ros, output_dir_ros, circRNA_tool_abbreviations_r)
    # Process combinations
    for combine in combination:
        combine_dir = f'{output_base_dir}/{combine}'
        calculate_statistics(combine_dir, circular_bed, output_base_dir, name=combine)

    # Paths for combination statistics
    union_path = f'{output_base_dir}/union_statistics.csv'
    rosette_path = f'{output_base_dir}/rosette_statistics.csv'
    intersect_path = f'{output_base_dir}/intersect_statistics.csv'
    output_path = f'{output_base_dir}/combination_stats.xlsx'

    # Generate combination stats
    generate_combination_stats(union_path, rosette_path, intersect_path, output_path, circRNA_real_combinations, circRNA_ros_combinations_mapping)

    # Move generated bed folders to results directory
    for folder_name in ['union', 'rosette', 'intersect']:
        src = os.path.join(filtering_base_dir, folder_name)
        dst = os.path.join(output_base_dir, folder_name)

        if os.path.exists(dst):
            shutil.rmtree(dst)  # Remove existing destination folder if it exists

        if os.path.exists(src):
            shutil.move(src, dst)

def write_statistics_to_csv(statistics, output_dir):
    """
    Write statistics to a CSV file.

    Parameters:
        statistics (dict): Dictionary containing statistics for each tool.
        output_dir (str): Path to the output CSV file.

    Returns:
        None
    """
    with open(output_dir, 'w', newline='') as csvfile:
        fieldnames = ['Tool', 'Common Circles', 'False Positives', 'False Negatives', 'Precision', 'Recall', 'F-score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for tool, stats in statistics.items():
            writer.writerow({
                'Tool': tool,
                'Common Circles': stats.get('Common Circles', 0),
                'False Positives': stats.get('False Positives', 0),
                'False Negatives': stats.get('False Negatives', 0),
                'Precision': stats.get('Precision', 0),
                'Recall': stats.get('Recall', 0),
                'F-score': stats.get('F-score', 0)
            })

def calculate_statistics(circular_dir, circular_bed, output_dir, name):
    """
    Calculate statistics by comparing circles in circular_directory with known_circle_bed_file.

    Parameters:
        circular_dir (str): Directory containing circular detection results for each tool.
        circular_bed (str): Path to the known circular BED file (ground truth).
        output_dir (str): Path to the directory where output files will be stored.

    Return:
        None.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Dynamically list the BED files in the circular_directory
    files_to_analyze = [f for f in os.listdir(circular_dir) if f.endswith('.bed')]
    
    statistics = {}
    total_known_circles = len(load_bed_file(circular_bed))

    for filename in files_to_analyze:
        file_path = os.path.join(circular_dir, filename)
        
        if not os.path.exists(file_path):
            continue
        
        tool_name = os.path.splitext(filename)[0]
        statistics[tool_name] = {}
        
        common_circles = compare_bed_files(circular_bed, file_path)
        total_detected = len(load_bed_file(file_path))
        false_positives = total_detected - len(common_circles)
        false_negatives = total_known_circles - len(common_circles)
        
        statistics[tool_name]['Common Circles'] = len(common_circles)
        statistics[tool_name]['False Positives'] = false_positives
        statistics[tool_name]['False Negatives'] = false_negatives

        # Calculate Precision, Recall, and F-score
        if len(common_circles) > 0:
            precision = len(common_circles) / total_detected
            recall = len(common_circles) / total_known_circles
            fscore = (2 * precision * recall) / (precision + recall) if not isclose(precision + recall, 0, abs_tol=1e-9) else 0
        else:
            precision = 0
            recall = 0
            fscore = 0

        statistics[tool_name]['Precision'] = precision
        statistics[tool_name]['Recall'] = recall
        statistics[tool_name]['F-score'] = fscore

    # Write all statistics to CSV once at the end
    csv_file_path = os.path.join(output_dir, f'{name}_statistics.csv')
    write_statistics_to_csv(statistics, csv_file_path)


import pandas as pd

def generate_combination_stats(union_path, rosette_path, intersect_path, output_path, real_combinations, ros_combinations_mapping):
    """
    Calculates metrics (Common Circles, False Positives, False Negatives, Precision, Recall, F-score)
    from union, rosette, and intersect circle detection results and saves them into an Excel file.

    Parameters:
        union_path (str): Path to the CSV file containing union results per combination.
        rosette_path (str): Path to the CSV file containing rosette results per combination.
        intersect_path (str): Path to the CSV file containing intersection results per combination.
        output_path (str): Path to the Excel file where the results will be written.
        real_combinations (list): List of actual tool combinations to evaluate (e.g., ['CMF', 'CME']).
        ros_combinations_mapping (dict): Dictionary mapping combinations to their corresponding rosette name (used for grouping).

    Returns:
        None
    """

    # Load input data
    union_df = pd.read_csv(union_path)
    rosette_df = pd.read_csv(rosette_path)
    intersect_df = pd.read_csv(intersect_path)

    metrics_dict = {}  # Optional: to store metrics in memory

    with pd.ExcelWriter(output_path) as writer: 
        # Calculate Common Circles, False Positives, False Negatives
        for metric in ['Common Circles', 'False Positives', 'False Negatives']:
            temp_df = pd.DataFrame()
            temp_df['Combination'] = real_combinations

            # Map to rosette names if applicable
            temp_df['Rosette_names'] = temp_df['Combination'].map(ros_combinations_mapping).fillna(temp_df['Combination'])

            # Merge Union metric
            temp_df = temp_df.merge(
                union_df[['Tool', metric]],
                left_on='Combination',
                right_on='Tool',
                how='left'
            ).rename(columns={metric: 'Union'}).drop(columns=['Tool'])

            # Merge Rosette metric
            temp_df = temp_df.merge(
                rosette_df[['Tool', metric]],
                left_on='Rosette_names',
                right_on='Tool',
                how='left'
            ).rename(columns={metric: 'Rosette'}).drop(columns=['Tool'])

            # For specific tools, set Rosette = Union (special cases)
            temp_df.loc[temp_df['Combination'].isin(['CE', 'CM', 'CF', 'EB', 'EM']), 'Rosette'] = temp_df['Union']

            # Merge Intersect metric
            temp_df = temp_df.merge(
                intersect_df[['Tool', metric]],
                left_on='Combination',
                right_on='Tool',
                how='left'
            ).rename(columns={metric: 'Intersect'}).drop(columns=['Tool'])

            # Calculate raw False Positives (FP candidates)
            temp_df['Double'] = (temp_df['Rosette'] - temp_df['Intersect']).abs()
            temp_df['Unique'] = (temp_df['Union'] - temp_df['Rosette']).abs()

            # Special adjustment if metric is False Negatives
            if metric == 'False Negatives':
                # Example adjustment (replace 1000 with your actual total if needed)
                temp_df['Double'] = 1000 - temp_df['Double']
                temp_df['Unique'] = 1000 - temp_df['Unique']

            temp_df = temp_df.fillna(0)

            # Store metric summary
            metrics_dict[metric] = temp_df[['Combination', 'Double', 'Unique']]

            # Export intermediate results to Excel (raw numeric data)
            temp_df[['Combination', 'Union', 'Rosette', 'Intersect', 'Double', 'Unique']].to_excel(
                writer,
                sheet_name=metric,
                index=False
            )

        # After writing raw metrics, read them back for further calculations
    df_common_circs = pd.read_excel(output_path, sheet_name='Common Circles')
    df_FP = pd.read_excel(output_path, sheet_name='False Positives')
    df_FN = pd.read_excel(output_path, sheet_name='False Negatives')

    with pd.ExcelWriter(output_path, mode='a', if_sheet_exists='replace') as writer:
        # Calculate Precision, Recall, F-score
        for metric in ['Precision', 'Recall', 'F-score']:
            temp_df = pd.DataFrame()
            temp_df['Combination'] = df_common_circs['Combination'].values

            for combination in ['Union', 'Rosette', 'Intersect', 'Double', 'Unique']:
                tp = df_common_circs[combination].astype(float).values
                fp = df_FP[combination].astype(float).values
                fn = df_FN[combination].astype(float).values

                if metric == 'Precision':
                    with np.errstate(divide='ignore', invalid='ignore'):
                        values = tp / (tp + fp)
                    temp_df[combination] = values

                elif metric == 'Recall':
                    with np.errstate(divide='ignore', invalid='ignore'):
                        values = tp / (tp + fn)
                    temp_df[combination] = values

                elif metric == 'F-score':
                    with np.errstate(divide='ignore', invalid='ignore'):
                        precision = tp / (tp + fp)
                        recall = tp / (tp + fn)
                        f1 = 2 * (precision * recall) / (precision + recall)
                    temp_df[combination] = f1

            temp_df = temp_df.fillna(0)

            # Create copies for display (to replace some values with '-')
            temp_df_disp = temp_df.copy()
            df_common_circs_disp = df_common_circs.copy()
            df_FP_disp = df_FP.copy()
            df_FN_disp = df_FN.copy()

            # Apply "-" replacement only on the copies for display purposes
            for df_ in [df_common_circs_disp, df_FP_disp, df_FN_disp]:
                df_.loc[:14, 'Rosette'] = "-"
                df_.loc[:4, 'Intersect'] = "-"
                df_.loc[:14, 'Double'] = "-"
                df_.loc[:4, 'Unique'] = "-"

            temp_df_disp.loc[:14, 'Rosette'] = "-"
            temp_df_disp.loc[:4, 'Intersect'] = "-"
            temp_df_disp.loc[:14, 'Double'] = "-"
            temp_df_disp.loc[:4, 'Unique'] = "-"

            # Write back all to Excel (overwrite sheets with display versions)
            df_common_circs_disp.to_excel(writer, sheet_name='Common Circles', index=False)
            df_FP_disp.to_excel(writer, sheet_name='False Positives', index=False)
            df_FN_disp.to_excel(writer, sheet_name='False Negatives', index=False)
            temp_df_disp[['Combination', 'Union', 'Rosette', 'Intersect', 'Double', 'Unique']].to_excel(
                writer,
                sheet_name=metric,
                index=False
            )


def process_filtering(circular_bed, filter_name, combination, tools, tools_abbreviations, tool_abbreviations_r, real_combinations, ros_combinations_mapping, filtering_base_dir, output_base_dir):
    filtering_base_dir = filtering_base_dir + filter_name
    output_base_dir = output_base_dir + filter_name
    # Analyze circular combinations
    analyze_circular_combination(filtering_base_dir, circular_bed, output_base_dir, tools, tools_abbreviations)

    # Perform rosette operation
    circular_dir_ros = f'{filtering_base_dir}/intersect'
    output_dir_ros = f'{filtering_base_dir}/rosette'
    rosette(circular_bed, circular_dir_ros, output_dir_ros, tool_abbreviations_r)

    # Process combinations
    for combine in tools:
        combine_dir = f'{filtering_base_dir}/{combine}'
        calculate_statistics(combine_dir, circular_bed, output_base_dir, combination)

        # Paths for combination statistics
        union_path = f'{output_base_dir}/union_statistics.csv'
        rosette_path = f'{output_base_dir}/rosette_statistics.csv'
        intersect_path = f'{output_base_dir}/intersect_statistics.csv'
        output_path = f'{output_base_dir}/combination_stats.xlsx'

        # Generate combination stats
        generate_combination_stats(union_path, rosette_path, intersect_path, output_path, real_combinations, ros_combinations_mapping)


def plot_combinations(base_path):
    """
    Reads F-score, Precision, and Recall data from combination_stats.xlsx for different filtering methods,
    processes the data, and creates strip plots comparing the metrics.
    
    Parameters:
    base_path (str): Base directory where filtering method subfolders are located.
    """
    # Define filtering options and file paths
    filtering_options = ['unfilter', 'filter-split', 'filter-duplicates', 'filter']
    file_paths = {f: f"{base_path}/{f}/combination_stats.xlsx" for f in filtering_options}
    
    # Define columns of interest
    columns_of_interest = ["Combination", "Union", "Rosette", "Intersect", "Double", "Unique"]
    metrics = ["F-score", "Precision", "Recall"]  # Metrics to plot
    
    # Define color palette
    palette = {'unfilter': '#d46014', 
               'filter-split': '#ddcd3d', 
               'filter-duplicates': '#064b76ff', 
               'filter': '#63bdf6ff'}

    # Iterate over each metric
    for metric in metrics:
        data_list = []
        
        # Read and process data
        for filtering, file_path in file_paths.items():
            df = pd.read_excel(file_path, sheet_name=metric)  # Ensure correct sheet name
            df = df[columns_of_interest]
            df = df[df["Combination"].str.split("_").str.len() > 2]  # Filter by combination length
            df.columns = df.columns.astype(str)  # Convert column names to strings
            df_melted = df.melt(id_vars=["Combination"], var_name="Metric", value_name=metric)
            df_melted[metric] = pd.to_numeric(df_melted[metric], errors="coerce")  # Convert to numeric
            df_melted["Filtering"] = filtering  # Add filtering method column
            data_list.append(df_melted)
        
        # Combine all data
        final_df = pd.concat(data_list, ignore_index=True).dropna(subset=[metric])
        final_df["Metric"] = pd.Categorical(final_df["Metric"], 
                                            categories=["Union", "Rosette", "Intersect", "Double", "Unique"], 
                                            ordered=True)

        # Create figure
        plt.figure(figsize=(9, 4))

        # Stripplot with jitter
        sns.stripplot(x="Metric", y=metric, hue="Filtering", data=final_df, 
                      dodge=True, size=5, alpha=0.7, palette=palette)

        # Add median and MAD lines
        for m_index, m in enumerate(columns_of_interest[1:]):
            for f_index, f in enumerate(filtering_options):
                x_pos = m_index + (f_index - 2) * 0.2 + 0.1
                df_sub = final_df[(final_df['Metric'] == m) & (final_df['Filtering'] == f)]
                if not df_sub.empty:
                    median = np.median(df_sub[metric].values)
                    MAD = median_abs_deviation(df_sub[metric].values)

                    w = 0.07
                    plt.plot([x_pos - w, x_pos + w], [median, median], color="#232323", zorder=10)
                    plt.plot([x_pos + w, x_pos + w], [median + MAD, median - MAD], color="#bcbcbc", zorder=10)

        # Adjust labels and legend
        plt.xlabel('')
        plt.ylabel(metric, fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylim(0, 1.05)
        plt.legend(labels=['unfilter', 'filter-split', 'filter-duplicates', 'filter'], loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
        #plt.legend([], frameon=False)  
        plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
        sns.despine()
        plt.tight_layout()

        # Save and show plot
        save_path = os.path.join(base_path, f'combinations/combination_{metric.lower()}_stripplot.png')

        # Ensure the directory exists
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        # Save the plot
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

def compare_combinations(base_path):
    """
    Performs statistical comparisons of F-score, Precision, and Recall between different combination methods
    across four filtering strategies using Kruskal-Wallis and Dunn’s test.

    Parameters:
    base_path (str): Base directory where filtering method subfolders are located.

    Returns:
    None (prints test results)
    """

    # Define filtering methods and corresponding file paths
    filtering_options = ['unfilter', 'filter_split', 'filter_duplicates', 'filter']
    file_paths = {f: f"{base_path}/{f}/combination_stats.xlsx" for f in filtering_options}

    # Define metrics and combination methods
    metrics = ["F-score", "Precision", "Recall"]
    columns_of_interest = ["Combination", "Union", "Rosette", "Intersect", "Double", "Unique"]

    # Iterate over each filtering method
    for filtering, file_path in file_paths.items():
        print(f"\n--- Statistical Analysis for Filtering Method: {filtering} ---\n")
        
        # Iterate over each metric
        for metric in metrics:
            print(f"\n--- Statistical Analysis for {metric} ---")
            
            # Collect data for the current filtering method
            all_data = []

            try:
                # Read the data for the current filtering method
                df = pd.read_excel(file_path, sheet_name=metric)
                df = df[columns_of_interest]
                df = df[df["Combination"].str.split("_").str.len() > 2]  # Keep only valid combinations
                
                # Melt dataframe for statistical comparison
                df_melted = df.melt(id_vars=["Combination"], var_name="Method", value_name=metric)
                df_melted["Filtering"] = filtering  # Track filtering method
                all_data.append(df_melted)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                continue

            # Combine data for the current filtering method
            final_df = pd.concat(all_data, ignore_index=True).dropna(subset=[metric])

            # Perform Kruskal-Wallis test to check for overall differences within this filtering method
            grouped_data = [final_df[final_df["Method"] == m][metric].values for m in columns_of_interest[1:]]
            kruskal_stat, p_value = kruskal(*grouped_data)

            print(f"Kruskal-Wallis test for {metric}: H = {kruskal_stat:.4f}, p = {p_value:.4f}")

            # If significant, proceed with Dunn’s test
            if p_value < 0.05:
                print("Significant differences detected! Performing Dunn’s post-hoc test...\n")
                posthoc_results = sp.posthoc_dunn(final_df, val_col=metric, group_col="Method", p_adjust="bonferroni")
                print(posthoc_results)
            else:
                print("No significant differences found.")
            

# Function to merge files (Circle-Seq or any other type specified)
def merge_files(directories, data, file_type="Circle-Seq"):
    # Create a list to store dataframes
    file_dfs = []

    for dir_path in directories:
        # Find files matching the specified type
        files = [f for f in os.listdir(dir_path) if f.startswith(file_type)]
        
        # Read and append the files to the list
        for file in files:
            file_path = os.path.join(dir_path, file)
            df = pd.read_csv(file_path, sep='\t', header=None)  # Assuming tab-separated .bed files
            # Select first three columns
            df = df.iloc[:, :3]
            # Remove 'chr' from the first column and create the new first column
            df[0] = df[0].str.replace('chr', '', regex=False)
            
            # Update columns with the correct format for .bed: Chromosome, Start, End
            df.columns = ['chrom', 'start', 'end']
            
            # Create the 'circle' column by concatenating the chromosome, start, and end
            df['circle'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
            
            # Add the 'circle' column to the list of dataframes
            file_dfs.append(df['circle'])  # Add the 'circle' column to the list

    # Concatenate all the 'circle' columns into one DataFrame and drop duplicates
    if file_dfs:
        merged_df = pd.concat(file_dfs, axis=0, ignore_index=True).drop_duplicates()
        # Rename the column only for saving (optional)
        merged_df = merged_df.to_frame(name="Circles")
        # Define the output path for the merged file in the "unfilter" directory
        output_path = f"results/{data}/real/{file_type}/{file_type}.bed"  # Save it inside unfilter
        # Save the merged DataFrame to a .bed file
        merged_df.to_csv(output_path, sep='\t', index=False, header=True)
        print(f"Merged {file_type} file saved in 'unfilter' folder: {output_path}")
    else:
        print(f"No valid {file_type} files to merge.")

def save_bed_files(circles_df, output_dir, filename):
    # Assuming circles_df is a DataFrame, extract necessary columns for BED format
    os.makedirs(output_dir, exist_ok=True)
    circles_df.to_csv(output_dir + '/' + filename, index=False, sep=',')

def analyze_and_save_circles(main_matrix_path, extra_matrix_path, ros_combinations_mapping, eccDNA_tools_abbreviations, output_base):
    # Load both matrices
    df_main = pd.read_csv(main_matrix_path)
    df_main.rename(columns={"Circles": "Circle"}, inplace=True)
    df_extra = pd.read_csv(extra_matrix_path)

    # Convert tool names to abbreviations
    tool_abbreviations = {v: k for k, v in eccDNA_tools_abbreviations.items()}

    # Define output directories
    output_dirs = { 
        "unique": os.path.join(output_base, "unique"),
        "intersect": os.path.join(output_base, "intersect"),
        "union": os.path.join(output_base, "union"),
        "rosette": os.path.join(output_base, "rosette"),
        "double": os.path.join(output_base, "double"),
    }

    # Iterate through comparisons
    for comp_name, comp_str in ros_combinations_mapping.items():
        tools = comp_name.split('_')  # Example: ['CE', 'CM', 'CF']
        
        # Map abbreviations to tool names
        tool_columns = [tool_abbreviations[tool] for tool in tools if tool in tool_abbreviations]
        
        # Compute metrics
        df_main['sum'] = df_main[tool_columns].sum(axis=1)

        # Categorize circles based on sum of tool columns
        unique_circles = df_main[df_main['sum'] == 1]["Circle"].tolist()
        intersect_circles = df_main[df_main['sum'] == len(tool_columns)]["Circle"].tolist()
        union_circles = df_main[df_main['sum'] >= 1]["Circle"].tolist()
        rosette_circles = list(set(union_circles) - set(unique_circles))
        double_circles = list(set(rosette_circles) - set(intersect_circles))

        # Merge these categorized circles with the extra matrix
        unique_circles_merged = df_extra[df_extra['Circle'].isin(unique_circles)]
        intersect_circles_merged = df_extra[df_extra['Circle'].isin(intersect_circles)]
        union_circles_merged = df_extra[df_extra['Circle'].isin(union_circles)]
        rosette_circles_merged = df_extra[df_extra['Circle'].isin(rosette_circles)]
        double_circles_merged = df_extra[df_extra['Circle'].isin(double_circles)]

        # Save the merged dataframes into separate files (you can adjust the file format)
        save_bed_files(unique_circles_merged, output_dirs["unique"], f"{comp_name}.bed")
        save_bed_files(intersect_circles_merged, output_dirs["intersect"], f"{comp_name}.bed")
        save_bed_files(union_circles_merged, output_dirs["union"], f"{comp_name}.bed")
        save_bed_files(rosette_circles_merged, output_dirs["rosette"], f"{comp_name}.bed")
        save_bed_files(double_circles_merged, output_dirs["double"], f"{comp_name}.bed")