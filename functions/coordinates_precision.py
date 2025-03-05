import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from scipy import stats
from functions.load_bed_file import load_bed_file
from functions.compare_circular_data import are_circles_equal
import csv

def compare_bed_files_threshold(file1_list, file2_list, thresholds):
    """
    Compare multiple bed files with another bed file and find common circles for different thresholds.

    Parameters:
        file1_list (list): List of paths to bed files for comparison.
        file2_list (list): List of paths to bed files for comparison (detection bed files). Should contain only one path.
        thresholds (list): List of thresholds for comparison.

    Returns:
        dict: Dictionary containing lists of common circle counts for each threshold for each file1.
    """
    if len(file1_list) != 1:
        raise ValueError("file1_list should contain exactly one file path.")

    common_counts = {file2: {threshold: [] for threshold in thresholds} for file2 in file2_list}
    circle_data1 = load_bed_file(file1_list[0])

    for file2 in file2_list:
        circle_data2 = load_bed_file(file2)
        for threshold in thresholds:
            common_circles = []
            for circle1 in circle_data1:
                for circle2 in circle_data2:
                    if are_circles_equal(circle1, circle2, threshold):
                        common_circles.append(circle1)
                        break
            common_counts[file2][threshold].append(len(common_circles))

    return common_counts

def plot_coordinates_precision(tools, coverages, circular_bed, circular_dir, output_dir, circle_type):
    """
    Function to plot common circle counts for different thresholds across multiple tools.

    Parameters:
        tools (list): List of tool names.
        coverages (list): List of coverage values.
        circular_bed (str): Path to the circular bed file.
        circular_dir (str): Directory containing circular detection results for each tool.
        output_dir (str): Directory where results and statistics will be saved.
        circle_type (str): Circle type.
    """
    # Ensure the save directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Define thresholds (exclude 0 to avoid division by zero)
    thresholds = [1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 50]
    
    # Define the custom color palette
    custom_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']
    
    sqrt_a_values = {tool: [] for tool in tools}  # Store sqrt(a) values for each tool

    # Create a figure with subplots arranged horizontally
    num_tools = len(tools)
    fig, axes = plt.subplots(1, num_tools, figsize=(8 * num_tools, 6), sharey=True)
    
    # If there's only one subplot, axes is not an array, so we need to make it an array
    if num_tools == 1:
        axes = [axes]

    # Prepare a list to store the printed results for CSV
    print_results = []

    # Loop over each tool and plot the common circle counts in the corresponding subplot
    for i, tool in enumerate(tools):
        # Generate paths for the current folder
        folder_path = f'{circular_dir}/{tool}/'
        circular_detection_list = [f'{folder_path}cov{cov}_{tool}.bed' for cov in coverages]

        # Compare bed files for different thresholds
        common_counts = compare_bed_files_threshold([circular_bed], circular_detection_list, thresholds)

        # Plot common counts for different thresholds for each file
        for j, (file1, counts_dict) in enumerate(common_counts.items()):
            counts = list(counts_dict.values())  # Convert dictionary values to a list

            # Flatten counts array in case it is 2D
            counts = np.ravel(counts)  # Ensure 'counts' is 1D

            if len(counts) != len(thresholds):
                print(f"Error: counts length ({len(counts)}) does not match thresholds length ({len(thresholds)})")
                continue

            # Prepare transformed x (1/x)
            inv_thresholds = 1 / np.array(thresholds)

            try:
                # Perform linear regression on transformed data
                slope, intercept, r_value, p_value, std_err = stats.linregress(inv_thresholds, counts)
                a = -slope  # 'a' in the original equation
                b = intercept  # 'b' in the original equation

                # Calculate the square root of 'a'
                a_sqrt = np.sqrt(a)

                # Append sqrt(a) value for the current tool
                sqrt_a_values[tool].append(a_sqrt)

                # Print the parameters 'a', 'sqrt(a)' along with the corresponding coverage
                result = {
                    'Tool': tool,
                    'Coverage': coverages[j],
                    'b': b,
                    'a': a,
                    'sqrt(a)': a_sqrt
                }
                print_results.append(result)  # Add result to the list

                # Extract label from the filename
                label = os.path.basename(file1).split('_')[0]

                # Plot the actual data in the corresponding subplot
                axes[i].plot(thresholds, counts, marker='o', color=custom_palette[j % len(custom_palette)], label=label)

            except ValueError as e:
                print(f"Error fitting curve for {tool}: {e}")
                continue

        # Customize each subplot
        axes[i].set_xlabel('Threshold (bp)', fontsize=12)
        axes[i].set_ylabel(f'# True {circle_type}', fontsize=12)
        axes[i].set_title(f'{tool}', fontsize=14)
        axes[i].grid(True, linestyle='--', alpha=0.7)
        axes[i].set_ylim(0, 1000)  # Limitar el eje Y a 1000
        sns.despine(top=True, right=True)
    
    # Add a single legend with numeric coverage values
    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', markersize=8) 
            for j, color in enumerate(custom_palette)]
    labels = coverages
    fig.legend(handles, labels, title="Coverage", frameon=False, fontsize=10, 
            loc='upper left', bbox_to_anchor=(0.95, 0.65), ncol=1)

    # Set a common title for the entire figure
    fig.suptitle(f'Coordinates precision of {circle_type} detection software', fontsize=16)

    # Adjust layout and save the plot
    plt.tight_layout(rect=[0, 0, 0.95, 1])  # Adjust layout to make room for the suptitle and legend
    save_path = os.path.join(output_dir, 'common_circles_count_all_tools.png')
    plt.savefig(save_path, dpi=300)

    # Show the combined plot
    plt.show()

    # Save the results in a CSV file with rounded values
    csv_file_path = os.path.join(output_dir, 'precision_analysis_results.csv')
    header = ['Tool', 'Coverage', 'sqrt(a)', 'b']

    # Open the CSV file for writing
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)

        # Write each result from print_results, rounding the values to 2 decimal places
        for result in print_results:
            rounded_result = [result['Tool'], result['Coverage'], round(result['sqrt(a)'], 2), round(result['b'], 2)]
            writer.writerow(rounded_result)

    # Clear the figure for the next iteration
    plt.clf()  # Clear the figure to avoid overlap

    ### NEW GRAPH FOR sqrt(a) vs b ###

    plt.figure(figsize=(12,4))

    for i, tool in enumerate(tools):
        folder_path = f'{circular_dir}/{tool}/'
        circular_detection_list = [f'{folder_path}cov{cov}_{tool}.bed' for cov in coverages]

        common_counts = compare_bed_files_threshold([circular_bed], circular_detection_list, thresholds)

        sqrt_a_values = []
        b_values = []

        for j, (file1, counts_dict) in enumerate(common_counts.items()):
            counts = np.ravel(list(counts_dict.values()))

            if len(counts) != len(thresholds):
                print(f"Error: counts length ({len(counts)}) does not match thresholds length ({len(thresholds)})")
                continue

            inv_thresholds = 1 / np.array(thresholds)

            try:
                slope, intercept, _, _, _ = stats.linregress(inv_thresholds, counts)
                a = -slope
                b = intercept

                sqrt_a_values.append(np.sqrt(a))
                b_values.append(b)
            except ValueError as e:
                print(f"Error fitting curve for {tool}: {e}")
                continue

        sorted_indices = np.argsort(coverages)
        sqrt_a_values = np.array(sqrt_a_values)[sorted_indices]
        b_values = np.array(b_values)[sorted_indices]
        coverages_sorted = np.array(coverages)[sorted_indices]

        # Get the color for the current tool
        tool_color = custom_palette[i]

        # Apply a color gradient for the line segments based on coverage
        for j in range(1, len(sqrt_a_values)):
            # Define line segment's color intensity based on coverage (gradually increasing)
            coverage_ratio = (coverages_sorted[j] - min(coverages_sorted)) / (max(coverages_sorted) - min(coverages_sorted))

            # Set the alpha value based on coverage: start more transparent and become more opaque
            alpha_value = 0.3 + 0.7 * coverage_ratio  # Alpha goes from 0.3 (more transparent) to 1 (fully opaque)

            # Plot the line segment with the adjusted transparency (alpha)
            # No marker is specified, so it should only plot the line
            line = plt.plot(sqrt_a_values[j-1:j+1], b_values[j-1:j+1], color=tool_color, linewidth=3, alpha=alpha_value)

            # Add the label for the legend once per tool
            if j == 1:  # Only add the legend label for the first line segment
                line[0].set_label(tool)

    
    plt.ylim(0, 1050)  
    
    # Create a custom legend with circles (using Line2D for the legend)
    legend_elements = [mlines.Line2D([0], [0], color=tool_color, label=tool) for tool, tool_color in zip(tools, custom_palette)]

    # Add the custom legend to the plot
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
    
    # Set larger tick labels
    plt.xticks(fontsize=16)  # Increase x-tick size
    plt.yticks(fontsize=16)  # Increase y-tick size

    # Display the plot
    plt.xlabel(r'$\sqrt{a}$', fontsize=16)
    plt.ylabel('b', fontsize=16)
    #plt.title(f'Precision Analysis of {circle_type} Detection Tools', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)
    sns.despine()

    # Save the plot
    save_path = os.path.join(output_dir, 'precision_analysis.png')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(save_path, dpi=300)
    plt.show()