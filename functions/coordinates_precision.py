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
    This function evaluates the precision of coordinate detection for circular DNA or RNA across various tools and coverage levels.
    It uses a known simulated BED file (circular_bed) as the reference and compares it to detected circles.
    
    For each tool and coverage level:
    - The function computes overlap with the ground truth at various thresholds.
    - Fits a linear regression to assess the asymptotic maximum number of detected circles (intercept at 1/threshold → 0).
    - Computes precision-like metrics:
        - count1_over_b: count at T=1 over the estimated maximum (b)
        - count2_over_b: count at T=2 over the estimated maximum (b)
        - count1_over_count2: ratio between T=1 and T=2
    
    These are plotted across tools and coverages to visualize precision trends.

    Parameters:
        tools (list): Names of tools used to detect circular DNA/RNA.
        coverages (list): Coverage levels to evaluate (e.g., [5, 10, 30]).
        circular_bed (str): Path to the simulated BED file with true circles.
        circular_dir (str): Directory containing BED files for each tool/coverage.
        output_dir (str): Directory to save output plots.
        circle_type (str): Label for the circle type (used in axis labels).
    """
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define thresholds for comparing circle overlaps (e.g., for precision measurement)
    thresholds = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]

    # Define color palettes
    colorblind_palette = sns.color_palette('colorblind')  # not used below
    custom_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']

    # For each tool, precompute common counts between the ground truth and detected circles
    for tool in tools:
        folder_path = f'{circular_dir}/{tool}/'
        # Generate file paths for each coverage level
        circular_detection_list = [f'{folder_path}cov{cov}_{tool}.bed' for cov in coverages]
        # Compare each file against the ground truth using thresholds
        common_counts = compare_bed_files_threshold([circular_bed], circular_detection_list, thresholds)

    # Define the metrics to compute and plot
    for metric_name, x_values_func in {
        'count1_over_b': lambda count1, count2, b: count1 / b if b != 0 else np.nan,
        'count2_over_b': lambda count1, count2, b: count2 / b if b != 0 else np.nan,
        'count1_over_count2': lambda count1, count2, b: count1 / count2 if count2 != 0 else np.nan,
    }.items():

        plt.figure(figsize=(5, 4))  # Create a new figure for each metric

        for i, tool in enumerate(tools):
            folder_path = f'{circular_dir}/{tool}/'
            circular_detection_list = [f'{folder_path}cov{cov}_{tool}.bed' for cov in coverages]
            common_counts = compare_bed_files_threshold([circular_bed], circular_detection_list, thresholds)

            x_values = []  # Metric values (x-axis)
            b_values = []  # Max circle count estimates (y-axis)

            for j, (file1, counts_dict) in enumerate(common_counts.items()):
                counts_list = list(counts_dict.values())
                if len(counts_list) != len(thresholds):
                    continue  # Skip if the count list doesn't match the number of thresholds

                counts = np.ravel(counts_list)
                inv_thresholds = 1 / np.array(thresholds)  # Inverse thresholds for linear modeling

                try:
                    # Fit a line: counts = a*(1/T) + b => b is used as max estimated number of circles
                    slope, intercept, _, _, _ = stats.linregress(inv_thresholds, counts)
                    a = -slope
                    b = intercept
                    count1 = counts_dict[1][0]  # Count at threshold T=1
                    count2 = counts_dict[2][0]  # Count at threshold T=2

                    x_value = x_values_func(count1, count2, b)  # Compute selected metric
                    if not np.isnan(x_value):
                        x_values.append(x_value)
                        b_values.append(b)
                except Exception:
                    continue  # Skip failed regressions

            # Sort values according to coverage for smooth plotting
            sorted_indices = np.argsort(coverages)
            x_values = np.array(x_values)[sorted_indices]
            b_values = np.array(b_values)[sorted_indices]
            coverages_sorted = np.array(coverages)[sorted_indices]

            tool_color = custom_palette[i]

            # Plot line segments between each pair of points
            for j in range(1, len(x_values)):
                coverage_ratio = (coverages_sorted[j] - min(coverages_sorted)) / (max(coverages_sorted) - min(coverages_sorted))
                alpha_value = 0.3 + 0.7 * coverage_ratio  # Transparency based on coverage
                line = plt.plot(x_values[j-1:j+1], b_values[j-1:j+1], color=tool_color, linewidth=3, alpha=alpha_value)
                if j == 1:
                    line[0].set_label(tool)  # Add label only once per tool

        # Plot formatting
        plt.ylim(0, 1050)
        plt.xlim(0, 1.05)
        plt.axhline(y=1000, linestyle='--', linewidth=1, color='grey')  # Reference line at 1000
        plt.legend([], frameon=False, fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        # Axis labels and output filename per metric
        if metric_name == 'count1_over_b':
            xlabel = rf'$\frac{{\#\, {circle_type}_{{\mathrm{{T}}=1}}}}{{\#\, {circle_type}_{{\mathrm{{max}}}}}}$'
            filename = 'count1_over_b.png'
        elif metric_name == 'count2_over_b':
            xlabel = rf'$\frac{{\#\, {circle_type}_{{\mathrm{{T}}=2}}}}{{\#\, {circle_type}_{{\mathrm{{max}}}}}}$'
            filename = 'count2_over_b.png'
        elif metric_name == 'count1_over_count2':
            xlabel = rf'$\frac{{\#\, {circle_type}_{{T=1}}}}{{\#\, {circle_type}_{{T=2}}}}$'
            filename = 'count1_over_count2.png'
        else:
            xlabel = 'Metric'
            filename = f'{metric_name}.png'

        # Finalize and save plot
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel(rf'$\#\, {circle_type}_{{\mathrm{{max}}}}$', fontsize=16)
        plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
        sns.despine()
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.show()
