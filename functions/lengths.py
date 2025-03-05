import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import seaborn as sns
from scipy.stats import ks_2samp
from matplotlib.patches import Rectangle

def plot_length_distributions_kde(circular_bed, filter_directory, tools, output_base_name, circle_type):
    """
    Generate KDE plots comparing the distributions of filtered detections for multiple tools.

    Parameters:
        circular_bed (str): Path to the known circular BED file (ground truth).
        filter_directory (str): Directory containing filtered detection BED files.
        tools (list): List of tool names.
        output_base_name (str): Base name for output files.
        circle_type (str): Circle type.
    
    Returns:
        None
    """
     # Load lengths from known circles BED file
    known_lengths = np.loadtxt(circular_bed, usecols=(2,), dtype=int) - np.loadtxt(circular_bed, usecols=(1,), dtype=int)
    known_lengths_df = pd.DataFrame({'Length': known_lengths, 'Tool': 'Simulated'})
    
    # Prepare a DataFrame for all tool lengths (starting with known circles as "Simulated")
    all_data = known_lengths_df.copy()

    # Loop through each tool and extract filtered lengths
    for tool in tools:
        filtered_file = os.path.join(filter_directory, tool, f'cov30_{tool}.bed')
        if os.path.exists(filtered_file):
            filtered_lengths = np.loadtxt(filtered_file, usecols=(2,), dtype=int) - np.loadtxt(filtered_file, usecols=(1,), dtype=int)
            tool_df = pd.DataFrame({'Length': filtered_lengths, 'Tool': tool})
            all_data = pd.concat([all_data, tool_df], ignore_index=True)

    # Import the colorblind palette from seaborn
    color_palette = ['#d46014', '#ddcd3d', '#0972b3', '#4fb4f5', '#b54582']

    # Manually set the first color to black, followed by the colorblind palette
    custom_palette = ['b3b3b3'] + color_palette

    # 1. Plot for lengths < 10000 (main plot)
    g = sns.displot(
        data=all_data[all_data['Length'] < 10000], 
        x="Length", 
        hue="Tool", 
        kind="kde", 
        aspect=1, 
        bw_adjust=0.1,
        palette=custom_palette
    )

    # Adjust line properties for the simulated circles (make it black and dashed)
    for line in g.axes[0][0].lines:
        if line.get_label() == 'Simulated':
            line.set_color('black')
            line.set_linestyle('--')

    plt.xlim([200, 10000])  # Set x-axis limits
    plt.xlabel('Length')
    plt.ylabel('Density')
    plt.title(f'{circle_type} - Length Distribution Comparison')  # Title includes circle_type

    # Save the full-length plot with circle_type in the filename
    output_directory = output_base_name
    os.makedirs(output_directory, exist_ok=True)
    output_path = os.path.join(output_directory, f'{circle_type}_length_distributions_comparison.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    print(f"Full length distribution plot saved to {output_path}")

    # 2. Plot for lengths < 1000 (short length plot)
    g_short = sns.displot(
        data=all_data[all_data['Length'] < 1000], 
        x="Length", 
        hue="Tool", 
        kind="kde", 
        aspect=1, 
        bw_adjust=0.1,
        palette=custom_palette
    )

    # Adjust line properties for the simulated circles (make it black and dashed)
    for line in g_short.axes[0][0].lines:
        if line.get_label() == 'Simulated':
            line.set_color('black')
            line.set_linestyle('--')

    plt.xlabel('Length')
    plt.ylabel('Density')
    plt.title(f'{circle_type} - Length Distribution Comparison of Short Circles')  # Title includes circle_type
    plt.xlim(200, 1000)  # Set x-axis limits for short circles

    # Save the short-length plot with circle_type in the filename
    output_path_short = os.path.join(output_directory, f'{circle_type}_short_length_distributions_comparison.png')
    plt.savefig(output_path_short, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    print(f"Short length distribution plot saved to {output_path_short}")

    # 3. Perform Kolmogorov-Smirnov test for full length data (<10000)
    simulated_lengths = all_data[(all_data['Tool'] == 'Simulated') & (all_data['Length'] < 10000)]['Length']
    ks_results_full = []

    # Compare simulated lengths with each tool's filtered lengths using KS test
    for tool in tools:
        tool_lengths = all_data[(all_data['Tool'] == tool) & (all_data['Length'] < 10000)]['Length']
        if len(tool_lengths) > 0:  # Ensure there's data to compare
            ks_stat, p_value = ks_2samp(simulated_lengths, tool_lengths)
            ks_results_full.append((tool, ks_stat, p_value))
            print(f"KS test (full) between Simulated and {tool}: KS Statistic = {ks_stat}, p-value = {p_value}")

    # Save full-length KS test results to a CSV file with circle_type in the filename
    ks_df_full = pd.DataFrame(ks_results_full, columns=['Tool', 'KS Statistic', 'p-value'])
    ks_df_full['Circle Type'] = circle_type  # Add circle_type to the DataFrame
    ks_output_path_full = os.path.join(output_directory, f'{circle_type}_ks_test_results_full.csv')
    ks_df_full.to_csv(ks_output_path_full, index=False)

    print(f"Kolmogorov-Smirnov test results (full) saved to {ks_output_path_full}")

    # 4. Perform Kolmogorov-Smirnov test for short lengths (<1000)
    simulated_lengths_short = all_data[(all_data['Tool'] == 'Simulated') & (all_data['Length'] < 1000)]['Length']
    ks_results_short = []

    # Compare simulated short lengths with each tool's filtered short lengths using KS test
    for tool in tools:
        tool_lengths_short = all_data[(all_data['Tool'] == tool) & (all_data['Length'] < 1000)]['Length']
        if len(tool_lengths_short) > 0:  # Ensure there's data to compare
            ks_stat, p_value = ks_2samp(simulated_lengths_short, tool_lengths_short)
            ks_results_short.append((tool, ks_stat, p_value))
            print(f"KS test (short) between Simulated and {tool}: KS Statistic = {ks_stat}, p-value = {p_value}")

    # Save short-length KS test results to a CSV file with circle_type in the filename
    ks_df_short = pd.DataFrame(ks_results_short, columns=['Tool', 'KS Statistic', 'p-value'])
    ks_output_path_short = os.path.join(output_directory, f'{circle_type}_ks_test_results_short.csv')
    ks_df_short.to_csv(ks_output_path_short, index=False)

    print(f"Kolmogorov-Smirnov test results (short) saved to {ks_output_path_short}")

    return output_path, output_path_short, ks_output_path_full, ks_output_path_short



def plot_length_distributions(circular_bed, filter_directory, tools, output_base_name, circle_type, draw_rectangle=True, min_length=320, max_length=480):
    """
    Generate KDE plots comparing the relative and absolute length distributions of filtered detections for multiple tools.

    Parameters:
        circular_bed (str): Path to the known circular BED file (ground truth).
        filter_directory (str): Directory containing filtered detection BED files.
        tools (list): List of tool names.
        output_base_name (str): Base name for output files.
        circle_type (str): Circle type.
        draw_rectangle (bool): Whether to draw the colored rectangle for 320-480 bp range (default is True).
        min_length (int): Minimum length of the rectangle (default is 320).
        max_length (int): Maximum length of the rectangle (default is 480).

    Returns:
        None
    """
    # Import the colorblind palette from seaborn
    color_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']

    # Load lengths from known circles BED file (Simulated data)
    known_lengths = np.loadtxt(circular_bed, usecols=(2,), dtype=int) - np.loadtxt(circular_bed, usecols=(1,), dtype=int)
    known_lengths_df = pd.DataFrame({'Length': known_lengths, 'Tool': 'Simulated'})

    # Prepare a DataFrame for all tool lengths
    all_data = known_lengths_df.copy()

    for tool in tools:
        filtered_file = os.path.join(filter_directory, tool, f'cov30_{tool}.bed')
        if os.path.exists(filtered_file):
            try:
                filtered_lengths = np.loadtxt(filtered_file, usecols=(2,), dtype=int) - np.loadtxt(filtered_file, usecols=(1,), dtype=int)
                tool_df = pd.DataFrame({'Length': filtered_lengths, 'Tool': tool})
                all_data = pd.concat([all_data, tool_df], ignore_index=True)
            except Exception as e:
                print(f"Error reading {filtered_file}: {e}")
        else:
            print(f"Warning: File {filtered_file} not found. Skipping {tool}.")

    def plot_relative_distribution(data, length_range, bins, output_filename, title, draw_rectangle, min_length, max_length, is_short=True):
        """
        Plot relative length distributions with optional rectangle highlighting.

        Parameters:
            data (pd.DataFrame): DataFrame containing 'Length' and 'Tool' columns.
            length_range (tuple): Range of lengths to plot (min_length, max_length).
            bins (array): Bin edges for the length intervals.
            output_filename (str): Path to save the plot.
            title (str): Title of the plot.
            draw_rectangle (bool): Whether to draw a rectangle highlighting the length range.
            min_length (int): Minimum length for the rectangle.
            max_length (int): Maximum length for the rectangle.
            is_short (bool): Whether the data is for short lengths (default is True).

        Returns:
            None
        """
        # Filter data within the specified range
        data = data[(data['Length'] >= length_range[0]) & (data['Length'] < length_range[1])].copy()

        # Create intervals
        data['Length_Interval'] = pd.cut(data['Length'], bins, right=False)
        data['Length_Interval'] = data['Length_Interval'].apply(lambda x: int(x.left))

        # Count the number of circles per interval for each tool
        length_counts = data.groupby(['Tool', 'Length_Interval']).size().reset_index(name='Count')

        # Get the counts for Simulated lengths in the same intervals
        simulated_counts = length_counts[length_counts['Tool'] == 'Simulated'][['Length_Interval', 'Count']]
        simulated_counts.rename(columns={'Count': 'Simulated_Count'}, inplace=True)

        # Merge with the simulated counts to calculate the relative counts
        relative_data = length_counts.merge(simulated_counts, on='Length_Interval', how='left')
        relative_data['Relative'] = relative_data['Count'] / relative_data['Simulated_Count']

        # Remove simulated data from the plot
        relative_data = relative_data[relative_data['Tool'] != 'Simulated']

        # Plot relative length distributions
        plt.figure(figsize=(14, 5))
        g = sns.lineplot(data=relative_data, x="Length_Interval", y="Relative", hue="Tool", marker='o', palette=color_palette)

        # Customize the x-axis labels
        if is_short:
            # For short data, set ticks to [200, 400, 600, 800, 1000]
            g.set_xticks([200, 400, 600, 800, 1000])
            g.set_xticklabels(['200', '400', '600', '800', '1000'])
            plt.xlim(100, 1050) 
        else:
            # For long data, set ticks to [2000, 4000, 6000, 8000, 10000] and limit x-axis to [1000, 10000]
            g.set_xticks([2000, 4000, 6000, 8000, 10000])
            g.set_xticklabels(['2000', '4000', '6000', '8000', '10000'])
            plt.xlim(500, 10500)  # Cut the graph at 1000

        plt.ylim([0, 2])
        plt.axhline(1, color='#b3b3b3', linestyle='--', label='Ideal Detection')  # Original grey line for "Ideal Detection"
        plt.xlabel('Length (bp)', fontsize=16)
        plt.ylabel('# Relative Count', fontsize=16)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.6)

        # Draw the rectangle if enabled
        if draw_rectangle:
            plt.gca().add_patch(Rectangle((min_length, 0), max_length - min_length, 2, linewidth=0, edgecolor='grey', facecolor='lightgrey', alpha=0.3))  # Lighter grey

        # Save the plot
        sns.despine()
        plt.xticks(fontsize=16)  # Increase x-tick size
        plt.yticks(fontsize=16)  # Increase y-tick size
        os.makedirs(output_base_name, exist_ok=True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit the legend
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()

        print(f"Relative length distribution plot saved to {output_filename}")


    def plot_absolute_distribution(data, length_range, bins, output_filename, title, smoothing_window, min_length, max_length, draw_rectangle, is_short=True):
        """
        Plot absolute length distributions with optional smoothing and rectangle highlighting.

        Parameters:
            data (pd.DataFrame): DataFrame containing 'Length' and 'Tool' columns.
            length_range (tuple): Range of lengths to plot (min_length, max_length).
            bins (array): Bin edges for the length intervals.
            output_filename (str): Path to save the plot.
            title (str): Title of the plot.
            smoothing_window (int): Window size for smoothing.
            min_length (int): Minimum length for the rectangle.
            max_length (int): Maximum length for the rectangle.
            draw_rectangle (bool): Whether to draw a rectangle highlighting the length range.
            is_short (bool): Whether the data is for short lengths (default is True).

        Returns:
            None
        """
        # Filter data within the specified range
        data = data[(data['Length'] >= length_range[0]) & (data['Length'] < length_range[1])].copy()

        # Create intervals
        data['Length_Interval'] = pd.cut(data['Length'], bins, right=False)
        data['Length_Interval'] = data['Length_Interval'].apply(lambda x: int(x.left))

        # Count the number of circles per interval for each tool
        length_counts = data.groupby(['Tool', 'Length_Interval']).size().reset_index(name='Count')

        # Apply smoothing to each tool's counts
        length_counts['Count'] = length_counts.groupby('Tool')['Count'].transform(lambda x: x.rolling(window=smoothing_window, min_periods=1, center=True).mean())

        # Plot absolute length distributions
        plt.figure(figsize=(14, 5))

        # Map colors to tools
        color_mapping = {"Simulated": "#b3b3b3"}
        tool_list = length_counts['Tool'].unique()
        color_idx = 0  # Index for other tools

        for tool in tool_list:
            if tool != "Simulated":
                color_mapping[tool] = color_palette[color_idx]
                color_idx += 1  # Move to the next color in the palette

        g = sns.lineplot(data=length_counts, x="Length_Interval", y="Count", hue="Tool", marker='o', palette=color_mapping)

        # Customize the x-axis labels
        if is_short:
            # For short data, set ticks to [200, 400, 600, 800, 1000]
            g.set_xticks([200, 400, 600, 800, 1000])
            g.set_xticklabels(['200', '400', '600', '800', '1000'])
            plt.xlim(100, 1050) 
        else:
            # For long data, set ticks to [2000, 4000, 6000, 8000, 10000] and limit x-axis to [1000, 10000]
            g.set_xticks([0, 2000, 4000, 6000, 8000, 10000])
            g.set_xticklabels(['0', '2000', '4000', '6000', '8000', '10000'])
            plt.xlim(100, 10500)  # Cut the graph at 1000

        plt.xlabel('Length (bp)', fontsize=16)
        plt.ylabel('# Count', fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.6)

        # Draw the rectangle if enabled
        if draw_rectangle:
            plt.gca().add_patch(Rectangle((min_length, 0), max_length - min_length, max(length_counts['Count']), linewidth=0, edgecolor='grey', facecolor='lightgrey', alpha=0.3))

        # Ensure 'Simulated' is the first label in the legend
        handles, labels = plt.gca().get_legend_handles_labels()
        sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: x[1] != "Simulated")
        sorted_handles, sorted_labels = zip(*sorted_handles_labels)
        plt.legend(sorted_handles, sorted_labels, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)

        # Save the plot
        sns.despine()
        plt.xticks(fontsize=16)  # Increase x-tick size
        plt.yticks(fontsize=16)  # Increase y-tick size
        os.makedirs(output_base_name, exist_ok=True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit the legend
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()

        print(f"Absolute length distribution plot saved to {output_filename}")

    # Plot short lengths (< 1000 bp) with a colored rectangle between min_length-max_length (if draw_rectangle=True)
    plot_relative_distribution(
        all_data,
        length_range=(140, 1040),
        bins=np.arange(140, 1040, 20),
        output_filename=os.path.join(output_base_name, 'relative_length_distributions_short.png'),
        title=f'Relative Length Distribution for {circle_type} Detection Software (Short Lengths)',
        draw_rectangle=draw_rectangle,
        min_length=min_length,
        max_length=max_length,
        is_short=True  # Indicate this is for short lengths
    )

    # Plot long lengths (>= 1000 bp)
    plot_relative_distribution(
        all_data,
        length_range=(140, 11000),
        bins=np.arange(750, 11000, 500),
        output_filename=os.path.join(output_base_name, 'relative_length_distributions_long.png'),
        title=f'Relative Length Distribution for {circle_type} Detection Software (Long Lengths)',
        draw_rectangle=False,
        min_length=min_length,
        max_length=max_length,
        is_short=False  # Indicate this is for long lengths
    )

    # Plot short lengths (< 1000 bp) without the rectangle for absolute counts
    plot_absolute_distribution(
        all_data,
        length_range=(140, 1040),
        bins=np.arange(140, 1040, 20),
        output_filename=os.path.join(output_base_name, 'absolute_length_distributions_short.png'),
        title=f'Absolute Length Distribution for {circle_type} Detection Software (Short Lengths)',
        smoothing_window=5,  # Example smoothing window size
        min_length=min_length,  # Example min length
        max_length=max_length,  # Example max length
        draw_rectangle=draw_rectangle,  # Set to True if you want to draw the rectangle
        is_short=True  # Indicate this is for short lengths
    )

    # Plot long lengths (>= 1000 bp) without the rectangle for absolute counts
    plot_absolute_distribution(
        all_data,
        length_range=(140, 11000),
        bins=np.arange(500, 11000, 500),
        output_filename=os.path.join(output_base_name, 'absolute_length_distributions_long.png'),
        title=f'Absolute Length Distribution for {circle_type} Detection Software (Long Lengths)',
        smoothing_window=5,  # Example smoothing window size
        min_length=min_length,  # Example min length
        max_length=max_length,  # Example max length
        draw_rectangle=False,  # Set to True if you want to draw the rectangle
        is_short=False  # Indicate this is for long lengths
    )

# Function to apply smoothing using a sliding window
def apply_smoothing(data, window_size):
    """
    Apply smoothing to the data using a rolling window.

    Parameters:
        data (pd.Series): The data to be smoothed.
        window_size (int): The size of the rolling window.

    Returns:
        pd.Series: The smoothed data.
    """
    return data.rolling(window=window_size, min_periods=1, center=True).mean()