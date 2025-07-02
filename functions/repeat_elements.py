import os
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from functions.compare_circular_data import calculate_metrics

def annotate_repeat_elements(repeats_file, bed_files, output_dir, tools, true_or_false):
    """
    Annotates BED files with repeat element information and saves them.

    Parameters:
        repeats_file (str): Path to repeat elements data.
        bed_files (list): List of BED file paths to annotate.
        output_dir (str): Directory to save annotated files.
        tools (list): List of tool names for output filenames.
        true_or_false (str): Label for the data type ('truepositives', 'falsenegatives', 'falsepositives').

    This function annotates BED files with repeat types and saves them to the specified directory.
    """
    # Read the repeats file
    repeats_df = pd.read_csv(repeats_file, sep='\t', header=None, names=['chr', 'start', 'end', 'type'])
    
    # Adjust chromosome format in repeats_df to match bed_df
    repeats_df['chr'] = repeats_df['chr'].str.replace('chr', '')

    # Function to find the repeat type in a specific coordinate
    def find_repeat_type(chr_, pos, repeats_df):
        overlaps = repeats_df[(repeats_df['chr'] == chr_) & 
                              (repeats_df['end'] >= pos) & 
                              (repeats_df['start'] <= pos)]
        return ' '.join(overlaps['type'].unique())

    # Loop over each bed file to annotate and save results
    for bed_file in bed_files:
        # Read the bed file
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end'])
        
        # Annotate the bed file with repeat types
        bed_df['repeat_start'] = bed_df.apply(lambda row: find_repeat_type(row['chr'], row['start'], repeats_df), axis=1)
        bed_df['repeat_end'] = bed_df.apply(lambda row: find_repeat_type(row['chr'], row['end'], repeats_df), axis=1)

        # Determine if this is the circular DNA bed file
        if 'circularDNA' in bed_file:
            output_file_name = f"circularDNA_repeat_elements.bed"
        elif 'circularRNA' in bed_file:
            output_file_name = f"circularRNA_repeat_elements.bed"
        else:
            # Extract tool name and append true_or_false tag
            tool_name = os.path.basename(os.path.dirname(bed_file))
            if true_or_false == 'truepositives':
                basename ='TP'
            elif true_or_false == 'falsenegatives':
                basename ='FN'
            elif true_or_false == 'falsepositives':
                basename ='FP'
            
            output_file_name = f"{tool_name}_repeat_elements_{basename}.bed"
        
        # Create the output path with the tool name included
        output_path = os.path.join(output_dir, true_or_false, output_file_name)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Save the annotated bed file
        bed_df.to_csv(output_path, sep='\t', index=False, header=False)

def process_repeat_elements(annotated_files, tools, output_csv):
    """
    Process repeat elements from annotated BED files, counting occurrences of each repeat type across different tools.

    Parameters:
        annotated_files (list of str): List of paths to annotated BED files.
        tools (list of str): List of tool names to be included in the analysis.
        output_csv (str): Path to the output CSV file where results will be saved.

    Returns:
        None
    """
    # Define categories
    predefined_categories = {'SINE', 'LINE', 'LTR', 'DNA', 'Satellite', 'Ø'}
    
    # Dictionary to store counts of repeat elements
    repeat_counts = defaultdict(lambda: defaultdict(int))
    
    def combine_repeats(row):
        """
        Combine repeat start and end elements into a single string, handling empty values.
        
        Parameters:
            row (pd.Series): Row of the DataFrame containing repeat start and end elements.
        
        Returns:
            str: Combined and categorized repeat element string.
        """
        # Handle cases where either repeat_start or repeat_end is empty
        repeat_start = row['repeat_start'] if pd.notna(row['repeat_start']) and row['repeat_start'] != '' else 'Ø'
        repeat_end = row['repeat_end'] if pd.notna(row['repeat_end']) and row['repeat_end'] != '' else 'Ø'
        
        # Combine repeat elements and sort
        combined = '-'.join(sorted([repeat_start, repeat_end]))
        
        # Split combined string into parts
        combined_parts = combined.split('-')
        
        # Check if any part is not in predefined categories
        if any(part not in predefined_categories for part in combined_parts):
            combined_parts = [part for part in combined_parts if part in predefined_categories]
            combined = 'Other-' + '-'.join(combined_parts) if combined_parts else 'Other-Other'

        return combined
    
    for bed_file in annotated_files:
        # Extract tool name from the file path
        tool_name = os.path.basename(bed_file).split('_')[0]
        
        # Handle specific tool names for ecc_finder
        if 'ecc_finder' in bed_file:
            if 'bwa' in bed_file:
                tool_name = 'ecc_finder-bwa'
            elif 'minimap2' in bed_file:
                tool_name = 'ecc_finder-minimap2'
        elif 'Circle_finder' in bed_file:
            tool_name = 'Circle_finder'
        elif 'circRNA_finder' in bed_file:
            tool_name = 'circRNA_finder'
        elif 'find_circ' in bed_file:
            tool_name = 'find_circ'
        
        # Read the annotated bed file
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'repeat_start', 'repeat_end'])
        
        # Replace empty or NaN values in 'repeat_start' and 'repeat_end' with 'Ø'
        bed_df['repeat_start'] = bed_df['repeat_start'].replace('', 'Ø').fillna('Ø')
        bed_df['repeat_end'] = bed_df['repeat_end'].replace('', 'Ø').fillna('Ø')
        
        # Combine the repeat columns into a single repeat element
        bed_df['repeat_element'] = bed_df.apply(combine_repeats, axis=1)
        
        # Count occurrences of each repeat element
        repeat_counts[tool_name].update(bed_df['repeat_element'].value_counts().to_dict())
    
    # Create a unified DataFrame from the counts dictionary
    combined_counts = defaultdict(dict)
    for tool, counts in repeat_counts.items():
        for repeat_element, count in counts.items():
            combined_counts[repeat_element][tool] = count
    
    combined_df = pd.DataFrame.from_dict(combined_counts, orient='index').fillna(0).astype(int)
    
    # Ensure all tools are present
    all_tools = tools
    for tool in all_tools:
        if tool not in combined_df.columns:
            combined_df[tool] = 0
    
    # Add a 'Repeat Element' column
    combined_df = combined_df.reset_index().rename(columns={'index': 'Repeat Element'})
    
    # Ensure 'Ø-Ø' is included in the DataFrame
    if 'Ø-Ø' not in combined_df['Repeat Element'].values:
        # Create a DataFrame for 'Ø-Ø'
        zero_row = pd.DataFrame([[0] * len(combined_df.columns)], columns=combined_df.columns)
        zero_row['Repeat Element'] = 'Ø-Ø'
        
        # Append the zero_row DataFrame to the existing DataFrame
        combined_df = pd.concat([combined_df, zero_row], ignore_index=True)
    
    # Define the correct column order based on available columns
    column_order = ['Repeat Element'] + tools
    available_columns = [col for col in column_order if col in combined_df.columns]
    
    # Reorder columns and save to CSV
    combined_df = combined_df[available_columns]
    combined_df.to_csv(output_csv, index=False)


def calculate_repeat_element_metrics(tp_file, fn_file, fp_file, output_dir, tools):
    """
    Calculate precision, recall, and F-score for repeat elements across multiple tools.

    Parameters:
        tp_file (str): Path to the CSV file containing true positive counts.
        fn_file (str): Path to the CSV file containing false negative counts.
        fp_file (str): Path to the CSV file containing false positive counts.
        output_dir (str): Directory to save the output CSV files.
        tools (list of str): List of tool names to evaluate.

    Returns:
        None
    """
    # Read the CSV files
    tp_df = pd.read_csv(tp_file)
    fn_df = pd.read_csv(fn_file)
    fp_df = pd.read_csv(fp_file)

    # Initialize dataframes to store precision, recall, and fscore
    precision_df = pd.DataFrame(columns=['Repeat Element'] + tools)
    recall_df = pd.DataFrame(columns=['Repeat Element'] + tools)
    fscore_df = pd.DataFrame(columns=['Repeat Element'] + tools)

    # Set Repeat Element as the index for each dataframe
    tp_df.set_index('Repeat Element', inplace=True)
    fn_df.set_index('Repeat Element', inplace=True)
    fp_df.set_index('Repeat Element', inplace=True)

    # Fill NaN values with 0 for calculation purposes
    tp_df.fillna(0, inplace=True)
    fn_df.fillna(0, inplace=True)
    fp_df.fillna(0, inplace=True)

    # Calculate precision, recall, and fscore
    for tool in tools:
        precision_values = []
        recall_values = []
        fscore_values = []

        for repeat_element in tp_df.index:
            tp = tp_df.loc[repeat_element, tool].astype(int)
            fn = fn_df.loc[repeat_element, tool].astype(int)
            fp = fp_df.loc[repeat_element, tool].astype(int)

            precision, recall, fscore = calculate_metrics(tp, fp, fn)

            precision_values.append(precision)
            recall_values.append(recall)
            fscore_values.append(fscore)

        precision_df[tool] = precision_values
        recall_df[tool] = recall_values
        fscore_df[tool] = fscore_values

    # Add Repeat Element column to the dataframes
    precision_df['Repeat Element'] = tp_df.index
    recall_df['Repeat Element'] = tp_df.index
    fscore_df['Repeat Element'] = tp_df.index

    # Reorder columns to have Repeat Element first and remove 'circularDNA' column
    def filter_columns(df, tools_to_keep):
        # Reorder columns to include 'Repeat Element' and tools_to_keep
        columns = ['Repeat Element'] + tools_to_keep
        df = df[columns]
        return df

    tools_to_keep = tools[1:]  # Exclude the first tool, assuming 'circularDNA' is the first

    precision_df = filter_columns(precision_df, tools_to_keep)
    recall_df = filter_columns(recall_df, tools_to_keep)
    fscore_df = filter_columns(fscore_df, tools_to_keep)

    # Save the results to CSV files
    precision_file = f"{output_dir}/repeat_elements_precision.csv"
    recall_file = f"{output_dir}/repeat_elements_recall.csv"
    fscore_file = f"{output_dir}/repeat_elements_fscore.csv"

    precision_df.to_csv(precision_file, index=False)
    recall_df.to_csv(recall_file, index=False)
    fscore_df.to_csv(fscore_file, index=False)


def plot_stats_repeat_elements(csv_file, metric_name, output_dir):
    """
    Plots and saves two dot plots of repeat element statistics by tool:
    1. All repeat element combinations.
    2. Only identical repeat element combinations (e.g., LINE-LINE, SINE-SINE).

    Parameters:
        csv_file (str): Path to the CSV file containing the data.
        metric_name (str): Name of the metric to plot (e.g., 'Precision', 'Recall').
        output_dir (str): Directory to save the plots.

    Returns:
        None
    """
    # Load the CSV file
    df = pd.read_csv(csv_file)

    # Melt the dataframe for seaborn compatibility
    df_melted = df.melt(id_vars=['Repeat Element'], var_name='Tool', value_name=metric_name)

    # Define a colorblind-friendly palette
    colorblind_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']

    # --- PLOT 1: All Combinations ---
    fig, ax = plt.subplots(figsize=(12, 3.6))
    sns.stripplot(
        data=df_melted,
        x='Repeat Element',
        y=metric_name,
        hue='Tool',
        palette=colorblind_palette,
        jitter=True,  
        dodge=True,   
        s=8,          
        alpha=0.8,    
        ax=ax
    )

    # Shade every other x-axis label
    x_labels = df['Repeat Element'].unique()
    for i, label in enumerate(x_labels):
        if i % 2 == 0:  # Shade every other label
            ax.axvspan(i - 0.5, i + 0.5, color='lightgray', alpha=0.3)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.ylabel(metric_name.capitalize().replace("Fscore", "F-score"), fontsize=16)
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
    #plt.legend([], [], frameon=False)  
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/{metric_name.lower()}_all_combinations_dotplot.png", dpi=300)
    plt.show()
    plt.close()

    # --- PLOT 2: Only Identical Combinations ---
    identical_combinations = ['LINE-LINE', 'SINE-SINE', 'Satellite-Satellite', 'DNA-DNA', 'Other-Other', 'Ø-Ø']
    df_identical = df[df['Repeat Element'].isin(identical_combinations)]
    df_identical_melted = df_identical.melt(id_vars=['Repeat Element'], var_name='Tool', value_name=metric_name)
    df_identical_melted['Repeat Element'] = df_identical_melted['Repeat Element'].str.split('-').str[0]

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.stripplot(
        data=df_identical_melted,
        x='Repeat Element',
        y=metric_name,
        hue='Tool',
        palette=colorblind_palette,
        jitter=True,
        dodge=True,
        s=8,
        alpha=0.8,
        ax=ax
    )

    # Shade every other x-axis label
    x_labels = df_identical['Repeat Element'].unique()
    for i, label in enumerate(x_labels):
        if i % 2 == 0:
            ax.axvspan(i - 0.5, i + 0.5, color='lightgray', alpha=0.3)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.ylabel(metric_name.capitalize().replace("Fscore", "F-score"), fontsize=16)
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
    #plt.legend([], [], frameon=False)  
    plt.tight_layout()#rect=[0, 0, 0.85, 1])  
    plt.savefig(f"{output_dir}/{metric_name.lower()}_identical_combinations_dotplot.png", dpi=300)
    plt.show()
    plt.close()