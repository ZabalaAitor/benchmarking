import os
import csv
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from functions.compare_circular_data import calculate_metrics

gene_file = "genomic_elements/genomic_genes.gtf"
exon_file = "genomic_elements/exon_genes.gtf"
other_file = "genomic_elements/genomic_annotation.gtf"

def annotate_genomic_elements(bed_file, output_dir, gene_file, exon_file, other_file):
    """
    Annotates genomic elements in BED files using provided annotation files.

    Parameters:
        bed_file (str): Path to the BED file to annotate.
        output_dir (str): Path to the output annotated BED file.
        gene_file (str): Path to the file with gene annotations.
        exon_file (str): Path to the file with exon annotations.
        other_file (str): Path to the file with other annotations.

    Returns:
        None
    """
    # Load gene positions
    gene_data = {}
    with open(gene_file, "r") as gene_file:
        for line in gene_file:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chromosome = cols[0]
            start_pos = int(cols[3])
            end_pos = int(cols[4])
            if chromosome not in gene_data:
                gene_data[chromosome] = []
            gene_data[chromosome].append((start_pos, end_pos))

    # Load exon positions
    exon_data = {}
    with open(exon_file, "r") as exon_file:
        for line in exon_file:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chromosome = cols[0]
            start_pos = int(cols[3])
            end_pos = int(cols[4])
            if chromosome not in exon_data:
                exon_data[chromosome] = []
            exon_data[chromosome].append((start_pos, end_pos))

    # Load other positions with annotations
    other_data = {}
    with open(other_file, "r") as other_file:
        for line in other_file:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chromosome = cols[0]
            start_pos = int(cols[3])
            end_pos = int(cols[4])
            annotation = cols[2]  # Get the annotation from the second column
            if chromosome not in other_data:
                other_data[chromosome] = []
            other_data[chromosome].append((start_pos, end_pos, annotation))

    # Annotate BSJ bed file
    with open(bed_file, "r") as bsj_file, open(output_dir, "w") as annotated_file:
        writer = csv.writer(annotated_file, delimiter='\t')
        for line in bsj_file:
            cols = line.strip().split("\t")
            chromosome = cols[0].replace('chr', '')
            start_position = int(cols[1])
            end_position = int(cols[2])
            start_annotation = "intergenic"  # Default annotation is intergenic
            end_annotation = "intergenic"  # Default annotation is intergenic

            # Function to annotate a single position
            def annotate_position(chromosome, position):
                # Check if position is within other regions with annotations
                if chromosome in other_data:
                    for other_start, other_end, other_annotation in other_data[chromosome]:
                        if other_start <= position <= other_end:
                            return other_annotation

                # Check if position is within exon regions
                if chromosome in exon_data:
                    for exon_start, exon_end in exon_data[chromosome]:
                        if exon_start <= position <= exon_end:
                            return "exon"

                # Check if position is within gene regions
                if chromosome in gene_data:
                    for gene_start, gene_end in gene_data[chromosome]:
                        if gene_start <= position <= gene_end:
                            return "intron"

                return "intergenic"

            # Annotate both start and end positions
            start_annotation = annotate_position(chromosome, start_position)
            end_annotation = annotate_position(chromosome, end_position)

            # Write annotated BSJ to file
            writer.writerow([chromosome, start_position, end_position, start_annotation, end_annotation])

def annotate_bed_files_genomic_elements(bed_files, output_dir, tools, true_or_false):
    """
    Annotates a list of BED files with genomic elements and saves the annotated files.

    Parameters:
        bed_files (list): List of paths to BED files to annotate.
        output_dir (str): Directory to save annotated files.
        tools (list): List of tool names used to determine output file names.
        true_or_false (str): Label for the data type ('truepositives', 'falsenegatives', 'falsepositives').
        
    Returns:
        None
    """ 
    # Determine basename for output files
    if true_or_false == 'truepositives':
        basename = 'TP'
    elif true_or_false == 'falsenegatives':
        basename = 'FN'
    elif true_or_false == 'falsepositives':
        basename = 'FP'
    
    # Loop over each BED file to annotate
    for bed_file in bed_files:
        # Determine tool from the file path
        tool_name = os.path.basename(os.path.dirname(bed_file))
        # Create the output file path
        annotated_bsj_file = os.path.join(output_dir, true_or_false, f'{tool_name}_genomic_elements_{basename}.bed')
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(annotated_bsj_file), exist_ok=True)
        
        # Annotate the BED file
        annotate_genomic_elements(bed_file, annotated_bsj_file, gene_file, exon_file, other_file)

def process_genomic_elements(annotated_files, output_csv, tools):
    """
    Process annotated BED files to count genomic element combinations (genomic_start and genomic_end) 
    across multiple tools and aggregate the results into a CSV file.

    Parameters:
        annotated_files (list of str): List of file paths to annotated BED files. Each file is expected 
                                       to contain five columns: chr, start, end, genomic_start, genomic_end.
        output_csv (str): Path to the output CSV file where the combined count matrix will be saved.
        tools (list of str): List of tool identifiers used to map file names to their corresponding tool.

    Returns:
        None: The function saves the output as a CSV file and does not return any value.
    """
    # Dictionary to store counts of genomic elements
    genomic_counts = defaultdict(lambda: defaultdict(int))
    
    # Use the tools list to generate the identifiers
    identifiers = tools
    
    def get_tool_name(file_path):
        base_name = os.path.basename(file_path)
        # Check for presence of each identifier in the file name
        for identifier in identifiers:
            if identifier in base_name:
                return identifier
        return 'unknown'  # Default value if no identifier matches
    
    def combine_genomics(row):
        # Handle cases where either genomic_start or genomic_end is empty
        genomic_start = row['genomic_start'] if pd.notna(row['genomic_start']) and row['genomic_start'] != '' else 'Ø'
        genomic_end = row['genomic_end'] if pd.notna(row['genomic_end']) and row['genomic_end'] != '' else 'Ø'
        
        # Combine genomic elements and sort
        combined = '-'.join(sorted([genomic_start, genomic_end]))
        
        return combined
    
    for bed_file in annotated_files:
        # Determine tool name based on file name
        tool_name = get_tool_name(bed_file)
        
        # Read the annotated bed file
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'genomic_start', 'genomic_end'])
        
        # Replace empty or NaN values in 'genomic_start' and 'genomic_end' with 'Ø'
        bed_df['genomic_start'] = bed_df['genomic_start'].replace('', 'Ø').fillna('Ø')
        bed_df['genomic_end'] = bed_df['genomic_end'].replace('', 'Ø').fillna('Ø')
        
        # Combine the genomic columns into a single genomic element
        bed_df['genomic_element'] = bed_df.apply(combine_genomics, axis=1)
        
        # Count occurrences of each genomic element
        genomic_counts[tool_name].update(bed_df['genomic_element'].value_counts().to_dict())
    
    # Create a unified DataFrame from the counts dictionary
    combined_counts = defaultdict(dict)
    for tool, counts in genomic_counts.items():
        for genomic_element, count in counts.items():
            combined_counts[genomic_element][tool] = count
    
    combined_df = pd.DataFrame.from_dict(combined_counts, orient='index').fillna(0).astype(int)
    
    # Ensure all tools are present
    for tool in tools:
        if tool not in combined_df.columns:
            combined_df[tool] = 0
    
    # Add a 'Genomic Element' column
    combined_df = combined_df.reset_index().rename(columns={'index': 'Genomic Element'})
    
    # Define the correct column order based on available columns
    column_order = ['Genomic Element'] + tools
    available_columns = [col for col in column_order if col in combined_df.columns]
    
    # Reorder columns and save to CSV
    combined_df = combined_df[available_columns]
    combined_df.to_csv(output_csv, index=False)

def calculate_genomic_element_metrics(tp_file, fn_file, fp_file, output_dir, tools):
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

    # Set Genomic Element as the index for each dataframe
    tp_df.set_index('Genomic Element', inplace=True)
    fn_df.set_index('Genomic Element', inplace=True)
    fp_df.set_index('Genomic Element', inplace=True)

    # Combine all unique genomic elements from all dataframes
    all_genomic_elements = tp_df.index.union(fn_df.index).union(fp_df.index)

    # Reindex the dataframes to ensure all genomic elements are accounted for
    tp_df = tp_df.reindex(all_genomic_elements, fill_value=0)
    fn_df = fn_df.reindex(all_genomic_elements, fill_value=0)
    fp_df = fp_df.reindex(all_genomic_elements, fill_value=0)

    # Initialize dataframes to store precision, recall, and fscore
    precision_df = pd.DataFrame(columns=['Genomic Element'] + tools)
    recall_df = pd.DataFrame(columns=['Genomic Element'] + tools)
    fscore_df = pd.DataFrame(columns=['Genomic Element'] + tools)

    # Calculate precision, recall, and fscore
    for tool in tools:
        precision_values = []
        recall_values = []
        fscore_values = []

        for genomic_element in all_genomic_elements:
            tp = tp_df.loc[genomic_element, tool].astype(int)
            fn = fn_df.loc[genomic_element, tool].astype(int)
            fp = fp_df.loc[genomic_element, tool].astype(int)

            precision, recall, fscore = calculate_metrics(tp, fp, fn)

            precision_values.append(precision)
            recall_values.append(recall)
            fscore_values.append(fscore)

        precision_df[tool] = precision_values
        recall_df[tool] = recall_values
        fscore_df[tool] = fscore_values

    # Add Genomic Element column to the dataframes
    precision_df['Genomic Element'] = all_genomic_elements
    recall_df['Genomic Element'] = all_genomic_elements
    fscore_df['Genomic Element'] = all_genomic_elements

    # Save the results to CSV files
    precision_file = f"{output_dir}/genomic_elements_precision.csv"
    recall_file = f"{output_dir}/genomic_elements_recall.csv"
    fscore_file = f"{output_dir}/genomic_elements_fscore.csv"

    precision_df.to_csv(precision_file, index=False)
    recall_df.to_csv(recall_file, index=False)
    fscore_df.to_csv(fscore_file, index=False)

def plot_stats_genomic_elements(csv_file, metric_name, output_dir):
    """
    Generate and save dot plots to visualize genomic element statistics for different tools.

    This function creates two dot plots based on a CSV file containing metric values 
    (e.g., precision, recall, F-score) for genomic element combinations:
        1. A dot plot including all combinations of genomic elements (e.g., LINE-SINE, LINE-LTR).
        2. A dot plot specifically highlighting identical element combinations (e.g., LINE-LINE, SINE-SINE).

    Parameters:
        csv_file (str): Path to the input CSV file. It must include a column named 'Genomic Element' and
                        additional columns for each tool with corresponding metric values.
        metric_name (str): The name of the metric to plot (e.g., 'precision', 'recall', 'Fscore').
                           This is used for y-axis labeling and file naming.
        output_dir (str): Directory where the output plot image will be saved.

    Returns:
        None: The function saves plot(s) to the specified directory and displays them, but does not return any value.
    """
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Melt the dataframe for seaborn compatibility
    df_melted = df.melt(id_vars=['Genomic Element'], var_name='Tool', value_name=metric_name)
    
    # Define a colorblind-friendly palette
    custom_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']
    
    # --- PLOT 1: All Combinations ---
    plt.figure(figsize=(8, 4.3))
    ax = sns.stripplot(
        data=df_melted,
        x='Genomic Element',
        y=metric_name,
        hue='Tool',
        palette=custom_palette,
        jitter=True,
        dodge=True,
        s=8,
        alpha=0.8
    )
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.ylabel(metric_name.capitalize().replace("Fscore", "F-score"), fontsize=16)
    plt.xlabel('')
    
    # Adjust x-axis labels
    new_labels = [label.replace("five_prime_utr", "5'-UTR").replace("three_prime_utr", "3'-UTR") for label in df_melted['Genomic Element'].unique()]

    # Shade every other x-axis label
    x_labels = df['Genomic Element'].unique()
    for i, label in enumerate(x_labels):
        if i % 2 == 0:  # Shade every other label
            ax.axvspan(i - 0.5, i + 0.5, color='lightgray', alpha=0.3)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.ylabel(metric_name.capitalize().replace("Fscore", "F-score"), fontsize=16)
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=16)
    ax.set_xticklabels(new_labels)
    plt.yticks(fontsize=16)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)
    #plt.legend([], [], frameon=False)  
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/{metric_name.lower()}_all_combinations_dotplot.png", dpi=300)
    plt.show()
    plt.close()


    # --- PLOT 2: Only Identical Combinations ---

    # Define order based on genomic frequency
    genomic_order = ["intergenic", "intron", "exon", "3'-UTR", "5'-UTR"]
    
    identical_combinations = ['intron-intron', 'intergenic-intergenic', 'exon-exon', 'three_prime_utr-three_prime_utr', 'five_prime_utr-five_prime_utr']
    df_identical = df[df['Genomic Element'].isin(identical_combinations)]

    df_identical_melted = df_identical.melt(id_vars=['Genomic Element'], var_name='Tool', value_name=metric_name)
    df_identical_melted['Genomic Element'] = df_identical_melted['Genomic Element'].str.split('-').str[0]  # Simplify names
    
    # Replace UTR names and filter only existing elements
    df_identical_melted['Genomic Element'] = df_identical_melted['Genomic Element'].replace({
        "three_prime_utr": "3'-UTR",
        "five_prime_utr": "5'-UTR"
    })
    present_genomic_elements_identical = [g for g in genomic_order if g in df_identical_melted['Genomic Element'].unique()]

    plt.figure(figsize=(8, 4))
    ax = sns.stripplot(
        data=df_identical_melted,
        x='Genomic Element',
        y=metric_name,
        hue='Tool',
        palette=custom_palette,
        jitter=True,
        dodge=True,
        s=8,
        alpha=0.8,
        order=present_genomic_elements_identical  # Apply order dynamically
    )
  
    # Shade every other x-axis label
    x_labels = df_identical['Genomic Element'].unique()
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
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{metric_name.lower()}_identical_combinations_dotplot.png", dpi=300)
    plt.show()
    plt.close()
