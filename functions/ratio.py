import os
import pandas as pd
import pysam
from scipy.stats import binom, kruskal, median_abs_deviation
from statsmodels.stats.multitest import multipletests as multi
from matplotlib import pyplot as plt
import seaborn as sns
import scikit_posthocs as sp
import numpy as np
import matplotlib.gridspec as gridspec


def process_circle_matrix(matrix_file, bam_file, output_file, use_chr_prefix=False, N_offset=1, verbose=False):
    """
    Processes a circle matrix and extracts read counts from a BAM file.

    Parameters:
        matrix_file (str): Path to the input CSV containing circle regions.
        bam_file (str): Path to the BAM file.
        output_file (str): Path to save the output CSV.
        use_chr_prefix (bool): Whether to add 'chr' prefix to chromosome names. Default is False.

    Returns:
        None
    """
    # Open BAM file
    bamobject = pysam.AlignmentFile(bam_file, "rb")

    # Read the CSV
    df = pd.read_csv(matrix_file)

    # Prepare results
    results = []

    for index, row in df.iterrows():
        circle = row["Circles"]
        chrom, coords = circle.split(":")
        start, end = map(int, coords.split("-"))

        # Add 'chr' prefix if needed
        if use_chr_prefix:
            chrom = 'chr' + chrom

        try:
            # Count how many tools detected the circle (how many 1s)
            tools_detected = row.iloc[1:].sum()

            # Fetch CJ1 Reads (start-1, start, start+1)
            CJ1 = {read.query_name for read in bamobject.fetch(chrom, start - N_offset, start + N_offset)}
            CJ1_n = len(set(CJ1))

            # display(CJ1)

            if verbose:
                print('CJ1', CJ1_n)
                for read in bamobject.fetch(chrom, start - N_offset, start + N_offset):
                    read_str = read.query_name
                    rstart_real = int(read_str.split('|')[2].split('-')[0])
                    rend_real = int(read_str.split('|')[2].split('-')[1])
                    rstart_STAR = int(read.reference_start)
                    rend_STAR = int(read.reference_end)

                    print(read_str.split('|')[0], read.cigarstring, rstart_real, rend_real, rstart_STAR, rend_STAR, rstart_real - rstart_STAR, rend_real - rend_STAR)


            # Fetch CJ2 Reads (end-1, end, end+1)
            CJ2 = {read.query_name for read in bamobject.fetch(chrom, end - N_offset, end + N_offset)}
            CJ2_n = len(set(CJ2))

            if verbose:
                print('CJ2', CJ2_n)
                for read in bamobject.fetch(chrom, end - N_offset, end + N_offset):
                    read_str = read.query_name
                    rstart_real = int(read_str.split('|')[2].split('-')[0])
                    rend_real = int(read_str.split('|')[2].split('-')[1])
                    rstart_STAR = int(read.reference_start)
                    rend_STAR = int(read.reference_end)

                    print(read_str.split('|')[0], read.cigarstring, rstart_real, rend_real, rstart_STAR, rend_STAR, rstart_real - rstart_STAR, rend_real - rend_STAR)


            # Sum of CJ1 and CJ2
            CJ_unique_reads = len(CJ1 | CJ2)
            CJ_reads = CJ1_n + CJ2_n

            # Fetch total reads in the region (start-1 to end+1)
            total_reads = len({read.query_name for read in bamobject.fetch(chrom, start - N_offset, end + N_offset)})

            # Compute ratios
            ratio_CJ1 = CJ1_n / CJ_reads if CJ_reads > 0 else 0
            ratio_CJ2 = CJ2_n / CJ_reads if CJ_reads > 0 else 0

            # Determine max and min CJ
            max_CJ = max(ratio_CJ1, ratio_CJ2)
            min_CJ = min(ratio_CJ1, ratio_CJ2)

            # Compute CJ max and min difference
            diff_CJ = max_CJ - min_CJ

            # Compute diff CJ probability
            min_reads_CJ = min(CJ1_n, CJ2_n)
            p_diff_cj = binom.cdf(min_reads_CJ, CJ_reads, p=0.5) * 2
            if p_diff_cj > 1: # Adjust for cases where n is even and CJ reads distribution is equal
                p_diff_cj = 1

            # Compute additional ratios
            ratio = CJ_unique_reads / total_reads if total_reads > 0 else 0
            circle_length = end - start + 1

            # Store results
            results.append([circle, circle_length, tools_detected, CJ1_n, CJ2_n, CJ_unique_reads, CJ_reads,  min_CJ, max_CJ, diff_CJ, p_diff_cj, total_reads, ratio])

        except ValueError as e:
            # Handle the error if there's an invalid chromosome or other issue
            print(f"Skipping invalid chromosome: {chrom} or region {circle}. Error: {e}")
            continue

    # Close BAM file
    bamobject.close()

    # Create output DataFrame
    output_df = pd.DataFrame(results, columns=[
        "Circle", "Length", "Tools", "CJ1 Reads", "CJ2 Reads", "CJ Unique Reads", "CJ Reads", "min CJ", "max CJ", "diff CJ", "p diff CJ", "Total Reads", "Ratio"
    ])
    
    # Compute adjusted p-values
    adjusted_p_CJ = multi(output_df['p diff CJ'].values, method='fdr_bh')[1]

    # Insert the new column right after 'p_diff_CJ'
    p_diff_index = output_df.columns.get_loc('p diff CJ')
    output_df.insert(p_diff_index + 1, 'p diff CJ*', adjusted_p_CJ)

    # Save to CSV
    output_df.to_csv(output_file, index=False)


def circle_diff(matrix_dir, matrix_reads, tools, output_dir, group):
    # Define file paths
    matrix_path = os.path.join(matrix_dir)
    matrix_with_reads_path = os.path.join(matrix_reads)

    # Load data
    matrix_df = pd.read_csv(matrix_path)
    reads_df = pd.read_csv(matrix_with_reads_path)

    # Merge on circle ID
    merged_df = pd.merge(matrix_df, reads_df[['Circle', 'CJ1 Reads', 'CJ2 Reads', 'CJ Reads', 'diff CJ', 'p diff CJ', 'p diff CJ*']], left_on='Circles', right_on='Circle', how='inner')

    n_threshold = 9
    alpha = 0.05

    summary_stats = []  # List to store summary statistics
    significance_stats = []  # List to store significance statistics

    for tool in tools:
        tool_data = merged_df[merged_df[tool] == 1]
        
        filter_data_max = tool_data[tool_data['CJ Reads'] >= n_threshold]
        filter_data_min = tool_data[tool_data['CJ Reads'] < n_threshold]
 
        # Collecting the statistics
        summary_stats.append({
            'Tool': tool,
            'Percentage CJ Reads >= 9': len(filter_data_max) / len(tool_data),
            'Median ACJ (>=9)': filter_data_max['diff CJ'].median(),
            'Mean ACJ (>=9)': filter_data_max['diff CJ'].mean(),
            'Median ACJ (<9)': filter_data_min['diff CJ'].median(),
            'Mean ACJ (<9)': filter_data_min['diff CJ'].mean(),
            'Ratio p diff CJ < alpha (>=9)': len(filter_data_max[filter_data_max['p diff CJ'] < alpha]) / len(filter_data_max),
            'Ratio p diff CJ* < alpha (>=9)': len(filter_data_max[filter_data_max['p diff CJ*'] < alpha]) / len(filter_data_max),
            'Ratio p diff CJ < alpha (<9)': len(filter_data_min[filter_data_min['p diff CJ'] < alpha]) / len(filter_data_min),
            'Ratio p diff CJ* < alpha (<9)': len(filter_data_min[filter_data_min['p diff CJ*'] < alpha]) / len(filter_data_min),
        })

        # if tool == "Simulated":
        #     circles_cabrones = filter_data_max[filter_data_max['p diff CJ*'] < alpha]
        #     display(circles_cabrones)

    # Save summary statistics as CSV
    summary_df = pd.DataFrame(summary_stats)
    
    # Define output path
    output_csv_path = os.path.join(output_dir, 'statistics_summary.csv')
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)

    # Save the summary DataFrame to CSV
    summary_df.to_csv(output_csv_path)

    # Prepare data for plotting
    plot_data = []
    for tool in tools:
        tool_data = merged_df[merged_df[tool] == 1]
        for _, row in tool_data.iterrows():
            plot_data.append({
                'Tool': tool,
                'diff CJ': row['diff CJ']
            })

    # Define color palette
    palette = {
        'Simulated': '#b3b3b3',
        'unfilter': '#d46014',
        'filter-split': '#ddcd3d',
        'filter-duplicates': '#064b76ff',
        'filter': '#63bdf6ff'
    }

    # Create DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)
    plot_df
    # Plot
    plt.figure(figsize=(8.5, 4))
    ax = sns.boxplot(x="Tool", y="diff CJ", data=plot_df, palette=palette)

    # Axis labels and ticks
    plt.xlabel('')
    plt.ylabel('ΔCJ', fontsize=16)
    xtick_labels = tools
    plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.tight_layout()

    # Add counts on top of each box
    group_counts = plot_df.groupby('Tool').size()
    global_ymax = plot_df['diff CJ'].max()

    for i, tool in enumerate(tools):
        count = group_counts.get(tool, 0)
        if count > 0:
            y_text = 1.05
            ax.text(i, y_text, str(count), ha='center', va='bottom', fontsize=12)


    # Save the plot
    save_path = os.path.join(output_dir, f'{group}_diff_cj_plot.png')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    # Group data by tool for Kruskal-Wallis test
    grouped = [group["diff CJ"].values for _, group in plot_df.groupby("Tool") if not group.empty]

    # Perform Kruskal-Wallis test
    if len(grouped) > 1:
        stat, p_val = kruskal(*grouped)
        kruskal_results = [{
            'Test': 'Kruskal-Wallis (across Tools)',
            'H-statistic': stat,
            'p-value': p_val
        }]
        print(f"Kruskal-Wallis test across tools: H = {stat:.4f}, p = {p_val:.4e}")
    else:
        print("Not enough groups to perform Kruskal-Wallis test.")

def circle_diff_real(matrix_dir, tools, filtering_methods, data, n_threshold=9, alpha=0.05):

    plot_data = []
    kruskal_results = []
    dunn_results = []
    summary_stats = []
    significance_stats = []

    def load_matrix_data(filtering_method, tool):
        matrix_path = os.path.join(matrix_dir, filtering_method, data, 'matrix.csv')
        matrix_with_reads_path = os.path.join(matrix_dir, filtering_method, data, 'matrix_with_reads.csv')

        if not os.path.exists(matrix_path) or not os.path.exists(matrix_with_reads_path):
            print(f"File(s) missing for {filtering_method}: {matrix_path} or {matrix_with_reads_path}")
            return

        matrix_data = pd.read_csv(matrix_path)
        matrix_with_reads = pd.read_csv(matrix_with_reads_path)

        merged_data = pd.merge(
            matrix_data,
            matrix_with_reads[['Circle', 'CJ Reads', 'diff CJ', 'p diff CJ', 'p diff CJ*']],
            left_on='Circles',
            right_on='Circle',
            how='inner'
        )

        tool_data = merged_data[merged_data[tool] == 1]
        filter_data_max = tool_data[tool_data['CJ Reads'] >= n_threshold]
        filter_data_min = tool_data[tool_data['CJ Reads'] < n_threshold]

        percent_max = len(filter_data_max) / len(tool_data) if len(tool_data) > 0 else np.nan

        summary_stats.append({
            'Tool': tool,
            'Filtering': filtering_method,
            'CJ ≥ threshold (%)': f"{percent_max:.2%}",
            'CJ ≥ threshold (n)': len(filter_data_max),
            'CJ < threshold (n)': len(filter_data_min),
            'Median ΔCJ (≥)': filter_data_max['diff CJ'].median(),
            'Mean ΔCJ (≥)': filter_data_max['diff CJ'].mean(),
            'Median ΔCJ (<)': filter_data_min['diff CJ'].median(),
            'Mean ΔCJ (<)': filter_data_min['diff CJ'].mean(),
            f'Ratio p diff CJ < {alpha} (≥)': (filter_data_max['p diff CJ'] < alpha).mean(),
            f'Count p diff CJ < {alpha} (≥)': (filter_data_max['p diff CJ'] < alpha).sum(),
            f'Ratio p diff CJ < {alpha} (<)': (filter_data_min['p diff CJ'] < alpha).mean(),
            f'Count p diff CJ < {alpha} (<)': (filter_data_min['p diff CJ'] < alpha).sum(),
            f'Ratio p diff CJ* < {alpha} (≥)': (filter_data_max['p diff CJ*'] < alpha).mean(),
            f'Count p diff CJ* < {alpha} (≥)': (filter_data_max['p diff CJ*'] < alpha).sum(),
            f'Ratio p diff CJ* < {alpha} (<)': (filter_data_min['p diff CJ*'] < alpha).mean(),
            f'Count p diff CJ* < {alpha} (<)': (filter_data_min['p diff CJ*'] < alpha).sum()
        })

        for _, row in filter_data_max.iterrows():
            plot_data.append({
                'Tool': tool,
                'Filtering': filtering_method,
                'Circle': row['Circles'],
                'diff CJ': row['diff CJ']
            })

    for tool in tools:
        for filtering_method in filtering_methods:
            load_matrix_data(filtering_method, tool)

    # Save summary stats once after all iterations
    summary_save_path = os.path.join(matrix_dir, data, 'statistics_summary.csv')
    os.makedirs(os.path.dirname(summary_save_path), exist_ok=True)
    pd.DataFrame(summary_stats).to_csv(summary_save_path, index=False)

    plot_df = pd.DataFrame(plot_data)

    if not plot_df.empty:
        palette = {
            'unfilter': '#d46014',
            'filter-split': '#ddcd3d',
            'filter-duplicates': '#064b76ff',
            'filter': '#63bdf6ff'
        }

        plt.figure(figsize=(10, 4))
        ax = sns.boxplot(x="Tool", y="diff CJ", hue="Filtering", data=plot_df, palette=palette)
        plt.xlabel('')
        plt.ylabel('ΔCJ', fontsize=16)
        if data == 'Circle-Seq' or data == 'ATAC-seq':
            plt.xticks(ticks=range(len(tools)), labels=['CIRCexplorer2', 'Circle-Map', 'Circle_finder', 'ecc_finder\nbwa', 'ecc_finder\nminimap2'], fontsize=16)
        else:
            plt.xticks(ticks=range(len(tools)), labels=tools, fontsize=16)
        plt.yticks(fontsize=16)
        ax.legend_.remove()
        plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
        sns.despine()
        plt.tight_layout()

        group_counts = plot_df.groupby(['Tool', 'Filtering']).size().reset_index(name='count')
        for i, (tool, filtering) in enumerate(group_counts[['Tool', 'Filtering']].values):
            count = group_counts.loc[
                (group_counts['Tool'] == tool) & (group_counts['Filtering'] == filtering),
                'count'
            ].values[0]
            tool_idx = tools.index(tool)
            filtering_idx = list(palette.keys()).index(filtering)
            total_hue = len(palette)
            spacing_factor = 0.85
            x_position = tool_idx + ((filtering_idx - total_hue / 2) * spacing_factor / total_hue) + 0.1
            ax.text(x_position, 1.02, f'{count}', ha='center', va='bottom', fontsize=14, rotation=45)

        save_path = os.path.join(matrix_dir, data, 'diff CJ.png')
        os.makedirs(os.path.join(matrix_dir, data), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.show()

        for filtering_method in filtering_methods:
            df_sub = plot_df[plot_df['Filtering'] == filtering_method]
            if not df_sub.empty:
                groups = [df_sub[df_sub['Tool'] == tool]['diff CJ'].values for tool in tools]
                if len(groups) > 1:
                    stat, p_val = kruskal(*groups)
                    kruskal_results.append({
                        'Filtering Method': filtering_method,
                        'Test': 'Kruskal-Wallis (Tools)',
                        'H-statistic': stat,
                        'p-value': p_val
                    })
                    print(f"Kruskal-Wallis test for filtering method {filtering_method}: H = {stat:.4f}, p = {p_val:.4f}")

                    dunn_result = sp.posthoc_dunn(df_sub, val_col='diff CJ', group_col='Tool', p_adjust='bonferroni')
                    dunn_results.append({
                        'Filtering Method': filtering_method,
                        'Test': 'Dunn\'s test (Tools)',
                        'results': dunn_result
                    })
                    print(f"  Dunn's test for {filtering_method}: \n{dunn_result}")


        output_dir = os.path.join(matrix_dir, data)
        os.makedirs(output_dir, exist_ok=True)

        output_excel_path = os.path.join(output_dir, 'statistical_results.xlsx')
        with pd.ExcelWriter(output_excel_path) as writer:
            # Save Kruskal-Wallis results
            pd.DataFrame(kruskal_results).to_excel(writer, sheet_name='Kruskal-Wallis', index=False)
            
            # Save Dunn's test results
            for result in dunn_results:
                sheet_name = f"Dunn_{result['Filtering Method']}"
                result['results'].to_excel(writer, sheet_name=sheet_name[:31], index=True)  # Excel sheet name limit = 31 chars


    else:
        print("Plotting skipped due to empty DataFrame.")   


palette = {
    'unfilter': '#d46014',
    'filter-split': '#ddcd3d',
    'filter-duplicates': '#064b76ff',
    'filter': '#63bdf6ff'
}

def capitalize_xticks(ax):
    """Helper function to capitalize x-axis tick labels."""
    labels = [label.get_text().capitalize() for label in ax.get_xticklabels()]
    ax.set_xticklabels(labels)

def diff_cj_combinations(circle, method, ros_combinations_mapping, filtering_methods, combining_method):
    data = []
    differences_data = []
    rosette_diff_cj = {}
    output_dir = f"results/{circle}/real/{method}/"
    os.makedirs(output_dir, exist_ok=True)

    n_threshold = 9
    alpha = 0.05

    for file_key, other in ros_combinations_mapping.items():
        for filter_method in filtering_methods:
            for combination in combining_method:
                file_path = f"results/{circle}/real/{filter_method}/{method}/{combination}/{file_key}.bed"
                if not os.path.exists(file_path) or os.stat(file_path).st_size == 0:
                    print(f"File {file_path} does not exist or is empty.")
                    continue
                try:
                    df = pd.read_csv(file_path, sep=',')
                    if 'diff CJ' in df.columns and 'CJ Reads' in df.columns:
                        tool_data = df
                        filter_data_max = df[df['CJ Reads'] >= n_threshold]
                        filter_data_min = df[df['CJ Reads'] < n_threshold]

                        if df.empty:
                            continue

                        diff_CJ = filter_data_max['diff CJ'].mean()

                        if combination == "rosette":
                            rosette_diff_cj[(file_key, filter_method)] = diff_CJ
                        percent_max = len(filter_data_max) / len(tool_data) if len(tool_data) > 0 else np.nan
                        data.append({
                            'Combination': combination,
                            'Key': file_key,
                            'Filtering': filter_method,
                            'CJ ≥ threshold (%)': f"{percent_max:.2%}",
                            'CJ ≥ threshold (n)': len(filter_data_max),
                            'CJ < threshold (n)': len(filter_data_min),
                            'Median ΔCJ (≥)': filter_data_max['diff CJ'].median(),
                            'Mean ΔCJ (≥)': filter_data_max['diff CJ'].mean(),
                            'Median ΔCJ (<)': filter_data_min['diff CJ'].median(),
                            'Mean ΔCJ (<)': filter_data_min['diff CJ'].mean(),
                            f'Ratio p diff CJ < {alpha} (≥)': (filter_data_max['p diff CJ'] < alpha).mean(),
                            f'Count p diff CJ < {alpha} (≥)': (filter_data_max['p diff CJ'] < alpha).sum(),
                            f'Ratio p diff CJ < {alpha} (<)': (filter_data_min['p diff CJ'] < alpha).mean(),
                            f'Count p diff CJ < {alpha} (<)': (filter_data_min['p diff CJ'] < alpha).sum(),
                            f'Ratio p diff CJ* < {alpha} (≥)': (filter_data_max['p diff CJ*'] < alpha).mean(),
                            f'Count p diff CJ* < {alpha} (≥)': (filter_data_max['p diff CJ*'] < alpha).sum(),
                            f'Ratio p diff CJ* < {alpha} (<)': (filter_data_min['p diff CJ*'] < alpha).mean(),
                            f'Count p diff CJ* < {alpha} (<)': (filter_data_min['p diff CJ*'] < alpha).sum()
                        })

                        if combination != "rosette" and (file_key, filter_method) in rosette_diff_cj:
                            rosette_mean = rosette_diff_cj[(file_key, filter_method)]
                            diff_cj_difference = rosette_mean - diff_CJ
                            differences_data.append({
                                'Key': file_key,
                                'Combination': combination,
                                'Filtering': filter_method,
                                'ΔCJ Difference': diff_cj_difference
                            }) 

                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
                    
    # Convert data lists to DataFrames
    data_df = pd.DataFrame(data)

    differences_df = pd.DataFrame(differences_data)

    # === PLOT 1: Overall Diff CJ ===
    plt.figure(figsize=(7, 4))
    ax = sns.stripplot(data=data_df, x="Combination", y="Mean ΔCJ (≥)", hue="Filtering", jitter=True, dodge=True, palette=palette)
    plt.xlabel("")
    plt.ylabel("ΔCJ", fontsize=16)
    capitalize_xticks(ax)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    sns.despine()
    plt.legend([], frameon=False)
    plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)

    filters = list(palette.keys())
    for m_index, m in enumerate(combining_method):
        for f_index, f in enumerate(filters):
            x_pos = m_index + (f_index - len(filters) / 2) * 0.2 + 0.1
            df_sub = data_df[(data_df['Combination'] == m) & (data_df['Filtering'] == f)]
            if not df_sub.empty:
                median = np.median(df_sub["Mean ΔCJ (≥)"].values)
                MAD = median_abs_deviation(df_sub["Mean ΔCJ (≥)"].values)
                w = 0.07
                plt.plot([x_pos - w, x_pos + w], [median, median], color="#232323", zorder=10)
                plt.plot([x_pos + w, x_pos + w], [median + MAD, median - MAD], color="#bcbcbc", zorder=10)


    plt.tight_layout()
    plt.savefig(f"{output_dir}/diffCJ_overall.png", dpi=300)
    plt.show()

    # === PLOT 2: Difference from Rosette ===
    plt.figure(figsize=(7, 4))
    df_all = differences_df[differences_df['Combination'].isin(['union', 'intersect', 'double', 'unique'])]
    ax = sns.stripplot(data=df_all, x="Combination", y="ΔCJ Difference", hue="Filtering", jitter=True, dodge=True, palette=palette)
    plt.xlabel("")
    plt.ylabel("ΔCJ$_{\\mathrm{Rosette}}$ − ΔCJ$_{\\mathrm{Group}}$", fontsize=16)
    plt.axhline(y=0, color='#bcbcbc', linewidth=2,  linestyle='--')
    capitalize_xticks(ax)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    sns.despine()
    plt.legend([], frameon=False)
    plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)

    filters = list(palette.keys())
    for m_index, m in enumerate(combining_method):   
        for f_index, f in enumerate(filters):
            x_pos = m_index + (f_index - 2) * 0.2 + 0.1 -1
            df_sub = df_all[(df_all['Combination'] == m) & (df_all['Filtering'] == f)]
            if not df_sub.empty:
                median = np.median(df_sub["ΔCJ Difference"].values)
                MAD = median_abs_deviation(df_sub["ΔCJ Difference"].values)
                w = 0.07
                plt.plot([x_pos - w, x_pos + w], [median, median], color="#232323", zorder=10)
                plt.plot([x_pos + w, x_pos + w], [median + MAD, median - MAD], color="#bcbcbc", zorder=10)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/diffCJ_vs_rosette_all_filters.png", dpi=300)
    plt.show()
    
    # Save CSVs for further analysis
    data_df.to_csv(f"{output_dir}/diffCJ_stats.csv", index=False)
    differences_df.to_csv(f"{output_dir}/diffCJ_vs_rosette.csv", index=False)

    print("Second Table (Difference Between Rosette and Other Combinations):")
    display_cols = ['Key', 'Filtering', 'Combination', 'ΔCJ Difference']
    print(differences_df[display_cols].sort_values(by=["Key", "Filtering", "Combination"]))

def plot_diffCJ_scatterplot(list_dfs):
    """
    Creates a multi-panel scatterplot showing Mean ΔCJ (≥) vs number of tools (n_tools) 
    across combinations and methods, using filtering strategies as hue.
    
    Parameters:
        list_dfs (list of pd.DataFrame): List of DataFrames, one per method.
    """

    # Add 'n_tools' column
    for df in list_dfs:
        if 'Key' in df.columns:
            df['n_tools'] = df['Key'].apply(lambda x: len(str(x).split('_')))

    list_methods = ["Circle-Seq", "ATAC-seq", "CNT", "RNASE"]
    combinations = ['intersect', 'rosette', 'double', 'union', 'unique']
    n_rows = len(combinations)
    n_cols = len(list_methods)

    palette = {
        'unfilter': '#d46014',
        'filter-split': '#ddcd3d',
        'filter-duplicates': '#064b76ff',
        'filter': '#63bdf6ff'
    }

    fig = plt.figure(figsize=(20, 8))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, wspace=0.2, hspace=0.15)

    all_handles = []
    all_labels = []

    for c, method in enumerate(list_methods):
        df = list_dfs[c]

        for r, combination in enumerate(combinations):
            ax = fig.add_subplot(gs[r, c])
            subset = df[df['Combination'] == combination]

            if not all(col in subset.columns for col in ['Mean ΔCJ (≥)', 'n_tools', 'Filtering']):
                ax.text(0.5, 0.5, "Missing columns", ha='center', va='center', fontsize=12, color='red')
                ax.set_axis_off()
                continue

            sns.scatterplot(
                data=subset,
                x='Mean ΔCJ (≥)',
                y='n_tools',
                hue='Filtering',
                palette=palette,
                s=50,
                ax=ax
            )

            if r == 0:
                ax.set_title(method, fontsize=16)

            if c == 0:
                ax.set_ylabel(combination, fontsize=16, rotation=0, ha='right', va='center')
            else:
                ax.set_ylabel('')

            ax.set_yticks([3, 4, 5])
            ax.set_ylim([2.5, 5.5])
            ax.tick_params(axis='y', labelsize=16)

            if r == 4:
                ax.set_xlabel('Mean ΔCJ (≥)', fontsize=16)
                ax.tick_params(axis='x', labelsize=16)

                if not subset.empty:
                    x_min = subset['Mean ΔCJ (≥)'].min()
                    x_max = subset['Mean ΔCJ (≥)'].max()
                    xticks = np.linspace(x_min, x_max, num=4)
                    ax.set_xticks(xticks)
                    ax.set_xticklabels([f"{x:.2f}" for x in xticks])
            else:
                ax.set_xticks([])
                ax.set_xticklabels([])
                ax.set_xlabel('')

            ax.get_legend().remove()

            if not all_handles:
                all_handles, all_labels = ax.get_legend_handles_labels()

    fig.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("results/eccDNA/real/diffCJ_scatterplot.png", dpi=300, bbox_inches="tight")
    plt.show()


def compute_ratios(df, method_name):
    """
    Compute various ratios based on combination methods for a given dataset.

    Parameters:
    - df (pd.DataFrame): Input DataFrame with combinations and counts
    - method_name (str): Name of the method (used for labeling)

    Returns:
    - pd.DataFrame: DataFrame with computed ratios and method label
    """
    df['Ratio_Union_Rosette'] = float('nan')
    df['Ratio_Union_Intersect'] = float('nan')
    df['Ratio_Rosette_Double'] = float('nan')

    for key in df['Key'].unique():
        for filtering in df['Filtering'].unique():
            subset = df[(df['Key'] == key) & (df['Filtering'] == filtering)]
            try:
                n_union = subset.loc[subset['Combination'] == 'union', 'CJ ≥ threshold (n)'].values[0]
                n_rosette = subset.loc[subset['Combination'] == 'rosette', 'CJ ≥ threshold (n)'].values[0]
                n_intersect = subset.loc[subset['Combination'] == 'intersect', 'CJ ≥ threshold (n)'].values[0]
                n_double = subset.loc[subset['Combination'] == 'double', 'CJ ≥ threshold (n)'].values[0]

                ratio_union_rosette = n_union / n_rosette if n_rosette else float('nan')
                ratio_union_intersect = n_union / n_intersect if n_intersect else float('nan')
                ratio_rosette_double = n_rosette / n_double if n_double else float('nan')

                mask = (df['Key'] == key) & (df['Filtering'] == filtering)
                df.loc[mask, 'Ratio_Union_Rosette'] = ratio_union_rosette
                df.loc[mask, 'Ratio_Union_Intersect'] = ratio_union_intersect
                df.loc[mask, 'Ratio_Rosette_Double'] = ratio_rosette_double

            except IndexError:
                print(f"⚠️ Missing combination for Key={key}, Filtering={filtering} in {method_name}")

    df_cut = df.drop_duplicates(subset=['Key', 'Filtering'])[
        ['Key', 'Filtering', 'Ratio_Union_Rosette', 'Ratio_Union_Intersect', 'Ratio_Rosette_Double']
    ].copy()

    df_cut['Method'] = method_name
    df_cut['Ratio_Rosette_Intersect'] = df_cut['Ratio_Union_Intersect'] / df_cut['Ratio_Union_Rosette']

    return df_cut


def process_and_plot_ratios(paths):
    """
    Process ratio metrics from multiple methods and create violin plots.
    Saves individual plots and combined_ratios.csv in each method's folder.

    Parameters:
    - paths (list of tuples): List of (csv_path, method_name) pairs

    Returns:
    - pd.DataFrame: Combined DataFrame with ratio metrics for all methods
    """
    df_combined = pd.DataFrame()

    # Define color palette
    palette = {
        'unfilter': '#d46014',
        'filter_split': '#ddcd3d',
        'filter_duplicates': '#064b76ff',
        'filter': '#63bdf6ff'
    }

    for path, method in paths:
        df = pd.read_csv(path)
        if 'Key' in df.columns:
            df['n_tools'] = df['Key'].apply(lambda x: len(str(x).split('_')))
        df_cut = compute_ratios(df, method)
        df_combined = pd.concat([df_combined, df_cut], ignore_index=True)

        # Create plot
        plt.figure(figsize=(7, 4))
        ax = plt.gca()
        sns.violinplot(data=df_cut, x='Filtering', y='Ratio_Rosette_Intersect', palette=palette, ax=ax)
        ax.set_ylabel("Rosette / Intersect", fontsize=16)
        ax.set_xlabel('', fontsize=16)
        ax.set_xticklabels(["unfilter", "filter-split", "filter-duplicates", "filter"], fontsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_title('', fontsize=16)
        ax.set_ylim(bottom=0)

        sns.despine()
        plt.tight_layout()

        # Save to same directory as input CSV
        method_dir = os.path.dirname(path)
        os.makedirs(method_dir, exist_ok=True)

        plot_path = os.path.join(method_dir, f"{method}_rosette_intersect_violinplot.png")
        csv_path = os.path.join(method_dir, "combined_ratios.csv")

        plt.savefig(plot_path, dpi=300)
        plt.show()
        plt.close()

        # Save CSV for the individual method
        df_cut.to_csv(csv_path, index=False)

    return df_combined

