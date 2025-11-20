import os
import pandas as pd
import pysam
import pyBigWig
from scipy.stats import binom, kruskal, median_abs_deviation
from statsmodels.stats.multitest import multipletests as multi
from matplotlib import pyplot as plt
import seaborn as sns
import scikit_posthocs as sp
import numpy as np
import math
import matplotlib.gridspec as gridspec

def process_circle_matrix(matrix_file, bam_file, output_file,
                                   mappability_bw=None, map_window=30,
                                   use_chr_prefix=False, N_offset=1,
                                   verbose=False):
    """ 
    Processes a circle matrix and extracts read counts from a BAM file, 
    adjusting CJ imbalance by mappability and read quality. 
    
    Parameters: 

    matrix_file (str): CSV with circle regions 
    bam_file (str): BAM file path 
    output_file (str): CSV output path 
    mappability_bw (str): BigWig file for per-base mappability (optional) 
    map_window (int): Window size to average mappability 
    use_chr_prefix (bool): Add 'chr' prefix to chromosome 
    names N_offset (int): Offset for read counting around CJ 
    verbose (bool): print debug info 
    """

    def safe_bw_values(bw, chrom, start, end, chrom_len):
        bw_chrom_prefix = any(c.startswith("chr") for c in bw.chroms())
        bw_chrom = chrom
        if bw_chrom_prefix and not chrom.startswith('chr'):
            bw_chrom = 'chr' + chrom
        if bw_chrom not in bw.chroms():
            return None
        start = max(0, start)
        end = min(end, chrom_len-1)
        if end <= start:
            return None
        return bw.values(bw_chrom, start, end, numpy=True)

    def binomial_two_sided(k_success, k_total, p_success):
        """Exact two-sided binomial test."""
        if k_total == 0:
            return 1.0
        pmf_obs = binom.pmf(k_success, k_total, p_success)
        pmfs = binom.pmf(range(0, k_total+1), k_total, p_success)
        eps = 1e-12
        return float(pmfs[pmfs <= (pmf_obs + eps)].sum())

    bam = pysam.AlignmentFile(bam_file, "rb")
    bw = None
    chrom_dict = None
    if mappability_bw:
        bw = pyBigWig.open(mappability_bw)
        chrom_dict = bw.chroms()
        if not use_chr_prefix:
            chrom_dict = {k.replace("chr",""): v for k,v in chrom_dict.items()}

    df = pd.read_csv(matrix_file)
    results = []

    for idx, row in df.iterrows():
        circle = row['Circles']
        try:
            chrom, coords = circle.split(":")
            start, end = map(int, coords.split("-"))
            if start > end:
                start, end = end, start
        except ValueError:
            continue

        if use_chr_prefix and not chrom.startswith('chr'):
            chrom = 'chr'+chrom

        try:
            tools_detected = row.iloc[1:].sum()

            # Fetch reads around junctions
            CJ1_reads = list(bam.fetch(chrom, start-N_offset, start+N_offset))
            CJ2_reads = list(bam.fetch(chrom, end-N_offset, end+N_offset))

            CJ1_count = len(CJ1_reads)
            CJ2_count = len(CJ2_reads)
            CJ_reads = CJ1_count + CJ2_count
            diff_CJ = abs(CJ1_count - CJ2_count)

            ratio_CJ1 = CJ1_count / CJ_reads if CJ_reads > 0 else 0
            ratio_CJ2 = CJ2_count / CJ_reads if CJ_reads > 0 else 0
            rel_diff_CJ = abs(ratio_CJ1 - ratio_CJ2)

            CJ_unique_reads = len(set(r.query_name for r in CJ1_reads + CJ2_reads))
            total_reads = len({r.query_name for r in bam.fetch(chrom, start-N_offset, end+N_offset)})
            ratio = CJ_unique_reads / total_reads if total_reads > 0 else 0

            # --- MAPQ values (flat lists) ---
            MAPQ_CJ1 = [r.mapping_quality for r in CJ1_reads]
            MAPQ_CJ2 = [r.mapping_quality for r in CJ2_reads]

            # --- Weighted counts based on MAPQ ---
            weights_CJ1 = [1 - 10**(- (q) / 10) for q in MAPQ_CJ1]
            weights_CJ2 = [1 - 10**(- (q) / 10) for q in MAPQ_CJ2]

            kL_w = sum(weights_CJ1)
            kR_w = sum(weights_CJ2)
            n_w = kL_w + kR_w
            diff_weighted = abs(kL_w - kR_w)
            ratio_wL = kL_w / n_w if n_w > 0 else 0
            ratio_wR = kR_w / n_w if n_w > 0 else 0
            rel_diff_weighted = abs(ratio_wL - ratio_wR)

            # Mean MAPQ for additional tracking
            MAPQ_mean_CJ1 = np.mean(MAPQ_CJ1) if MAPQ_CJ1 else np.nan
            MAPQ_mean_CJ2 = np.mean(MAPQ_CJ2) if MAPQ_CJ2 else np.nan
            MAPQ_diff = abs(MAPQ_mean_CJ1 - MAPQ_mean_CJ2)

            # --- Mappability ---
            p_L = 0.5
            M_L, M_R = None, None
            if bw and chrom in chrom_dict:
                chrom_len = chrom_dict[chrom]
                vals_L = safe_bw_values(bw, chrom, start, start+map_window, chrom_len)
                vals_R = safe_bw_values(bw, chrom, end-map_window, end, chrom_len)
                if vals_L is not None and vals_R is not None:
                    M_L = float(pd.Series(vals_L).mean(skipna=True))
                    M_R = float(pd.Series(vals_R).mean(skipna=True))
                    if (M_L + M_R) > 0:
                        p_L = M_L / (M_L + M_R)

            # --- Binomial tests ---
            if CJ_reads == 0:
                p_raw = np.nan
            else:
                p_raw = binomial_two_sided(CJ1_count, CJ_reads, 0.5)
            kL_obs = int(math.ceil(kL_w))
            kR_obs = int(math.ceil(kR_w))
            k_total = kL_obs + kR_obs
            if k_total == 0:
                p_weighted = np.nan
            else:
                p_weighted = binomial_two_sided(kL_obs, k_total, p_L)

            results.append([
                circle, end-start+1, tools_detected,
                CJ1_count, CJ2_count, CJ_reads, diff_CJ,
                total_reads, ratio,
                ratio_CJ1, ratio_CJ2, rel_diff_CJ,
                kL_w, kR_w, n_w, diff_weighted,
                ratio_wL, ratio_wR, rel_diff_weighted,
                p_raw, p_weighted, 
                MAPQ_CJ1, MAPQ_CJ2,
                MAPQ_mean_CJ1, MAPQ_mean_CJ2, MAPQ_diff,
                M_L, M_R, p_L
            ])
        except ValueError as e:
            if verbose:
                print(f"Skipping invalid chromosome: {chrom} or region {circle}. Error: {e}")
            continue
        
    bam.close()
    if bw:
        bw.close()

    # --- Build DataFrame ---
    columns = [
        "Circle","Length","Tools",
        "Reads_CJ1","Reads_CJ2","Reads_CJ","diff_CJ",
        "Total_Reads","Ratio",
        "Ratio_CJ1","Ratio_CJ2","reldiff_CJ",
        "kL_weighted","kR_weighted","n_weighted","diff_weighted",
        "KL_ratio","KR_ratio","reldiff_weighted",
        "Binom_p_raw","Binom_p_weighted",
        "MAPQ_CJ1","MAPQ_CJ2",
        "MAPQ_mean_CJ1","MAPQ_mean_CJ2","MAPQ_diff",
        "Mappability_L","Mappability_R","p_L"
    ]

    output_df = pd.DataFrame(results, columns=columns)

    # --- FDR corrections ---
    # Initialize FDR columns
    output_df['Binom_p_raw_FDR'] = np.nan
    output_df['Binom_p_weighted_FDR'] = np.nan

    # Compute FDR only on valid p-values
    mask_raw = output_df['Binom_p_raw'].notna()
    if mask_raw.any():
        output_df.loc[mask_raw, 'Binom_p_raw_FDR'] = multi(output_df.loc[mask_raw, 'Binom_p_raw'], method='fdr_bh')[1]

    mask_weighted = output_df['Binom_p_weighted'].notna()
    if mask_weighted.any():
        output_df.loc[mask_weighted, 'Binom_p_weighted_FDR'] = multi(output_df.loc[mask_weighted, 'Binom_p_weighted'], method='fdr_bh')[1]

    # --- Save to CSV ---
    output_df.to_csv(output_file, index=False)


def circle_diff(matrix_dir, matrix_reads, tools, output_dir, group):
    """
    Analyzes circle detection differences across tools based on read counts and statistics.

    Parameters:
        matrix_dir (str): Path to the CSV file containing circle presence matrix.
        matrix_reads (str): Path to the CSV file containing reads and differential CJ statistics.
        tools (list of str): List of tool names to analyze.
        output_dir (str): Directory to save summary CSV and plots.
        group (str): Group identifier used for naming output files.

    Returns:
        None
    """
    # Define file paths
    matrix_path = os.path.join(matrix_dir)
    matrix_with_reads_path = os.path.join(matrix_reads)

    # Load data
    matrix_df = pd.read_csv(matrix_path)
    reads_df = pd.read_csv(matrix_with_reads_path)

    # Merge on circle ID
    merged_df = pd.merge(
        matrix_df,
        reads_df[[
            "Circle","Length","Tools",
            "Reads_CJ1","Reads_CJ2","Reads_CJ","diff_CJ",
            "Total_Reads","Ratio",
            "Ratio_CJ1","Ratio_CJ2","reldiff_CJ",
            "kL_weighted","kR_weighted","n_weighted","diff_weighted",
            "KL_ratio","KR_ratio","reldiff_weighted",
            "Binom_p_raw","Binom_p_weighted",
            "Binom_p_raw_FDR","Binom_p_weighted_FDR",
            "MAPQ_CJ1","MAPQ_CJ2",
            "Mappability_L","Mappability_R","p_L"
        ]],
        left_on='Circles',
        right_on='Circle',
        how='inner'
    )

    n_threshold = 9
    alpha = 0.05

    summary_stats = []  # List to store summary statistics

    for tool in tools:
        tool_data = merged_df[merged_df[tool] == 1]

        # If no circles detected by this tool, produce NaNs/zeros safely
        n_tool = len(tool_data)
        if n_tool == 0:
            summary_stats.append({
                'Tool': tool,
                'Non_Informative (Unweighted)': np.nan,
                'Non_Informative_% (Unweighted)': np.nan,
                'CJ_Reads_<9 (Unweighted)': np.nan,
                'CJ_Reads_<9_% (Unweighted)': np.nan,
                'CJ_Reads_≥9 (Unweighted)': np.nan,
                'CJ_Reads_≥9_% (Unweighted)': np.nan,
                'ACJ_Median_<9 (Unweighted)': np.nan,
                'ACJ_Mean_<9 (Unweighted)': np.nan,
                'ACJ_Median_≥9 (Unweighted)': np.nan,
                'ACJ_Mean_≥9 (Unweighted)': np.nan,
                'p<CJ_alpha_<9_% (Unweighted)': np.nan,
                'pFDR<CJ_alpha_<9_% (Unweighted)': np.nan,
                'p<CJ_alpha_≥9_% (Unweighted)': np.nan,
                'pFDR<CJ_alpha_≥9_% (Unweighted)': np.nan,
                'p<CJ_alpha_<9_% (Unweighted)': np.nan,
                'pFDR<CJ_alpha_<9_% (Unweighted)': np.nan,
                'Non_Informative (Weighted)': np.nan,
                'Non_Informative_% (Weighted)': np.nan,
                'CJ_Reads_<9 (Weighted)': np.nan,
                'CJ_Reads_<9_% (Weighted)': np.nan,
                'CJ_Reads_≥9 (Weighted)': np.nan,
                'CJ_Reads_≥9_% (Weighted)': np.nan,
                'ACJ_Median_<9 (Weighted)': np.nan,
                'ACJ_Mean_<9 (Weighted)': np.nan,
                'ACJ_Median_≥9 (Weighted)': np.nan,
                'ACJ_Mean_≥9 (Weighted)': np.nan,
                'p<CJ_alpha_<9_% (Weighted)': np.nan,
                'pFDR<CJ_alpha_<9_% (Weighted)': np.nan,
                'p<CJ_alpha_≥9_% (Weighted)': np.nan,
                'pFDR<CJ_alpha_≥9_% (Weighted)': np.nan,
                'p<CJ_alpha_<9_% (Weighted)': np.nan,
                'pFDR<CJ_alpha_<9_% (Weighted)': np.nan,
            })
            continue

        filter_data_nan = tool_data[tool_data['Reads_CJ'] == 0]
        filter_data_max = tool_data[tool_data['Reads_CJ'] >= n_threshold]
        filter_data_min = tool_data[(tool_data['Reads_CJ'] > 0) & (tool_data['Reads_CJ'] < n_threshold)]
        filter_data_nan_weighted = tool_data[tool_data['n_weighted'] == 0]
        filter_data_max_weighted = tool_data[tool_data['n_weighted'] >= n_threshold]
        filter_data_min_weighted = tool_data[(tool_data['n_weighted'] > 0) & (tool_data['n_weighted'] < n_threshold)]

        # safe proportion helper
        def safe_prop(num, denom):
            return (num / denom) if denom > 0 else np.nan

        # Unweighted p-value checks use Binom_p_raw and its FDR
        p_raw_ge_alpha = safe_prop(
            len(filter_data_max[filter_data_max['Binom_p_raw'] < alpha]),
            len(filter_data_max)
        )
        p_raw_ge_alpha_fdr = safe_prop(
            len(filter_data_max[filter_data_max['Binom_p_raw_FDR'] < alpha]),
            len(filter_data_max)
        )
        p_raw_lt_alpha = safe_prop(
            len(filter_data_min[filter_data_min['Binom_p_raw'] < alpha]),
            len(filter_data_min)
        )
        p_raw_lt_alpha_fdr = safe_prop(
            len(filter_data_min[filter_data_min['Binom_p_raw_FDR'] < alpha]),
            len(filter_data_min)
        )

        # Weighted p-value checks use Binom_p_weighted and its FDR
        p_w_ge_alpha = safe_prop(
            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted'] < alpha]),
            len(filter_data_max_weighted)
        )
        p_w_ge_alpha_fdr = safe_prop(
            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted_FDR'] < alpha]),
            len(filter_data_max_weighted)
        )
        p_w_lt_alpha = safe_prop(
            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted'] < alpha]),
            len(filter_data_min_weighted)
        )
        p_w_lt_alpha_fdr = safe_prop(
            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted_FDR'] < alpha]),
            len(filter_data_min_weighted)
        )

        summary_stats.append({
            'Tool': tool,

            # Unweighted data
            'Non_Informative (Unweighted)': len(filter_data_nan),
            'Non_Informative_% (Unweighted)': len(filter_data_nan) / len(tool_data),
            'CJ_Reads_<9 (Unweighted)': len(filter_data_min),
            'CJ_Reads_<9_% (Unweighted)': len(filter_data_min) / len(tool_data),
            'CJ_Reads_≥9 (Unweighted)': len(filter_data_max),
            'CJ_Reads_≥9_% (Unweighted)': len(filter_data_max) / len(tool_data),
            'ACJ_Median_<9 (Unweighted)': filter_data_min['reldiff_CJ'].median() if len(filter_data_min)>0 else np.nan,
            'ACJ_Mean_<9 (Unweighted)': filter_data_min['reldiff_CJ'].mean() if len(filter_data_min)>0 else np.nan,
            'ACJ_Median_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].median() if len(filter_data_max)>0 else np.nan,
            'ACJ_Mean_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].mean() if len(filter_data_max)>0 else np.nan,
            'p<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha,
            'pFDR<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha_fdr,
            'p<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha,
            'pFDR<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha_fdr,

            # Weighted data
            'Non_Informative (Weighted)': len(filter_data_nan_weighted),
            'Non_Informative_% (Weighted)': len(filter_data_nan_weighted) / len(tool_data),
            'CJ_Reads_<9 (Weighted)': len(filter_data_min_weighted),
            'CJ_Reads_<9_% (Weighted)': len(filter_data_min_weighted) / len(tool_data),
            'CJ_Reads_≥9 (Weighted)': len(filter_data_max_weighted),
            'CJ_Reads_≥9_% (Weighted)': len(filter_data_max_weighted) / len(tool_data),
            'ACJ_Median_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].median() if len(filter_data_min_weighted)>0 else np.nan,
            'ACJ_Mean_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].mean() if len(filter_data_min_weighted)>0 else np.nan,
            'ACJ_Median_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].median() if len(filter_data_max_weighted)>0 else np.nan,
            'ACJ_Mean_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].mean() if len(filter_data_max_weighted)>0 else np.nan,
            'p<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha,
            'pFDR<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha_fdr,
            'p<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha,
            'pFDR<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha_fdr,
        })

    # Save summary statistics as CSV
    summary_df = pd.DataFrame(summary_stats)

    # Define output path and save
    output_csv_path = os.path.join(output_dir, f'statistics_summary_{group}.csv')
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    summary_df.to_csv(output_csv_path, index=False)

    # Prepare data for plotting
    plot_data = []
    for tool in tools:
        tool_data = merged_df[merged_df[tool] == 1]
        for _, row in tool_data.iterrows():
            plot_data.append({
                'Tool': tool,
                'diff CJ': row.get('reldiff_CJ', np.nan),
                'diff_weighted': row.get('reldiff_weighted', np.nan)
            })

    # Define color palette (keeps your original colors)
    palette = {
        'Simulated': '#b3b3b3',
        'unfilter': '#d46014',
        'filter-split': '#ddcd3d',
        'filter-duplicates': '#064b76ff',
        'filter': '#63bdf6ff'
    }

    # Create DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # If plot_df is empty, skip plotting gracefully
    if plot_df.empty:
        print("No data to plot.")
        return

    # Plot ΔCJ (unweighted)
    plt.figure(figsize=(8.5, 4))
    ax = sns.boxplot(x="Tool", y="diff CJ", data=plot_df, palette=palette)
    plt.xlabel('')
    plt.ylabel('ΔCJ', fontsize=16)
    xtick_labels = tools
    plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.tight_layout()

    # Add counts on top of each box
    # Compute counts for ΔCJw
    group_counts_weighted = plot_df.groupby("Tool")["diff_weighted"].count().to_dict()

    # Add counts on top of each box
    for i, tool in enumerate(tools):
        count = group_counts_weighted.get(tool, 0)
        if count > 0:
            y_text = 1.05
            ax.text(i, y_text, str(count), ha='center', va='bottom', fontsize=12)

    # Save the plot
    save_path = os.path.join(output_dir, f'{group}_diff_cj.png')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    # Plot ΔCJw (weighted)
    plt.figure(figsize=(8.5, 4))
    ax = sns.boxplot(x="Tool", y="diff_weighted", data=plot_df, palette=palette)
    plt.xlabel('')
    plt.ylabel('ΔCJ', fontsize=16)
    plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
    sns.despine()
    plt.tight_layout()

    # Compute counts for ΔCJw
    group_counts_weighted = plot_df.groupby("Tool")["diff_weighted"].count().to_dict()

    # Add counts on top of each box
    for i, tool in enumerate(tools):
        count = group_counts_weighted.get(tool, 0)
        if count > 0:
            y_text = 1.05
            ax.text(i, y_text, str(count), ha='center', va='bottom', fontsize=12)


    # Save the plot
    save_path = os.path.join(output_dir, f'{group}_diff_cj_weighted.png')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    # Group data by tool for Kruskal-Wallis test (use diff_weighted)
    grouped = [group["diff_weighted"].values for _, group in plot_df.groupby("Tool") if not group.empty]

    # Perform Kruskal-Wallis test
    if len(grouped) > 1:
        stat, p_val = kruskal(*grouped)
        print(f"Kruskal-Wallis test across: H = {stat:.4f}, p = {p_val:.4e}")
    else:
        print("Not enough groups to perform Kruskal-Wallis test.")

def plot_mappability_and_mapq(matrix_file, read_file, output_dir):
    """
    Plot violin plots of Mappability and MAPQ for True Positives, False Negatives, and False Positives per tool.

    Parameters:
        matrix_file (str): CSV file with columns
        read_file (str): CSV with circle info
        output_dir (str): Directory to save the plots.
    """

    # === Read data ===
    matrix = pd.read_csv(matrix_file)
    true_df = pd.read_csv(read_file)

    # Normalize column names
    if 'Circles' in matrix.columns:
        matrix.rename(columns={'Circles': 'Circle'}, inplace=True)
    elif 'Circle' not in matrix.columns:
        matrix.rename(columns={matrix.columns[0]: 'Circle'}, inplace=True)

    # Ensure numeric 0/1
    for col in matrix.columns[1:]:
        matrix[col] = pd.to_numeric(matrix[col], errors='coerce').fillna(0).astype(int)

    tools = [c for c in matrix.columns if c not in ['Circle', 'Simulated']]

    # === Colorblind-friendly palette ===
    colorblind_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']

    # === Helper: gather data by condition ===
    def gather_data(condition, base_df, category):
        mapp_data = []
        mapq_data = []
        mapp_min_data = []
        mapp_max_data = []
        mapq_min_data = []
        mapq_max_data = []

        for tool in tools:
            subset = matrix.query(condition(tool))
            merged = subset.merge(base_df, on='Circle', how='inner')

            for _, entry in merged.iterrows():

                # --- ORIGINAL VALUES (KEEP) ---
                mapp_data.append({'Tool': tool, 'Value': entry['Mappability_L'], 'Category': category})
                mapp_data.append({'Tool': tool, 'Value': entry['Mappability_R'], 'Category': category})

                mapq_data.append({'Tool': tool, 'Value': entry['MAPQ_mean_CJ1'], 'Category': category})
                mapq_data.append({'Tool': tool, 'Value': entry['MAPQ_mean_CJ2'], 'Category': category})

                # === NEW VALUES: MIN & MAX ===
                mapp_L = entry['Mappability_L']
                mapp_R = entry['Mappability_R']
                mapq_L = entry['MAPQ_mean_CJ1']
                mapq_R = entry['MAPQ_mean_CJ2']

                mapp_min_data.append({'Tool': tool, 'Value': min(mapp_L, mapp_R), 'Category': category})
                mapp_max_data.append({'Tool': tool, 'Value': max(mapp_L, mapp_R), 'Category': category})

                mapq_min_data.append({'Tool': tool, 'Value': min(mapq_L, mapq_R), 'Category': category})
                mapq_max_data.append({'Tool': tool, 'Value': max(mapq_L, mapq_R), 'Category': category})

        return (
            mapp_data, mapq_data,
            mapp_min_data, mapp_max_data,
            mapq_min_data, mapq_max_data
        )

    # === Define conditions ===
    tp_condition = lambda tool: f"(Simulated == 1 and `{tool}` == 1)"
    fn_condition = lambda tool: f"(Simulated == 1 and `{tool}` == 0)"
    fp_condition = lambda tool: f"(Simulated == 0 and `{tool}` == 1)"

    # === Collect data ===
    (
        tp_mapp, tp_mapq,
        tp_mapp_min, tp_mapp_max,
        tp_mapq_min, tp_mapq_max
    ) = gather_data(tp_condition, true_df, 'True Positive')

    (
        fn_mapp, fn_mapq,
        fn_mapp_min, fn_mapp_max,
        fn_mapq_min, fn_mapq_max
    ) = gather_data(fn_condition, true_df, 'False Negative')

    (
        fp_mapp, fp_mapq,
        fp_mapp_min, fp_mapp_max,
        fp_mapq_min, fp_mapq_max
    ) = gather_data(fp_condition, true_df, 'False Positive')

    # === Combine ORIGINAL ===
    all_mapp = pd.DataFrame(tp_mapp + fn_mapp + fp_mapp)
    all_mapq = pd.DataFrame(tp_mapq + fn_mapq + fp_mapq)

    # === Combine NEW (MIN/MAX) ===
    all_mapp_min = pd.DataFrame(tp_mapp_min + fn_mapp_min + fp_mapp_min)
    all_mapp_max = pd.DataFrame(tp_mapp_max + fn_mapp_max + fp_mapp_max)

    all_mapq_min = pd.DataFrame(tp_mapq_min + fn_mapq_min + fp_mapq_min)
    all_mapq_max = pd.DataFrame(tp_mapq_max + fn_mapq_max + fp_mapq_max)

    # === Plotting ===
    def normalize_tool_name(tool):
        if 'ecc_finder' in tool:
            return tool.replace('-', '\n')
        return tool

    def plot_violin(data, y_label, title_prefix):
        data = data.copy()
        data['Tool'] = data['Tool'].apply(normalize_tool_name)

        # Get all tools present in the full dataset
        all_tools = list(pd.unique(data['Tool']))

        for category in ['True Positive', 'False Negative', 'False Positive']:
            subset = data[data['Category'] == category]

            plt.figure(figsize=(10, 4))
            ax = sns.violinplot(
                data=subset,
                x='Tool',
                y='Value',
                order=all_tools,                     
                inner='box',
                palette=colorblind_palette,
                cut=0
            )

            # Remove top and right spines
            sns.despine(top=True, right=True)

            # Count values per tool (missing tools get 0)
            tool_counts = subset['Tool'].value_counts().to_dict()

            # Determine y-range for label placement
            ymax = subset['Value'].max() if not subset.empty else 1

            # Add counts (including 0) above each violin
            for i, tool in enumerate(all_tools):
                count = tool_counts.get(tool, 0)
                ax.text(
                    i, ymax + (0.05 * ymax),
                    str(count),
                    ha='center', va='bottom',
                    fontsize=14
                )

            plt.ylabel(y_label, fontsize=16)
            plt.xlabel("")
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.tight_layout()
            ax.set_axisbelow(True)
            plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
            plt.savefig(f"{output_dir}/{title_prefix}_{category.replace(' ', '_')}.png", dpi=300)
            plt.show()


    # === ORIGINAL PLOTS ===
    plot_violin(all_mapp, "Mappability", "Distribution_of_Mappability_per_Tool")
    plot_violin(all_mapq, "MAPQ Mean", "Distribution_of_MAPQ_Mean_per_Tool")

    plot_violin(all_mapp_min, "Min Mappability", "Distribution_of_Mappability_Min_per_Tool")
    plot_violin(all_mapp_max, "Max Mappability", "Distribution_of_Mappability_Max_per_Tool")

    plot_violin(all_mapq_min, "Min MAPQ Mean", "Distribution_of_MAPQ_Min_per_Tool")
    plot_violin(all_mapq_max, "Max MAPQ Mean", "Distribution_of_MAPQ_Max_per_Tool")

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

def circle_diff_real(matrix_dir, tools, filtering_methods, data, n_threshold=9, alpha=0.05):
    """
    Perform comparative analysis of detected circular DNA/RNA across tools and filtering strategies.

    This function loads detection matrices and associated CJ Read counts for various combinations of 
    tools and filtering methods, computes summary statistics, categorizes circles by CJ Read threshold, 
    generates comparative boxplots, and performs non-parametric statistical tests (Kruskal-Wallis and Dunn's test) 
    to evaluate differences in ΔCJ (difference in circular junction reads) distributions.

    Summary statistics include the number and percentage of high-confidence circles (above threshold), 
    median/mean ΔCJ values, and the proportion of statistically significant changes.

    Parameters:
        matrix_dir (str): Path to the base directory containing 'matrix.csv' and 'matrix_with_reads.csv' files 
                        for each filtering method and dataset.
        tools (list): List of tool names to include in the analysis (e.g., ['Circle-Map', 'CIRCexplorer2']).
        filtering_methods (list): List of filtering method names (e.g., ['filter', 'unfilter']).
        data (str): Name of the dataset (used in file paths and output plot titles).
        n_threshold (int, optional): Threshold for CJ Read count to define high-confidence circles (default is 9).
        alpha (float, optional): Significance level for p-value thresholds (default is 0.05).

    Returns:
        None

    Outputs:
        CSV file with summary statistics for each tool and filtering method.
        PNG boxplot visualizing ΔCJ distributions across tools and filters.
        Excel file containing Kruskal-Wallis and Dunn's test results for ΔCJ comparisons.
    """

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
            matrix_with_reads[[
                "Circle","Length","Tools",
                "Reads_CJ1","Reads_CJ2","Reads_CJ","diff_CJ",
                "Total_Reads","Ratio",
                "Ratio_CJ1","Ratio_CJ2","reldiff_CJ",
                "kL_weighted","kR_weighted","n_weighted","diff_weighted",
                "KL_ratio","KR_ratio","reldiff_weighted",
                "Binom_p_raw","Binom_p_weighted",
                "Binom_p_raw_FDR","Binom_p_weighted_FDR",
                "MAPQ_CJ1","MAPQ_CJ2",
                "Mappability_L","Mappability_R","p_L"]],
            left_on='Circles',
            right_on='Circle',
            how='inner'
        )

        tool_data = merged_data[merged_data[tool] == 1]
        filter_data_nan = tool_data[tool_data['Reads_CJ'] == 0]
        filter_data_max = tool_data[tool_data['Reads_CJ'] >= n_threshold]
        filter_data_min = tool_data[(tool_data['Reads_CJ'] > 0) & (tool_data['Reads_CJ'] < n_threshold)]
        filter_data_nan_weighted = tool_data[tool_data['n_weighted'] == 0]
        filter_data_max_weighted = tool_data[tool_data['n_weighted'] >= n_threshold]
        filter_data_min_weighted = tool_data[(tool_data['n_weighted'] > 0) & (tool_data['n_weighted'] < n_threshold)]

        # safe proportion helper
        def safe_prop(num, denom):
            return (num / denom) if denom > 0 else np.nan

        # Unweighted p-value checks use Binom_p_raw and its FDR
        p_raw_ge_alpha = safe_prop(
            len(filter_data_max[filter_data_max['Binom_p_raw'] < alpha]),
            len(filter_data_max)
        )
        p_raw_ge_alpha_fdr = safe_prop(
            len(filter_data_max[filter_data_max['Binom_p_raw_FDR'] < alpha]),
            len(filter_data_max)
        )
        p_raw_lt_alpha = safe_prop(
            len(filter_data_min[filter_data_min['Binom_p_raw'] < alpha]),
            len(filter_data_min)
        )
        p_raw_lt_alpha_fdr = safe_prop(
            len(filter_data_min[filter_data_min['Binom_p_raw_FDR'] < alpha]),
            len(filter_data_min)
        )

        # Weighted p-value checks use Binom_p_weighted and its FDR
        p_w_ge_alpha = safe_prop(
            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted'] < alpha]),
            len(filter_data_max_weighted)
        )
        p_w_ge_alpha_fdr = safe_prop(
            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted_FDR'] < alpha]),
            len(filter_data_max_weighted)
        )
        p_w_lt_alpha = safe_prop(
            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted'] < alpha]),
            len(filter_data_min_weighted)
        )
        p_w_lt_alpha_fdr = safe_prop(
            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted_FDR'] < alpha]),
            len(filter_data_min_weighted)
        )

        summary_stats.append({
            'Tool': tool,
            'Filtering': filtering_method,

            # Unweighted data
            'Non_Informative (Unweighted)': len(filter_data_nan),
            'Non_Informative_% (Unweighted)': len(filter_data_nan) / len(tool_data),
            'CJ_Reads_<9 (Unweighted)': len(filter_data_min),
            'CJ_Reads_<9_% (Unweighted)': len(filter_data_min) / len(tool_data),
            'CJ_Reads_≥9 (Unweighted)': len(filter_data_max),
            'CJ_Reads_≥9_% (Unweighted)': len(filter_data_max) / len(tool_data),
            'ACJ_Median_<9 (Unweighted)': filter_data_min['reldiff_CJ'].median() if len(filter_data_min)>0 else np.nan,
            'ACJ_Mean_<9 (Unweighted)': filter_data_min['reldiff_CJ'].mean() if len(filter_data_min)>0 else np.nan,
            'ACJ_Median_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].median() if len(filter_data_max)>0 else np.nan,
            'ACJ_Mean_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].mean() if len(filter_data_max)>0 else np.nan,
            'p<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha,
            'pFDR<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha_fdr,
            'p<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha,
            'pFDR<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha_fdr,

            # Weighted data
            'Non_Informative (Weighted)': len(filter_data_nan_weighted),
            'Non_Informative_% (Weighted)': len(filter_data_nan_weighted) / len(tool_data),
            'CJ_Reads_<9 (Weighted)': len(filter_data_min_weighted),
            'CJ_Reads_<9_% (Weighted)': len(filter_data_min_weighted) / len(tool_data),
            'CJ_Reads_≥9 (Weighted)': len(filter_data_max_weighted),
            'CJ_Reads_≥9_% (Weighted)': len(filter_data_max_weighted) / len(tool_data),
            'ACJ_Median_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].median() if len(filter_data_min_weighted)>0 else np.nan,
            'ACJ_Mean_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].mean() if len(filter_data_min_weighted)>0 else np.nan,
            'ACJ_Median_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].median() if len(filter_data_max_weighted)>0 else np.nan,
            'ACJ_Mean_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].mean() if len(filter_data_max_weighted)>0 else np.nan,
            'p<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha,
            'pFDR<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha_fdr,
            'p<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha,
            'pFDR<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha_fdr,
        })

        for _, row in filter_data_max.iterrows():
            plot_data.append({
                'Tool': tool,
                'Filtering': filtering_method,
                'Circle': row['Circles'],
                'diff CJ': row['reldiff_CJ'],
                'diff_weighted': np.nan  # placeholder, or merge later
            })

        for _, row in filter_data_max_weighted.iterrows():
            plot_data.append({
                'Tool': tool,
                'Filtering': filtering_method,
                'Circle': row['Circles'],
                'diff CJ': np.nan,  # placeholder
                'diff_weighted': row['reldiff_weighted']
            })

        # Merge unweighted and weighted by 'Circles'
        df_unw = filter_data_max[['Circles', 'reldiff_CJ']].rename(columns={'reldiff_CJ': 'diff CJ'})
        df_w = filter_data_max_weighted[['Circles', 'reldiff_weighted']].rename(columns={'reldiff_weighted': 'diff_weighted'})

        merged_plot_data = pd.merge(df_unw, df_w, on='Circles', how='outer')
        merged_plot_data['Tool'] = tool
        merged_plot_data['Filtering'] = filtering_method

        plot_data.extend(merged_plot_data.to_dict(orient='records'))

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

        # Prepare tool tick labels
        tools_tick = [tool.replace('finder-', 'finder\n') for tool in tools]
        # --- Unweighted boxplot ---
        plot_df_unw = plot_df.dropna(subset=['diff CJ'])
        plt.figure(figsize=(10, 4))
        ax = sns.boxplot(x="Tool", y="diff CJ", hue="Filtering", data=plot_df_unw, palette=palette)
        plt.xlabel('')
        plt.ylabel('ΔCJ', fontsize=16)
        plt.xticks(ticks=range(len(tools)), labels=tools_tick, fontsize=16)
        plt.yticks(fontsize=16)
        ax.legend_.remove()
        plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
        sns.despine()
        plt.tight_layout()

        # annotate counts
        group_counts_unw = plot_df_unw.groupby(['Tool', 'Filtering'])['Circle'].nunique().reset_index(name='count')
        for _, row in group_counts_unw.iterrows():
            tool = row['Tool']
            filtering = row['Filtering']
            count = row['count']
            tool_idx = tools.index(tool)
            filtering_idx = list(palette.keys()).index(filtering)
            total_hue = len(palette)
            spacing_factor = 0.85
            x_position = tool_idx + ((filtering_idx - total_hue / 2) * spacing_factor / total_hue) + 0.1
            ax.text(x_position, 1.02, f'{count}', ha='center', va='bottom', fontsize=14, rotation=45)

        plt.savefig(os.path.join(matrix_dir, data, 'diffCJ.png'), dpi=300, bbox_inches="tight")
        plt.show()

        for filtering_method in filtering_methods:
            df_sub = plot_df.dropna(subset=['diff CJ'])
            df_sub = df_sub[df_sub['Filtering'] == filtering_method]
            if not df_sub.empty:
                groups = [df_sub[df_sub['Tool'] == tool]['diff CJ'].values for tool in tools]

                # Remove empty or NaN-only groups
                groups_clean = [g[~np.isnan(g)] for g in groups if len(g[~np.isnan(g)]) > 0]
                # Keep only groups with >1 unique value
                valid_groups = [g for g in groups_clean if len(np.unique(g)) > 1]

                if len(valid_groups) > 1:
                    stat, p_val = kruskal(*valid_groups)
                    print(f"Kruskal-Wallis test for filtering method {filtering_method}: H = {stat:.4f}, p = {p_val:.4f}")

        for filtering_method in filtering_methods:
            df_sub = plot_df[plot_df['Filtering'] == filtering_method]
            if not df_sub.empty:
                groups = [df_sub[df_sub['Tool'] == tool]['diff CJ'].values for tool in tools]
                if len(groups) > 1:
                    dunn_result = sp.posthoc_dunn(df_sub, val_col='diff CJ', group_col='Tool', p_adjust='bonferroni')
                    dunn_results.append({
                        'Filtering Method': filtering_method,
                        'Test': 'Dunn\'s test (Tools)',
                        'results': dunn_result
                    })
                    print(f"  Dunn's test for {filtering_method}: \n{dunn_result}")



        output_dir = os.path.join(matrix_dir, data)
        os.makedirs(output_dir, exist_ok=True)

        output_excel_path = os.path.join(output_dir, 'statistical_results_weighted.xlsx')
        with pd.ExcelWriter(output_excel_path) as writer:
            # Save Kruskal-Wallis results
            pd.DataFrame(kruskal_results).to_excel(writer, sheet_name='Kruskal-Wallis', index=False)
            
            # Save Dunn's test results
            for result in dunn_results:
                sheet_name = f"Dunn_{result['Filtering Method']}"
                result['results'].to_excel(writer, sheet_name=sheet_name[:31], index=True)  # Excel sheet name limit = 31 chars

        # --- Weighted boxplot ---
        plot_df_w = plot_df.dropna(subset=['diff_weighted'])
        plt.figure(figsize=(10, 4))
        ax = sns.boxplot(x="Tool", y="diff_weighted", hue="Filtering", data=plot_df_w, palette=palette)
        plt.xlabel('')
        plt.ylabel('ΔCJ', fontsize=16)
        plt.xticks(ticks=range(len(tools)), labels=tools_tick, fontsize=16)
        plt.yticks(fontsize=16)
        ax.legend_.remove()
        plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)
        sns.despine()
        plt.tight_layout()

        # annotate counts (weighted)
        group_counts_w = plot_df_w.groupby(['Tool', 'Filtering'])['Circle'].nunique().reset_index(name='count')
        for _, row in group_counts_w.iterrows():
            tool = row['Tool']
            filtering = row['Filtering']
            count = row['count']
            tool_idx = tools.index(tool)
            filtering_idx = list(palette.keys()).index(filtering)
            total_hue = len(palette)
            spacing_factor = 0.85
            x_position = tool_idx + ((filtering_idx - total_hue / 2) * spacing_factor / total_hue) + 0.1
            ax.text(x_position, 1.02, f'{count}', ha='center', va='bottom', fontsize=14, rotation=45)


        plt.savefig(os.path.join(matrix_dir, data, 'diff_weighted.png'), dpi=300, bbox_inches="tight")
        plt.show()

        for filtering_method in filtering_methods:
            df_sub = plot_df.dropna(subset=['diff_weighted'])
            df_sub = df_sub[df_sub['Filtering'] == filtering_method]
            if not df_sub.empty:
                groups = [df_sub[df_sub['Tool'] == tool]['diff_weighted'].values for tool in tools]

                # Remove empty or NaN-only groups
                groups_clean = [g[~np.isnan(g)] for g in groups if len(g[~np.isnan(g)]) > 0]
                # Keep only groups with >1 unique value
                valid_groups = [g for g in groups_clean if len(np.unique(g)) > 1]

                if len(valid_groups) > 1:
                    stat, p_val = kruskal(*valid_groups)
                    kruskal_results.append({
                        'Filtering Method': filtering_method,
                        'Test': 'Kruskal-Wallis (Tools, Weighted)',
                        'H-statistic': stat,
                        'p-value': p_val
                    })
        
        for filtering_method in filtering_methods:
            df_sub = plot_df[plot_df['Filtering'] == filtering_method]
            if not df_sub.empty:
                groups = [df_sub[df_sub['Tool'] == tool]['diff_weighted'].values for tool in tools]
                if len(groups) > 1:
                    dunn_result = sp.posthoc_dunn(df_sub, val_col='diff_weighted', group_col='Tool', p_adjust='bonferroni')
                    dunn_results.append({
                        'Filtering Method': filtering_method,
                        'Test': 'Dunn\'s test (Tools)',
                        'results': dunn_result
                    })
                    print(f"  Dunn's test for {filtering_method}: \n{dunn_result}")


        output_dir = os.path.join(matrix_dir, data)
        os.makedirs(output_dir, exist_ok=True)

        output_excel_path = os.path.join(output_dir, 'statistical_results_weighted.xlsx')
        with pd.ExcelWriter(output_excel_path) as writer:
            # Save Kruskal-Wallis results
            pd.DataFrame(kruskal_results).to_excel(writer, sheet_name='Kruskal-Wallis', index=False)
            
            # Save Dunn's test results
            for result in dunn_results:
                sheet_name = f"Dunn_{result['Filtering Method']}"
                result['results'].to_excel(writer, sheet_name=sheet_name[:31], index=True)  # Excel sheet name limit = 31 chars

    else:
        print("Plotting skipped due to empty DataFrame.")      


def diff_cj_combinations(circle, method, ros_combinations_mapping, filtering_methods, combining_method):
    """
    Analyze ΔCJ (difference in circular junction reads) across combinations and compare against 'rosette'.

    Parameters:
        circle: Type of circle (e.g., 'circDNA')
        method: Method name (e.g., 'ciriquant')
        ros_combinations_mapping: Mapping of sample keys to group metadata
        filtering_methods: List of filtering methods (e.g., ['bwa', 'minimap2'])
        combining_method: List of combination strategies (e.g., ['union', 'rosette'])
    """
    data = []
    differences_data = []
    rosette_diff_cj = {}
    rosette_diff_cj_weighted = {}
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
                    if 'reldiff_CJ' in df.columns and 'Reads_CJ' in df.columns:
                        tool_data = df
                        filter_data_nan = tool_data[tool_data['Reads_CJ'] == 0]
                        filter_data_max = tool_data[tool_data['Reads_CJ'] >= n_threshold]
                        filter_data_min = tool_data[(tool_data['Reads_CJ'] > 0) & (tool_data['Reads_CJ'] < n_threshold)]
                        filter_data_nan_weighted = tool_data[tool_data['n_weighted'] == 0]
                        filter_data_max_weighted = tool_data[tool_data['n_weighted'] >= n_threshold]
                        filter_data_min_weighted = tool_data[(tool_data['n_weighted'] > 0) & (tool_data['n_weighted'] < n_threshold)]

                        if df.empty:
                            continue

                        diff_CJ = filter_data_max['reldiff_CJ'].mean()
                        diff_CJ_weighted = filter_data_max_weighted['reldiff_weighted'].mean()

                        if combination == "rosette":
                            rosette_diff_cj[(file_key, filter_method)] = diff_CJ
                            rosette_diff_cj_weighted[(file_key, filter_method)] = diff_CJ_weighted

                        percent_max = len(filter_data_max) / len(tool_data) if len(tool_data) > 0 else np.nan

                        # safe proportion helper
                        def safe_prop(num, denom):
                            return (num / denom) if denom > 0 else np.nan

                        # Unweighted p-value checks use Binom_p_raw and its FDR
                        p_raw_ge_alpha = safe_prop(
                            len(filter_data_max[filter_data_max['Binom_p_raw'] < alpha]),
                            len(filter_data_max)
                        )
                        p_raw_ge_alpha_fdr = safe_prop(
                            len(filter_data_max[filter_data_max['Binom_p_raw_FDR'] < alpha]),
                            len(filter_data_max)
                        )
                        p_raw_lt_alpha = safe_prop(
                            len(filter_data_min[filter_data_min['Binom_p_raw'] < alpha]),
                            len(filter_data_min)
                        )
                        p_raw_lt_alpha_fdr = safe_prop(
                            len(filter_data_min[filter_data_min['Binom_p_raw_FDR'] < alpha]),
                            len(filter_data_min)
                        )

                        # Weighted p-value checks use Binom_p_weighted and its FDR
                        p_w_ge_alpha = safe_prop(
                            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted'] < alpha]),
                            len(filter_data_max_weighted)
                        )
                        p_w_ge_alpha_fdr = safe_prop(
                            len(filter_data_max_weighted[filter_data_max_weighted['Binom_p_weighted_FDR'] < alpha]),
                            len(filter_data_max_weighted)
                        )
                        p_w_lt_alpha = safe_prop(
                            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted'] < alpha]),
                            len(filter_data_min_weighted)
                        )
                        p_w_lt_alpha_fdr = safe_prop(
                            len(filter_data_min_weighted[filter_data_min_weighted['Binom_p_weighted_FDR'] < alpha]),
                            len(filter_data_min_weighted)
                        )

                        data.append({
                            "Combination": combination,
                            'Key': file_key,
                            'Filtering': filter_method,
                            # Unweighted data
                            'Non_Informative (Unweighted)': len(filter_data_nan),
                            'Non_Informative_% (Unweighted)': len(filter_data_nan) / len(tool_data),
                            'CJ_Reads_<9 (Unweighted)': len(filter_data_min),
                            'CJ_Reads_<9_% (Unweighted)': len(filter_data_min) / len(tool_data),
                            'CJ_Reads_≥9 (Unweighted)': len(filter_data_max),
                            'CJ_Reads_≥9_% (Unweighted)': len(filter_data_max) / len(tool_data),
                            'ACJ_Median_<9 (Unweighted)': filter_data_min['reldiff_CJ'].median() if len(filter_data_min)>0 else np.nan,
                            'ACJ_Mean_<9 (Unweighted)': filter_data_min['reldiff_CJ'].mean() if len(filter_data_min)>0 else np.nan,
                            'ACJ_Median_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].median() if len(filter_data_max)>0 else np.nan,
                            'ACJ_Mean_≥9 (Unweighted)': filter_data_max['reldiff_CJ'].mean() if len(filter_data_max)>0 else np.nan,
                            'p<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha,
                            'pFDR<CJ_alpha_<9_% (Unweighted)': p_raw_lt_alpha_fdr,
                            'p<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha,
                            'pFDR<CJ_alpha_≥9_% (Unweighted)': p_raw_ge_alpha_fdr,
                            

                            # Weighted data
                            'Non_Informative (Weighted)': len(filter_data_nan_weighted),
                            'Non_Informative_% (Weighted)': len(filter_data_nan_weighted) / len(tool_data),
                            'CJ_Reads_<9 (Weighted)': len(filter_data_min_weighted),
                            'CJ_Reads_<9_% (Weighted)': len(filter_data_min_weighted) / len(tool_data),
                            'CJ_Reads_≥9 (Weighted)': len(filter_data_max_weighted),
                            'CJ_Reads_≥9_% (Weighted)': len(filter_data_max_weighted) / len(tool_data),
                            'ACJ_Median_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].median() if len(filter_data_min_weighted)>0 else np.nan,
                            'ACJ_Mean_<9 (Weighted)': filter_data_min_weighted['reldiff_weighted'].mean() if len(filter_data_min_weighted)>0 else np.nan,
                            'ACJ_Median_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].median() if len(filter_data_max_weighted)>0 else np.nan,
                            'ACJ_Mean_≥9 (Weighted)': filter_data_max_weighted['reldiff_weighted'].mean() if len(filter_data_max_weighted)>0 else np.nan,
                            'p<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha,
                            'pFDR<CJ_alpha_<9_% (Weighted)': p_w_lt_alpha_fdr,
                            'p<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha,
                            'pFDR<CJ_alpha_≥9_% (Weighted)': p_w_ge_alpha_fdr,
                        })

                        if combination != "rosette" and (file_key, filter_method) in rosette_diff_cj:
                            rosette_mean = rosette_diff_cj[(file_key, filter_method)]
                            rosette_mean_weighted = rosette_diff_cj_weighted.get((file_key, filter_method), np.nan)

                            diff_cj_difference = rosette_mean - diff_CJ
                            diff_cj_difference_weighted = rosette_mean_weighted - diff_CJ_weighted

                            differences_data.append({
                                'Key': file_key,
                                'Combination': combination,
                                'Filtering': filter_method,
                                'ΔCJ Difference': diff_cj_difference,
                                'ΔCJ Difference Weighted': diff_cj_difference_weighted
                            })

                            
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
                    
    # Convert data lists to DataFrames
    data_df = pd.DataFrame(data)
    print(differences_data[:5])

    differences_df = pd.DataFrame(differences_data)

    # === PLOT 1: Overall Diff CJ ===
    plt.figure(figsize=(7, 4))
    ax = sns.stripplot(data=data_df, x="Combination", y="ACJ_Median_≥9 (Unweighted)", hue="Filtering", jitter=True, dodge=True, palette=palette)
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
                median = np.median(df_sub["ACJ_Median_≥9 (Unweighted)"].values)
                MAD = median_abs_deviation(df_sub["ACJ_Median_≥9 (Unweighted)"].values)
                w = 0.07
                plt.plot([x_pos - w, x_pos + w], [median, median], color="#232323", zorder=10)
                plt.plot([x_pos + w, x_pos + w], [median + MAD, median - MAD], color="#bcbcbc", zorder=10)


    plt.tight_layout()
    plt.savefig(f"{output_dir}/diffCJ_overall.png", dpi=300)
    plt.show()

    # == PLOT 2: Overall Diff CJ ===
    plt.figure(figsize=(7, 4))
    ax = sns.stripplot(data=data_df, x="Combination", y="ACJ_Median_≥9 (Weighted)", hue="Filtering", jitter=True, dodge=True, palette=palette)
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
                median = np.median(df_sub["ACJ_Median_≥9 (Weighted)"].values)
                MAD = median_abs_deviation(df_sub["ACJ_Median_≥9 (Weighted)"].values)
                w = 0.07
                plt.plot([x_pos - w, x_pos + w], [median, median], color="#232323", zorder=10)
                plt.plot([x_pos + w, x_pos + w], [median + MAD, median - MAD], color="#bcbcbc", zorder=10)


    plt.tight_layout()
    plt.savefig(f"{output_dir}/diffCJ_overall_weighted.png", dpi=300)
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
            x_pos = m_index + (f_index - len(filters) / 2) * 0.2 - 1 + 0.1
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
    
    # === PLOT 2: Difference from Rosette ===
    plt.figure(figsize=(7, 4))
    df_all = differences_df[differences_df['Combination'].isin(['union', 'intersect', 'double', 'unique'])]
    ax = sns.stripplot(data=df_all, x="Combination", y="ΔCJ Difference Weighted", hue="Filtering", jitter=True, dodge=True, palette=palette)
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
            x_pos = m_index + (f_index - len(filters) / 2) * 0.2 - 1 + 0.1
            df_sub = df_all[(df_all['Combination'] == m) & (df_all['Filtering'] == f)]
            if not df_sub.empty:
                median = np.median(df_sub["ΔCJ Difference Weighted"].values)
                MAD = median_abs_deviation(df_sub["ΔCJ Difference Weighted"].values)
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
    display_cols = ['Key', 'Filtering', 'Combination', 'ΔCJ Difference', 'ΔCJ Difference Weighted']
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

            if not all(col in subset.columns for col in ['ACJ_Mean_≥9 (Weighted)', 'n_tools', 'Filtering']):
                ax.text(0.5, 0.5, "Missing columns", ha='center', va='center', fontsize=12, color='red')
                ax.set_axis_off()
                continue

            sns.scatterplot(
                data=subset,
                x='ACJ_Mean_≥9 (Weighted)',
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
                    x_min = subset['ACJ_Mean_≥9 (Weighted)'].min()
                    x_max = subset['ACJ_Mean_≥9 (Weighted)'].max()
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
        df (pd.DataFrame): Input DataFrame with combinations and counts
        method_name (str): Name of the method (used for labeling)

    Returns:
        pd.DataFrame: DataFrame with computed ratios and method label
    """
    df['Ratio_Union_Rosette'] = float('nan')
    df['Ratio_Union_Intersect'] = float('nan')
    df['Ratio_Rosette_Double'] = float('nan')

    for key in df['Key'].unique():
        for filtering in df['Filtering'].unique():
            subset = df[(df['Key'] == key) & (df['Filtering'] == filtering)]
            try:
                n_union = subset.loc[subset['Combination'] == 'union', 'CJ_Reads_≥9 (Weighted)'].values[0]
                n_rosette = subset.loc[subset['Combination'] == 'rosette', 'CJ_Reads_≥9 (Weighted)'].values[0]
                n_intersect = subset.loc[subset['Combination'] == 'intersect', 'CJ_Reads_≥9 (Weighted)'].values[0]
                n_double = subset.loc[subset['Combination'] == 'double', 'CJ_Reads_≥9 (Weighted)'].values[0]

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
        paths (list of tuples): List of (csv_path, method_name) pairs

    Returns:
        pd.DataFrame: Combined DataFrame with ratio metrics for all methods
    """
    df_combined = pd.DataFrame()

    # Define color palette
    palette = {
        'unfilter': '#d46014',
        'filter-split': '#ddcd3d',
        'filter-duplicates': '#064b76ff',
        'filter': '#63bdf6ff'
    }

    violin_order = ["unfilter", "filter-split", "filter-duplicates", "filter"]

    for path, method in paths:
        df = pd.read_csv(path)
        if 'Key' in df.columns:
            df['n_tools'] = df['Key'].apply(lambda x: len(str(x).split('_')))
        df_cut = compute_ratios(df, method)
        df_combined = pd.concat([df_combined, df_cut], ignore_index=True)

        # Create plot
        plt.figure(figsize=(7, 4))
        ax = plt.gca()

        sns.violinplot(
            data=df_cut, 
            x='Filtering', 
            y='Ratio_Rosette_Intersect', 
            palette=palette, 
            order=violin_order,
            ax=ax,
            inner='box',
            cut=0
        )

        # Put grid behind the violins
        ax.set_axisbelow(True)
        plt.grid(axis="y", linestyle='--', linewidth=0.5, alpha=0.7)

        ax.set_ylabel("Rosette / Intersect", fontsize=16)
        ax.set_xlabel('', fontsize=16)
        ax.set_xticklabels(violin_order, fontsize=16)
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
