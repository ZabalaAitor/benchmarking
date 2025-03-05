import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

def plot_jaccard_index_comparisons(matrix_files, output_path, palette="colorblind", circle_type='eccDNA'):
    """
    Generate a scatter plot of Jaccard Index comparisons for different coverages and tool pairs.

    Args:
        matrix_files (list): List of paths to CSV files with Jaccard matrices for various coverages.
        output_path (str): Directory to save the generated plot.
        palette (str): Seaborn palette to use for coverages (default is 'colorblind').
        circle_type (str): Circle type.
    """
    def calculate_jaccard(series1, series2):
        """Calculate Jaccard index between two binary series."""
        intersection = (series1 & series2).sum()
        union = (series1 | series2).sum()
        return intersection / union if union > 0 else 0

    # Load matrices and normalize circle IDs
    all_matrices = {}
    coverages = []

    for file in matrix_files:
        coverage = int(os.path.basename(file).split('_')[1][3:-4])
        coverages.append(coverage)
        df = pd.read_csv(file).set_index('Circle')  # Assume 'Circle' column exists
        all_matrices[coverage] = df

    # Extract tool names from the first file
    tool_names = all_matrices[coverages[0]].columns

    # Compute Jaccard index for each pair of tools and each coverage
    jaccard_data = []
    for coverage in coverages:
        df = all_matrices[coverage]
        for tool1, tool2 in itertools.combinations(tool_names, 2):
            jaccard_values = [
                calculate_jaccard(df[tool1] > 0, df[tool2] > 0)
                for circle in df.index
            ]
            jaccard_data.append({
                "Coverage": coverage,
                "Tool Comparison": f"{tool1} vs {tool2}",
                "Jaccard Index": sum(jaccard_values) / len(jaccard_values)  # Average Jaccard index
            })

    # Convert to DataFrame
    jaccard_df = pd.DataFrame(jaccard_data)

    # Assign colors to each coverage using the specified palette
    sns_palette = sns.color_palette(palette, len(coverages))
    coverage_colors = {
        coverage: sns_palette[i]
        for i, coverage in enumerate(sorted(coverages))
    }

    # Plot
    plt.figure(figsize=(6, 4))
    for comparison in jaccard_df["Tool Comparison"].unique():
        for coverage in coverages:
            subset = jaccard_df[(jaccard_df["Tool Comparison"] == comparison) & (jaccard_df["Coverage"] == coverage)]
            plt.scatter(
                subset["Jaccard Index"], [comparison] * len(subset),
                color=coverage_colors[coverage],
                label=f"Coverage {coverage}" if comparison == jaccard_df["Tool Comparison"].unique()[0] else None,
                alpha=0.7
            )

    # Style updates
    plt.title(f"Pairwise Comparison of Tools for Jaccard Index\nfor {circle_type} Detection Software Across Coverages", fontsize=14)
    plt.xlabel("Jaccard Index", fontsize=12)
    plt.ylabel("", fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    sns.despine()

    # Add a single legend with numeric coverage values
    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', markersize=8) 
            for coverage, color in coverage_colors.items()]
    labels = list(coverage_colors.keys())
    plt.legend(handles, labels, title="Coverage", frameon=False, loc='upper left', 
            bbox_to_anchor=(1.05, 0.65), fontsize=10, ncol=1)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Make room for legend

    # Save plot
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(os.path.join(output_path, "jaccard_index_comparisons.png"), dpi=300)
    plt.show()