import os
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

def upset_plot(matrix_dir, output_dir, tools):
    """
    Creates an UpSet plot with row order matching `tools` and column order by cardinality.

    Parameters:
        matrix_dir (str): Path to the CSV file containing the data.
        output_dir (str): Path where the output plot image will be saved.
        tools (list): List of tool names to define row order (top to bottom).
    """
    # Load data
    df = pd.read_csv(matrix_dir)
    df = df[tools]  # ensure consistent tool columns

    # Build memberships
    memberships = [[tool for tool in tools if row[tool]] for _, row in df.iterrows()]
    data = from_memberships(memberships)

    # Reorder MultiIndex levels to match `tools`
    index_df = data.index.to_frame(index=False)
    data.index = pd.MultiIndex.from_frame(index_df[tools])  # Boolean, no categories!

    # Create UpSet plot
    upset = UpSet(
        data,
        subset_size="count",            # Sort bars by size
        sort_by="cardinality",          # Largest intersections on left
        show_counts=True,
        facecolor="lightsteelblue",
        element_size=45,
        intersection_plot_elements=5
    )

    fig = plt.figure(figsize=(12, 5))
    upset.plot(fig=fig)

    # Font tweaks
    plt.ylabel("Count", fontsize=13)
    plt.yticks(fontsize=13)
    for ax in fig.axes:
        for label in ax.get_xticklabels():
            label.set_fontsize(13)
        for label in ax.get_yticklabels():
            label.set_fontsize(13)
        for text in ax.texts:
            text.set_fontsize(13)

    # Save
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "upset_plot.png"), dpi=300, bbox_inches="tight")
    plt.show()

