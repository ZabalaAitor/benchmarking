import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.metrics import jaccard_score
from upsetplot import UpSet, from_memberships

def upset_plot(matrix_dir, output_dir, tools):
    """
    Creates an UpSet plot from a CSV file and saves it to the specified output path.

    Parameters:
        matrix_dir (str): Path to the CSV file containing the data.
        output_dir (str): Path where the output plot image will be saved.
        tools (list): List of tools names to consider for the UpSet plot.
    """
    # Load the data from CSV
    df = pd.read_csv(matrix_dir)
    memberships = [ [tool for tool in tools if row[tool]] for index, row in df.iterrows()]
    example = from_memberships(memberships)
    upset = UpSet(example, facecolor="lightsteelblue", element_size=35, intersection_plot_elements=5, show_counts=True, sort_by='cardinality', subset_size='count')
    upset.plot()
    plt.ylabel('Count', fontsize=12)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(output_dir + 'upset_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
