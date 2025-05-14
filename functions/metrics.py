import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_circular_detection(csv_path, title, save_directory=None, circle_type='eccDNA'):
    """
    Generates a line plot for the given coverage data from a CSV file.

    Parameters:
        csv_path (str): Path to the CSV file containing coverage data.
        title (str): Title for the plot.
        save_directory (str, optional): Directory where the plot should be saved. If not provided, the plot is not saved.
        circle_type (str): Circle type.
    """
    # Import CSV file
    df = pd.read_csv(csv_path)

    # Reshape data for plotting
    df_melted = df.melt(id_vars='Tool', var_name='Coverage', value_name='Value')

    # Remove the 'cov' prefix from the Coverage column
    df_melted['Coverage'] = df_melted['Coverage'].str.replace('cov', '')

    # Convert 'Coverage' to a numerical value for proper plotting
    df_melted['Coverage'] = pd.to_numeric(df_melted['Coverage'], errors='coerce')

    # Define the custom color palette
    custom_palette = ['#d46014', '#ddcd3d', '#064b76ff', '#63bdf6ff', '#b54582']

    # Plot
    plt.figure(figsize=(5,4))

    # Plot the data for each tool with the custom palette
    for i, tool in enumerate(df_melted['Tool'].unique()):
        tool_data = df_melted[df_melted['Tool'] == tool]
        plt.plot(tool_data['Coverage'], tool_data['Value'], marker='o', label=tool, color=custom_palette[i % len(custom_palette)])

    # Customize grid with dashed lines
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # Improve plot aesthetics
    sns.despine()
    plt.xlabel('Coverage', fontsize=16)
    plt.ylabel(f'{title}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit the legend

    # Set larger tick labels
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  

    # Set dynamic y-axis range based on the maximum value
    max_value = df_melted['Value'].max()
    if max_value < 1:
        plt.ylim(0, 1.05)  
    
    # Add the custom legend to the plot
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=16)

    # Save plot if save_directory is provided
    if save_directory:
        os.makedirs(save_directory, exist_ok=True)
        save_path = os.path.join(save_directory, f'{title.lower().replace(" ", "_")}_plot.png')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()
