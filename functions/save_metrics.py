import os
import pandas as pd

def load_stats(folder, method, benchmarking_base_path):
    """
    Loads precision, recall, fscore, truepositives, falsenegatives, and falsepositives CSVs
    based on the folder and method (tool) being processed.
    """
    
    stats_path = os.path.join(benchmarking_base_path, folder, 'statistics')

    true_positives = pd.read_csv(os.path.join(stats_path, "truepositives.csv"))
    false_negatives = pd.read_csv(os.path.join(stats_path, "falsenegatives.csv"))
    false_positives = pd.read_csv(os.path.join(stats_path, "falsepositives.csv"))
    precision = pd.read_csv(os.path.join(stats_path, "precision.csv"))
    recall = pd.read_csv(os.path.join(stats_path, "recall.csv"))
    fscore = pd.read_csv(os.path.join(stats_path, "fscore.csv"))

    for df in [true_positives, false_negatives, false_positives, precision, recall, fscore]:
        df.columns = df.columns.str.replace('cov', '')

    return true_positives, false_negatives, false_positives, precision, recall, fscore

def process_and_save_metrics(base_path, benchmarking_base_path, tools, cov_values, folders, output_path='analysis_results.xlsx'):
    """
    Process benchmarking data and save the results to an Excel file with multiple sheets.
    """
    counts_data = []
    truepositives_data = []
    falsenegatives_data = []
    falsepositives_data = []
    precision_data = []
    recall_data = []
    fscore_data = []

    for tool in tools:
        for cov in cov_values:
            row_counts = {"Tool": tool, "Coverage": cov}
            truepositives_row = {"Tool": tool, "Coverage": cov}
            falsenegatives_row = {"Tool": tool, "Coverage": cov}
            falsepositives_row = {"Tool": tool, "Coverage": cov}
            precision_row = {"Tool": tool, "Coverage": cov}
            recall_row = {"Tool": tool, "Coverage": cov}
            fscore_row = {"Tool": tool, "Coverage": cov}

            for folder in folders:
                # Count rows in the BED file
                bed_file = os.path.join(base_path, folder, tool, f"cov{cov}_{tool}.bed")
                row_count = 0
                if os.path.exists(bed_file):
                    with open(bed_file, "r") as f:
                        row_count = sum(1 for line in f)
                row_counts[folder.capitalize()] = row_count

                # Load statistics for current folder
                true_positives, false_negatives, false_positives, precision, recall, fscore = load_stats(folder, tool, benchmarking_base_path)

                if str(cov) in true_positives.columns:
                    tp_val = true_positives.loc[true_positives['Tool'] == tool, str(cov)]
                    if not tp_val.empty:
                        truepositives_row[folder.capitalize()] = tp_val.values[0]

                if str(cov) in false_negatives.columns:
                    fn_val = false_negatives.loc[false_negatives['Tool'] == tool, str(cov)]
                    if not fn_val.empty:
                        falsenegatives_row[folder.capitalize()] = fn_val.values[0]

                if str(cov) in false_positives.columns:
                    fp_val = false_positives.loc[false_positives['Tool'] == tool, str(cov)]
                    if not fp_val.empty:
                        falsepositives_row[folder.capitalize()] = fp_val.values[0]

                if str(cov) in precision.columns:
                    prec_val = precision.loc[precision['Tool'] == tool, str(cov)]
                    if not prec_val.empty:
                        precision_row[folder.capitalize()] = prec_val.values[0]

                if str(cov) in recall.columns:
                    rec_val = recall.loc[recall['Tool'] == tool, str(cov)]
                    if not rec_val.empty:
                        recall_row[folder.capitalize()] = rec_val.values[0]

                if str(cov) in fscore.columns:
                    fs_val = fscore.loc[fscore['Tool'] == tool, str(cov)]
                    if not fs_val.empty:
                        fscore_row[folder.capitalize()] = fs_val.values[0]

            counts_data.append(row_counts)
            truepositives_data.append(truepositives_row)
            falsenegatives_data.append(falsenegatives_row)
            falsepositives_data.append(falsepositives_row)
            precision_data.append(precision_row)
            recall_data.append(recall_row)
            fscore_data.append(fscore_row)

    # Save all metrics to Excel
    with pd.ExcelWriter(output_path) as writer:
        pd.DataFrame(counts_data).to_excel(writer, sheet_name="counts", index=False)
        pd.DataFrame(truepositives_data).to_excel(writer, sheet_name="truepositives", index=False)
        pd.DataFrame(falsenegatives_data).to_excel(writer, sheet_name="falsenegatives", index=False)
        pd.DataFrame(falsepositives_data).to_excel(writer, sheet_name="falsepositives", index=False)
        pd.DataFrame(precision_data).to_excel(writer, sheet_name="precision", index=False)
        pd.DataFrame(recall_data).to_excel(writer, sheet_name="recall", index=False)
        pd.DataFrame(fscore_data).to_excel(writer, sheet_name="fscore", index=False)