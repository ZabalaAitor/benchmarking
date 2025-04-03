import pysam
import pandas as pd

def process_circle_matrix(matrix_file, bam_file, output_file, use_chr_prefix=False):
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
            CJ1 = {read.query_name for read in bamobject.fetch(chrom, start - 1, start + 1)}
            CJ1_n = len(set(CJ1))

            # Fetch CJ2 Reads (end-1, end, end+1)
            CJ2 = {read.query_name for read in bamobject.fetch(chrom, end - 1, end + 1)}
            CJ2_n = len(set(CJ2))

            # Sum of CJ1 and CJ2
            CJ_unique_reads = len(CJ1 | CJ2)
            CJ_reads = CJ1_n + CJ2_n

            # Fetch total reads in the region (start-1 to end+1)
            total_reads = len({read.query_name for read in bamobject.fetch(chrom, start - 1, end + 1)})

            # Compute ratios
            ratio_CJ1 = CJ1_n / CJ_reads if CJ_reads > 0 else 0
            ratio_CJ2 = CJ2_n / CJ_reads if CJ_reads > 0 else 0

            # Determine max and min CJ
            max_CJ = max(ratio_CJ1, ratio_CJ2)
            min_CJ = min(ratio_CJ1, ratio_CJ2)

            # Compute CJ max and min difference
            diff_CJ = max_CJ - min_CJ

            # Compute additional ratios
            ratio = CJ_unique_reads / total_reads if total_reads > 0 else 0
            circle_length = end - start + 1

            # Store results
            results.append([circle, circle_length, tools_detected, CJ1_n, CJ2_n, CJ_unique_reads, CJ_reads,  min_CJ, max_CJ, diff_CJ, total_reads, ratio])

        except ValueError as e:
            # Handle the error if there's an invalid chromosome or other issue
            print(f"Skipping invalid chromosome: {chrom} or region {circle}. Error: {e}")
            continue

    # Close BAM file
    bamobject.close()

    # Create output DataFrame
    output_df = pd.DataFrame(results, columns=[
        "Circle", "Length", "Tools", "CJ1 Reads", "CJ2 Reads", "CJ Unique Reads", "CJ Reads", "min CJ", "max CJ", "diff CJ", "Total Reads", "Ratio"
    ])

    # Save to CSV
    output_df.to_csv(output_file, index=False)