# This function depends on the external tool `ecc_finder` for circular DNA detection.
# To use this functionality, you must manually clone the ecc_finder repository, as it is not included in this project:
#     git clone https://github.com/njaupan/ecc_finder.git
#
# Make sure to run any `ecc_finder.py` commands from within that directory,
# or provide the full path to the script when calling it from this pipeline.
#
# Reference:
# Panpan Z, Haoran P, Christel L, Etienne B, Marie M.
# ecc_finder: a robust and accurate tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data.
# Frontiers in Plant Science, 2021.
# https://www.frontiersin.org/articles/10.3389/fpls.2021.743742/

import os
import csv
import argparse

def run_ecc_finder(reference, reference_idx, pairs, output_dir, aligner):
    for query_fq1, query_fq2 in pairs:
        sample_name = os.path.basename(query_fq1).split('.')[0]
        output_prefix = os.path.join(output_dir, sample_name)
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)
        
        # Construct command based on aligner
        if aligner == 'bwa':
            cmd = [
                'python', 'ecc_finder.py', 'map-sr',
                reference, query_fq1, query_fq2,
                '-r', reference_idx,
                '-o', output_prefix
            ]
        elif aligner == 'minimap2':
            cmd = [
                'python', 'ecc_finder.py', 'map-sr',
                reference_idx, query_fq1, query_fq2,
                '-r', reference,
                '-o', output_prefix,
                '--aligner', 'minimap2'
            ]
        else:
            raise ValueError(f"Unsupported aligner: {aligner}")

        print('Running command:', ' '.join(cmd))
        os.system(' '.join(cmd))

def load_pairs_from_csv(samplesheet_path):
    pairs = []
    with open(samplesheet_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            fq1 = row['fq1'].strip()
            fq2 = row['fq2'].strip()
            pairs.append((fq1, fq2))
    return pairs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ecc_finder on multiple samples.")
    parser.add_argument('--reference', required=True, help="Path to the reference FASTA")
    parser.add_argument('--reference_idx', required=True, help="Path to the reference index")
    parser.add_argument('--samplesheet', required=True, help="CSV file with columns of trrmmed fq1,fq2 from nf-core/circdna")
    parser.add_argument('--output_dir', required=True, help="Directory to write output files")
    parser.add_argument('--aligner', choices=['bwa', 'minimap2'], required=True, help="Aligner to use")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    pairs = load_pairs_from_csv(args.samplesheet)
    run_ecc_finder(args.reference, args.reference_idx, pairs, args.output_dir, args.aligner)