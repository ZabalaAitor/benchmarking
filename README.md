# Comparative analysis of eccDNA and circRNA tools shows increased accuracy of tool combination

This repository contains the code and notebooks used in the study **"Comparative analysis of eccDNA and circRNA tools shows increased accuracy of tool combination"**.

## 📁 Repository Structure

- `benchmarking.ipynb`: Main notebook performing the analysis and generating plots.
- `functions/`: Python scripts with modular functions for:

## 📥 How to Use This Repository

### 1. Clone the Repository

```bash
git clone https://github.com/ZabalaAitor/benchmarking
cd benchmarking
```

### 2. Download input data from Zenodo:

🔗 https://zenodo.org/record/15783793

After downloading, extract the contents and organize them as follows:

- data/ → should contain all BED files with raw and filtered outputs from circRNA and eccDNA detection tools.

- repeatmasker/ → must contain the txt file required for repeat element annotation analysis.

- genomic_elements/ → must contain the GTF files and other genome annotation data required for genomic annotation analysis.

These three folders (data/, repeatmasker/ and genomic_elements/) are essential for running the pipeline.

## 📊 Results

The results of the benchmarking analysis are available at:  
🔗 [https://zenodo.org/records/15783795](https://zenodo.org/records/15783795)
