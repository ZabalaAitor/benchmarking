# Comparative analysis of circRNA and eccDNA tools shows increased accuracy of tool combination

This repository contains the code and notebooks used in the study **"Comparative analysis of circRNA and eccDNA tools shows increased accuracy of tool combination"**.

## 📁 Repository Structure

- `benchmarking.ipynb`: Main notebook performing the analysis and generating plots.
- `functions/`: Python scripts with modular functions for:

## 📥 How to Use This Repository

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/circular-tools-benchmarking.git
cd benchmarking
```

### 2. Download input data from Zenodo:

🔗 https://zenodo.org/record/1234567 

After downloading, extract the contents and organize them as follows:

- data/ → should contain all BED files with raw and filtered outputs from circRNA and eccDNA detection tools.

- genomic_elements/ → must contain the GTF files and other genome annotation data required for genomic annotation analysis.

These two folders (data/ and genomic_elements/) are essential for running the pipeline.
