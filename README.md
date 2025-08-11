# 16S_microbiome_analysis_MPhil-BACG-project


This repository contains scripts, outputs, and documentation related to my MPhil thesis project in the field of Biodata Analytics and Computational Genomics, with a focus on microbial community profiling and functional inference in traditional African fermented maize-based foods. The analysis was conducted using 16S rRNA gene amplicon sequencing data, with bioinformatics pipelines implemented in QIIME2 and downstream statistical analyses and visualisations performed in R.

## 📘 Project Overview

The project investigates the taxonomic composition, alpha and beta diversity, and functional capacities of bacterial communities associated with African fermented maize-based foods. Amplicon Sequence Variant (ASV)-based analysis was performed using DADA2 in QIIME2, followed by taxonomic classification with the SILVA 138 database. Predictive functional profiling was carried out using PICRUSt2, and extensive visualisations were generated using R (v4.4.2) and various tidyverse and microbiome analysis packages.

## 🧪 Main Analyses and Features

- **Taxonomic Profiling**: Barplots and heatmaps of bacterial composition at multiple taxonomic levels
- **Alpha Diversity**: Richness and evenness indices with visualizations (Shannon, Chao1, Observed ASVs)
- **Beta Diversity**: Ordination plots (PCoA) using Bray–Curtis dissimilarity metrics
- **Rarefaction Analysis**: Assessment of sequencing depth sufficiency
- **Functional Prediction**: Pathway and enzyme-level profiling using PICRUSt2 outputs (KEGG, EC, MetaCyc)
- **Comparative Analysis**: Insights into microbial and functional differences between Ga and Fante Kenkey
- **Publication-Ready Figures**: All figures are generated using `ggplot2`, `phyloseq`, `qiime2R`, and other R packages, and are available in high-resolution PNG format

## 📁 Repository Structure
📦16S_microbiome_analysis_MPhil-BACG-project/
├── data/           # Input datasets (QIIME2 artifacts, metadata, etc.)
├── scripts/        # R scripts for visualization and statistical analysis
├── figures/        # Output plots used in the thesis
├── outputs/        # CSV files of computed metrics and intermediate outputs
├── README.md       # Project overview and documentation
└── requirements.txt # R package dependencies (optional)




## 📊 Key Tools and Packages

- **QIIME2 (v2024.2)** – ASV inference, taxonomic classification
- **PICRUSt2 (v2.5.2)** – Functional metagenome prediction
- **R (v4.4.2)** – Statistical analysis and visualization
    - `phyloseq`, `qiime2R`, `microbiome`, `vegan`, `pheatmap`, `ggplot2`, `dplyr`, `tidyr`, `readr`

## 📌 How to Reproduce

To reproduce the analyses:

1. Clone the repository:
   ```bash
   git clone https://github.com/humphreyaddy/16S_microbiome_analysis_MPhil-BACG-project.git
   cd 16S_microbiome_analysis_MPhil-BACG-project

