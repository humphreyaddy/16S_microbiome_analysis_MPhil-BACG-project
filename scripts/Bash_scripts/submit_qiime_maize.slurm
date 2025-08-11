#!/bin/bash
#SBATCH --job-name=MQA_H
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/logs/qiime_maize_%j.out
#SBATCH --error=/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/logs/qiime_maize_%j.err

# -------------------------------
# Load Conda and activate QIIME2
# -------------------------------
echo "Activating QIIME2 environment..."
source /work/10451/humphrey_addy/miniconda3/etc/profile.d/conda.sh
conda activate qiime2

# -------------------------------
# Set working directory
# -------------------------------
WORKDIR="/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data"
cd "$WORKDIR"

echo "Running on host: $HOSTNAME"
echo "Working directly in: $WORKDIR"

# -------------------------------
# Run QIIME2 pipeline
# -------------------------------
echo "Starting QIIME2 pipeline..."
bash qiime_script_maize.sh

echo "QIIME2 job completed successfully."

