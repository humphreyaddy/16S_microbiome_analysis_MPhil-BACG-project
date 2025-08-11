#!/bin/bash
#SBATCH --job-name=MPA_H
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/picrust2_output/picrust_logs/picrust_maize_%j.out
#SBATCH --error=/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/picrust2_output/picrust_logs/picrust_maize_%j.err

# -------------------------------
# Load Conda and activate picrust2 
# -------------------------------
echo "Activating picrust2  environment..."
source /work/10451/humphrey_addy/miniconda3/etc/profile.d/conda.sh
conda activate picrust2 

# -------------------------------
# Set working directory
# -------------------------------
WORKDIR="/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/picrust2_output"
cd "$WORKDIR"

echo "Running on host: $HOSTNAME"
echo "Working directly in: $WORKDIR"

# -------------------------------
# Run picrust2 pipeline
# -------------------------------
echo "Starting picrust2 pipeline..."
bash picrust_maize_script.sh

echo "picrust2 job completed successfully."

