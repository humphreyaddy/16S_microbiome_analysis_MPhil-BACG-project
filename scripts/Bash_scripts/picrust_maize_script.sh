#! /bin/bash

set -x


# Set number of threads
THREADS=32

# Set working directory
INPUT_dir="/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/output/phyloseq/"
OUTPUT_dir="/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/picrust2_output"


picrust2_pipeline.py -s "$INPUT_dir"/dna-sequences.fasta -i "$INPUT_dir"/feature-table.biom -o "$OUTPUT_dir"/picrust2_final_results -p "$THREADS" --stratified


