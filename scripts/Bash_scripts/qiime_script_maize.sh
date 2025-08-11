#!/bin/bash

# =========================================
# My QIIME 2 Analysis for 16S PE Data
# =========================================

# Set number of threads
THREADS=16

# Set working directory
WORKDIR="/work/10451/humphrey_addy/ls6/16S_raw_data/maize_only_16S_data/"
cd "$WORKDIR"

# Create output directory
mkdir -p output


echo "Using $THREADS threads/jobs for some computationally intensive steps."

# -------------------------------------------------
# 1. Import Raw Paired-End Sequences
# -------------------------------------------------
echo "Importing raw paired-end sequences..."

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path output/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

# -------------------------------------------------
# 2. Demultiplex Sequences
# -------------------------------------------------
echo "Demultiplexing sequences..."

qiime demux summarize \
  --i-data output/demux-paired-end.qza \
  --o-visualization output/demux-paired-end.qzv

# -------------------------------------------------
# 3. Quality Control and Feature Table Construction
# -------------------------------------------------
echo "Performing quality control with DADA2..."

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs output/demux-paired-end.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250  # Match forward and reverse truncation
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-threads "$THREADS" \
  --o-representative-sequences output/rep-seqs.qza \
  --o-table output/table.qza \
  --o-denoising-stats output/denoising-stats.qza

# Summarize DADA2 stats
qiime metadata tabulate \
  --m-input-file output/denoising-stats.qza \
  --o-visualization output/denoising-stats.qzv

# -------------------------------------------------
# 4. Summarize FeatureTable and FeatureData
# -------------------------------------------------
echo "Summarizing FeatureTable and FeatureData..."

qiime feature-table summarize \
  --i-table output/table.qza \
  --o-visualization output/table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data output/rep-seqs.qza \
  --o-visualization output/rep-seqs.qzv

# -------------------------------------------------
# 5. Taxonomic Classification with Pretrained Classifier
# -------------------------------------------------
echo "Classifying taxonomy using pretrained SILVA classifier..."

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads output/rep-seqs.qza \
  --o-classification output/taxonomy.qza

# Summarize taxonomy
qiime metadata tabulate \
  --m-input-file output/taxonomy.qza \
  --o-visualization output/taxonomy.qzv

qiime taxa barplot \
  --i-table output/table.qza \
  --i-taxonomy output/taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization output/taxa-bar-plots.qzv

# -------------------------------------------------
# 6. Phylogenetic Tree Construction
# -------------------------------------------------
echo "Constructing phylogenetic tree..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences output/rep-seqs.qza \
  --p-n-threads "$THREADS" \
  --o-alignment output/aligned-rep-seqs.qza \
  --o-masked-alignment output/masked-aligned-rep-seqs.qza \
  --o-tree output/unrooted-tree.qza \
  --o-rooted-tree output/rooted-tree.qza

# -------------------------------------------------
# 7. Export to Phyloseq-Compatible Format
# -------------------------------------------------
echo "Exporting to Phyloseq-compatible format..."

mkdir -p output/phyloseq

# Export feature table
qiime tools export \
  --input-path output/table.qza \
  --output-path output/phyloseq

# Convert BIOM to TSV
biom convert \
  -i output/phyloseq/feature-table.biom \
  -o output/phyloseq/table.tsv \
  --to-tsv

# Clean up TSV
sed -i '1d' output/phyloseq/table.tsv
sed -i 's/#OTU ID//' output/phyloseq/table.tsv

# Export representative sequences
qiime tools export \
  --input-path output/rep-seqs.qza \
  --output-path output/phyloseq

echo "Results saved in: $WORKDIR/output"

