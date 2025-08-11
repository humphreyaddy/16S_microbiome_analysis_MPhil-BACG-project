# Load required packages
library(tidyverse)
library(networkD3)
library(RColorBrewer)
library(scales)

# === File Paths ===
taxonomy_file <- "/home/humphrey/Desktop/maize_data/TACC/qiime_tax_metadata.csv"
contrib_file <- "/home/humphrey/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"
ko_annot_file <- "/home/humphrey/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz"

# === Read Data ===

# Taxonomy metadata (genus-level)
tax_df <- read_csv(taxonomy_file, show_col_types = FALSE)

# KO contribution table
contrib_df <- read_tsv(contrib_file)

# KO annotation table (ID + Description)
ko_annot <- read_tsv(ko_annot_file, col_names = c("function", "Description"), show_col_types = FALSE)

# === Preprocessing ===

# Melt long
df_long <- contrib_df |>
  pivot_longer(-function, names_to = "Sample_Taxon", values_to = "taxon_function_abun") |>
  separate(Sample_Taxon, into = c("Sample", "Taxon"), sep = "\\|") |>
  filter(taxon_function_abun > 0)

# Join Taxonomy
df_long <- left_join(df_long, tax_df, by = c("Taxon" = "Feature.ID"))

# Use Genus or fallback to Taxon name
df_long$Genus <- ifelse(is.na(df_long$Genus), df_long$Taxon, df_long$Genus)

# Merge KO annotation
df_long <- left_join(df_long, ko_annot, by = "function")

# Combine ID + Description
df_long$FunctionLabel <- paste0(df_long$function, " - ", df_long$Description)

# === Aggregate Function Abundance ===
function_total <- df_long |>
  group_by(function) |>
  summarise(total_abundance = sum(taxon_function_abun), .groups = "drop") |>
  arrange(desc(total_abundance))

# Select top 10 most abundant KOs
top_functions <- function_total$function[1:10]

# Filter to top 10
df_top <- df_long |> filter(function %in% top_functions)

# === Aggregate per Taxon → Function → Sample ===
df_agg <- df_top |>
  group_by(Genus, FunctionLabel, Sample) |>
  summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop")

# === Create Nodes ===
taxa_nodes <- unique(df_agg$Genus)
function_nodes <- unique(df_agg$FunctionLabel)
sample_nodes <- unique(df_agg$Sample)

nodes <- tibble(name = c(taxa_nodes, function_nodes, sample_nodes))

# Node indices
taxa_idx <- setNames(seq_along(taxa_nodes) - 1, taxa_nodes)
function_idx <- setNames(seq_along(function_nodes) - 1 + length(taxa_nodes), function_nodes)
sample_idx <- setNames(seq_along(sample_nodes) - 1 + length(taxa_nodes) + length(function_nodes), sample_nodes)

# === Create Links ===

# Taxon → Function
links_tf <- df_agg |>
  group_by(Genus, FunctionLabel) |>
  summarise(value = sum(taxon_function_abun), .groups = "drop") |>
  mutate(source = taxa_idx[Genus],
         target = function_idx[FunctionLabel])

# Function → Sample
links_fs <- df_agg |>
  group_by(FunctionLabel, Sample) |>
  summarise(value = sum(taxon_function_abun), .groups = "drop") |>
  mutate(source = function_idx[FunctionLabel],
         target = sample_idx[Sample])

# Combine all links
links <- bind_rows(links_tf, links_fs)

# === Assign Consistent Colors by Taxon ===
n_taxa <- length(taxa_nodes)
taxon_colors <- hue_pal()(n_taxa)
names(taxon_colors) <- taxa_nodes

# Color for each node (taxon → function → sample)
node_colors <- rep("gray", nrow(nodes))
for (i in seq_along(taxa_nodes)) {
  # Color taxon node
  node_colors[taxa_idx[taxa_nodes[i]] + 1] <- taxon_colors[i]
  
  # Color function nodes connected to taxon
  related_functions <- df_agg |> filter(Genus == taxa_nodes[i]) |> pull(FunctionLabel) |> unique()
  for (func in related_functions) {
    func_index <- function_idx[func] + 1
    node_colors[func_index] <- taxon_colors[i]
  }
  
  # Color sample nodes connected to taxon (optional: same color through)
  related_samples <- df_agg |> filter(Genus == taxa_nodes[i]) |> pull(Sample) |> unique()
  for (samp in related_samples) {
    samp_index <- sample_idx[samp] + 1
    node_colors[samp_index] <- taxon_colors[i]
  }
}

# === Plot Sankey ===
sankeyNetwork(Links = links,
              Nodes = nodes,
              Source = "source",
              Target = "target",
              Value = "value",
              NodeID = "name",
              fontSize = 12,
              nodeWidth = 30,
              colourScale = JS(paste0("d3.scaleOrdinal().range([\"", paste(node_colors, collapse = "\",\""), "\"])")),
              LinkGroup = NULL,
              NodeGroup = NULL)
