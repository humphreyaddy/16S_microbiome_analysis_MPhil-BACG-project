# Load required libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(readr)
library(ggfittext)
library(grid)

#setwd("/home/humphrey/Desktop/maize_data/TACC/humphrey_picrust2_results/R visualisation/new/one/")


# Step 1: Load KO prediction matrix
ko_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ko_data <- read_tsv(ko_file)

# Step 2: Load KO annotations
ko_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz"
ko_annot <- read_tsv(ko_annot_file, col_names = FALSE)
colnames(ko_annot) <- c("ko_id", "description")
ko_annot$ko_id <- str_replace(ko_annot$ko_id, "ko:", "")

# Step 3: Clean KO IDs and keep only KO IDs (not descriptions) as identifiers
ko_data <- ko_data %>%
  rename(ko_id = `function`) %>%
  mutate(ko_id = str_replace(ko_id, "ko:", ""))

# Step 4: Aggregate by KO ID
ko_data_aggregated <- ko_data %>%
  group_by(ko_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Step 5: Convert to matrix (KO x Sample)
ko_mat <- ko_data_aggregated %>%
  column_to_rownames("ko_id") %>%
  as.matrix()

# Step 6: Transpose matrix so samples are rows
ko_mat_t <- t(ko_mat)

# Step 7: Select top 10 KO IDs by total abundance across all samples
top_ko_ids <- names(sort(colSums(ko_mat_t), decreasing = TRUE)[1:10])

# Step 8: Subset matrix to top 10 KOs
top_ko_mat <- ko_mat_t[, top_ko_ids]

# Step 9: Normalize each sample (row-wise z-score)
scaled_mat <- t(scale(t(top_ko_mat)))

# Step 10: Transpose for heatmap display (KO IDs as rows)
transposed_mat <- t(scaled_mat)

# Step 11: Annotate KO IDs with descriptions
top_ko_annot <- ko_annot %>% filter(ko_id %in% rownames(transposed_mat))
desc_labels <- setNames(top_ko_annot$description, top_ko_annot$ko_id)

# Step 12: Create custom row labels (for legends or tooltips)
rownames(transposed_mat) <- paste0(rownames(transposed_mat), " | ", desc_labels[rownames(transposed_mat)])

# Step 13: Plot and export as SVG
svg("Heatmap_Top15_KO_IDs.svg", width = 12, height = 8)
pheatmap(transposed_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         display_numbers = matrix(sprintf("%.2f", transposed_mat), 
                                  nrow = nrow(transposed_mat), 
                                  dimnames = dimnames(transposed_mat)),
         number_color = "black",
         fontsize_number = 8,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Top 10 KO IDs by Total Abundance")
dev.off()




# Step 1: Load KO prediction matrix
ko_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ko_data <- read_tsv(ko_file)

# Step 2: Load KO name mapping
ko_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz"
ko_annot <- read_tsv(ko_annot_file, col_names = FALSE)
colnames(ko_annot) <- c("ko_id", "description")
ko_annot$ko_id <- str_replace(ko_annot$ko_id, "ko:", "")

# Step 3: Reshape KO data (long format)
ko_long <- ko_data %>%
  rename(ko_id = `function`) %>%
  mutate(ko_id = str_replace(ko_id, "ko:", "")) %>%
  pivot_longer(cols = -ko_id, names_to = "Sample", values_to = "Abundance")

# Step 4: Identify Top 10 KO IDs by Total Abundance
top_ko_ids <- ko_long %>%
  group_by(ko_id) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 15) %>%
  pull(ko_id)

# Step 5: Filter to top 10 KO IDs and annotate
ko_top_long <- ko_long %>%
  filter(ko_id %in% top_ko_ids) %>%
  left_join(ko_annot, by = "ko_id") %>%
  mutate(description = ifelse(is.na(description), ko_id, paste0(ko_id, ": ", description)))  # fallback

# Step 6: Aggregate for stacked bar plot
ko_top_long_sum <- ko_top_long %>%
  group_by(Sample, description) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  arrange(Sample, desc(Abundance)) %>%
  group_by(Sample) %>%
  mutate(cum_sum = cumsum(Abundance) - Abundance / 2) %>%
  ungroup()

# Step 7: Compute sample totals for bar-top annotation
bar_sums <- ko_top_long_sum %>%
  group_by(Sample) %>%
  summarise(Total = sum(Abundance), .groups = "drop")

# Step 8: Set custom color palette
#bar_colors <- brewer.pal(n = 15, "Paired")
bar_colors <- c(
  # Original 8 from Dark2
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666",
  
  # 7 handpicked complementary and visually distinct colors
  "#8DA0CB",  # pastel blue
  "#FC8D62",  # coral
  "#FFD92F",  # golden yellow
  "#A6D854",  # lime green
  "#E78AC3",  # pinkish purple
  "#A6CEE3",  # soft sky blue
  "#B3B3B3"   # soft gray
)






ko_top_long_sum$description <- factor(ko_top_long_sum$description,
                                      levels = rev(unique(ko_top_long_sum$description)))

# Step 9: Build stacked barplot
p <- ggplot(ko_top_long_sum, aes(x = Sample, y = Abundance, fill = description)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cum_sum, label = round(Abundance, 1)),
            size = 3, fontface = "bold", color = "black")+
  geom_text(data = bar_sums,
            aes(x = Sample, y = Total, label = round(Total, 1)),
            vjust = -0.5, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = bar_colors) +
  theme_minimal() +
  labs(title = "Top 10 KO IDs Across Samples",
       y = "Predicted Abundance",
       x = "Sample",
       fill = "KO Function (ID: Description)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Step 10: Save as SVG at 330 DPI
ggsave("1Top_10_KO_IDs_Across_Samples.svg", plot = p,
       dpi = 330, width = 10, height = 6, units = "in")





# Step 1: Reload KO prediction data
ko_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ko_data <- readr::read_tsv(ko_file)

# Step 2: Reload KO name annotation
ko_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz"
ko_annot <- readr::read_tsv(ko_annot_file, col_names = FALSE)
colnames(ko_annot) <- c("ko_id", "description")
ko_annot$ko_id <- stringr::str_replace(ko_annot$ko_id, "ko:", "")

# Step 3: Annotate KO table
ko_annotated <- ko_data %>%
  dplyr::rename(ko_id = `function`) %>%
  dplyr::mutate(ko_id = stringr::str_replace(ko_id, "ko:", "")) %>%
  dplyr::left_join(ko_annot, by = "ko_id")

# Step 4: Reshape to long format
ko_long <- ko_annotated %>%
  tidyr::pivot_longer(cols = -c(ko_id, description),
                      names_to = "Sample",
                      values_to = "Abundance")

# Step 5: Identify top 10 most abundant KO functions overall
top_kos <- ko_long %>%
  dplyr::group_by(description) %>%
  dplyr::summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::arrange(desc(total_abundance)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::pull(description)

# Step 6: Filter data for top KOs and handle NA
ko_top_long <- ko_long %>%
  dplyr::mutate(description = ifelse(is.na(description), "Unannotated", description)) %>%
  dplyr::filter(description %in% top_kos)

# Step 7: Generate color palette for top 10 KOs
#bar_colors <- c(RColorBrewer::brewer.pal(8, "Paired"),
                RColorBrewer::brewer.pal(12, "Set3"))[1:length(top_kos)]


bar_colors <- c(
  # Original 8 from Dark2
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666",
  
  # 7 handpicked complementary and visually distinct colors
  "#8DA0CB",  # pastel blue
  "#FC8D62",  # coral
  "#FFD92F",  # golden yellow
  "#A6D854",  # lime green
  "#E78AC3",  # pinkish purple
  "#A6CEE3",  # soft sky blue
  "#B3B3B3"   # soft gray
)



# Step 8: Aggregate KO values per sample and description
ko_top_long_sum <- ko_top_long %>%
  dplyr::group_by(Sample, description) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

# Step 9: Compute label positions within stacked bars
ko_top_long_sum <- ko_top_long_sum %>%
  dplyr::arrange(Sample, desc(description)) %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(cum_sum = cumsum(Abundance) - Abundance / 2) %>%
  dplyr::ungroup()

# Step 10: Calculate total per sample
bar_sums <- ko_top_long_sum %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(Total = sum(Abundance), .groups = "drop")

# Step 11: Build the annotated plot
p <- ggplot2::ggplot(ko_top_long_sum, aes(x = Sample, y = Abundance, fill = description)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cum_sum, label = ifelse(Abundance > 5, round(Abundance, 1), "")),
            color = "white", size = 3, fontface = "bold", vjust = 0.5) +
  geom_text(data = bar_sums,
            aes(x = Sample, y = Total, label = round(Total, 1)),
            vjust = -0.5, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = bar_colors) +
  theme_minimal() +
  labs(title = "Top 10 Predicted KO Functions Across Samples",
       y = "Predicted Abundance",
       x = "Sample",
       fill = "KO Function") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
p
# Step 12: Save the plot
ggplot2::ggsave("Top 10 Predicted KO Functions Across Samples.svg", plot = p,
                dpi = 330, width = 10, height = 6, units = "in")






library(stringr)
# Step 1: Load KO contribution table
ko_contrib_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"
ko_contrib <- read_tsv(ko_contrib_file)

# Step 2: Load KO annotations
ko_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz"
ko_annot <- read_tsv(ko_annot_file, col_names = FALSE)
colnames(ko_annot) <- c("ko_id", "description")
ko_annot$ko_id <- str_replace(ko_annot$ko_id, "ko:", "")

# Step 3: Set KO ID of interest
target_ko_id <- "K06147"

# Step 4: Normalize KO IDs in the contribution table
ko_contrib <- ko_contrib %>%
  mutate(`function` = str_replace(`function`, "ko:", ""))

# Step 5: Filter data for the target KO
target_data <- ko_contrib %>%
  filter(`function` == target_ko_id)

# Step 6: Retrieve the KO description
target_description <- ko_annot %>%
  filter(ko_id == target_ko_id) %>%
  pull(description) %>%
  .[1]

# Step 7: Aggregate top contributors
contrib_summary <- target_data %>%
  group_by(taxon) %>%
  summarise(Total_Contribution = sum(norm_taxon_function_contrib, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_Contribution)) %>%
  slice_head(n = 15)

# Step 8: Plot
# Define custom 15-color palette
bar_colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#8DA0CB", "#FC8D62", "#FFD92F", "#A6D854",
  "#E78AC3", "#A6CEE3", "#B3B3B3"
)

# Plot using custom palette
p <- ggplot(contrib_summary,
            aes(x = reorder(taxon, Total_Contribution),
                y = Total_Contribution,
                fill = taxon)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(Total_Contribution, 2)),
            hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = bar_colors) +  # Use custom palette
  theme_minimal() +
  labs(title = paste("Top Contributors to", target_ko_id, "-", target_description),
       x = "Taxon", y = "Total Contribution Score") +
  theme(plot.title = element_text(hjust = 0.5))

# Show the plot
print(p)


# Step 9: Save the plot
filename <- paste0("Top_Contributors_to_", target_ko_id, "_", str_replace_all(target_description, "[^a-zA-Z0-9]+", "_"), ".svg")
ggsave(filename, plot = p, dpi = 330, width = 10, height = 6, units = "in")






# ----------- Part A: Total Predicted KO Functions per Sample -----------

# Step 1: Load KO prediction table
ko_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ko_data <- read_tsv(ko_file)

# Step 2: Remove KO ID column
ko_only <- ko_data %>% select(-`function`)

# Step 3: Count non-zero KO functions per sample
ko_richness <- colSums(ko_only > 0) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  rename(Num_KOs = ".")

# Optional: store sample order for consistency
sample_order <- ko_richness$Sample

# Define consistent sample color mapping
sample_colors <- setNames(
  RColorBrewer::brewer.pal(n = length(sample_order), "Accent")[1:length(sample_order)],
  sample_order
)

scale_fill_manual(values = sample_colors)

p_ko <- ggplot(ko_richness, aes(x = Sample, y = Num_KOs, fill = Sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Num_KOs), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = sample_colors) +
  labs(title = "Total Number of Predicted KO Functions per Sample",
       y = "Number of Functions", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 5: Save KO richness plot
ggsave("Total Number of Predicted KO Functions per Sample.svg", plot = p_ko,
       dpi = 330, width = 10, height = 6, units = "in")
#ggsave("Total Number of Predicted KO Functions per Sample.pdf", plot = p_ko,
#       dpi = 330, width = 10, height = 6, units = "in")

# Optional: export table
write_csv(ko_richness, "KO_Richness_Per_Sample.csv")

# ----------- Part B: NSTI Score per Sample -----------

# Step 1: Load NSTI scores
nsti_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/weighted_nsti.tsv.gz"
nsti <- read_tsv(nsti_file)

# Step 2: Ensure consistent sample order
nsti$sample <- factor(nsti$sample, levels = sample_order)

# Define consistent sample color mapping
sample_colors <- setNames(
  RColorBrewer::brewer.pal(n = length(sample_order), "Accent")[1:length(sample_order)],
  sample_order
)

scale_fill_manual(values = sample_colors)

p_nsti <- ggplot(nsti, aes(x = sample, y = weighted_NSTI, fill = sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(weighted_NSTI, 3)), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = sample_colors) +
  labs(title = "NSTI Score per Sample",
       x = "Sample", y = "Mean NSTI (Prediction Accuracy)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 4: Save NSTI plot
ggsave("NSTI Score per Sample.svg", plot = p_nsti,
       dpi = 330, width = 10, height = 6, units = "in")
#ggsave("NSTI Score per Sample.pdf", plot = p_nsti,
#       dpi = 330, width = 10, height = 6, units = "in")

# Optional: export table
write_csv(nsti, "NSTI_Scores_Per_Sample.csv")

# ----------- Optional: Display plots if running interactively -----------
print(p_ko)
print(p_nsti)






# Step 1: Load KO prediction matrix
ko_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ko_data <- read_tsv(ko_file)

# Step 2: Prepare KO matrix (KO x Sample)
ko_mat <- ko_data %>%
  column_to_rownames("function") %>%
  as.matrix()

# Step 3: Transpose (samples as rows)
ko_mat_t <- t(ko_mat)

# Step 4: Convert to relative abundance (row-wise normalization)
ko_rel <- ko_mat_t / rowSums(ko_mat_t)

# Step 5: Compute Bray-Curtis distance matrix
bray_dist <- vegdist(ko_rel, method = "bray")

# Step 6: Perform PCoA
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(pcoa_df)

# Step 7: Assign consistent colors using Accent palette
sample_names <- pcoa_df$Sample
num_samples <- length(sample_names)
accent_colors <- brewer.pal(n = max(3, min(8, num_samples)), name = "Accent")
names(accent_colors) <- sample_names

# Step 8: Plot PCoA
p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Sample, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel(show.legend = FALSE, size = 4, fontface = "bold") +
  scale_color_manual(values = accent_colors) +
  theme_minimal() +
  labs(title = "PCoA of KO-Based Functional Profiles",
       x = "PCoA Axis 1",
       y = "PCoA Axis 2") +
  theme(legend.position = "none")

# Step 9: Save as high-resolution SVG
ggsave("PCoA of KO-Based Functional Profiles.svg",
       plot = p_pcoa, dpi = 330, width = 10, height = 6, units = "in")

# Step 10: Display
p_pcoa





###Enzyme Commission numbers####)

# Step 1: Load EC prediction table
ec_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ec_data <- read_tsv(ec_file)

# Step 2: Load EC annotation table
ec_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ec_name.txt.gz"
ec_annot <- read_tsv(ec_annot_file, col_names = FALSE)
colnames(ec_annot) <- c("ec_id", "description")

# Step 3: Clean EC data and aggregate by ec_id
ec_data_clean <- ec_data %>%
  rename(ec_id = `function`) %>%
  group_by(ec_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Step 4: Convert to matrix (EC x Sample)
ec_mat <- ec_data_clean %>%
  column_to_rownames("ec_id") %>%
  as.matrix()

# Step 5: Transpose (Samples as rows)
ec_mat_t <- t(ec_mat)

# Step 6: Identify top 10 EC IDs by total abundance across all samples
top_ec_ids <- names(sort(colSums(ec_mat_t), decreasing = TRUE)[1:10])

# Step 7: Subset matrix to top 10 ECs
top_ec_mat <- ec_mat_t[, top_ec_ids]

# Step 8: Row-wise z-score normalization
scaled_ec_mat <- t(scale(t(top_ec_mat)))

# Step 9: Annotate EC IDs with descriptions
top_ec_annot <- ec_annot %>% filter(ec_id %in% top_ec_ids)
desc_labels <- setNames(top_ec_annot$description, top_ec_annot$ec_id)

# Step 10: Rename columns to "EC | Description"
colnames(scaled_ec_mat) <- paste0(colnames(scaled_ec_mat), " | ", desc_labels[colnames(scaled_ec_mat)])

# Step 12 (Updated): Transpose matrix for plotting (Samples as rows, EC IDs as columns)
transposed_ec_mat <- t(scaled_ec_mat)

# Prepare annotation labels for transposed matrix
annotation_labels <- matrix(
  sprintf("%.2f", transposed_ec_mat),
  nrow = nrow(transposed_ec_mat),
  dimnames = dimnames(transposed_ec_mat)
)

# Save transposed heatmap as SVG
svg("Heatmap_Top15_EC_IDs_TRANSPOSED.svg", width = 12, height = 8)
pheatmap(transposed_ec_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = annotation_labels,
         number_color = "black",
         color = colorRampPalette(brewer.pal(9, "PuBuGn"))(100),
         main = "Top 10 Predicted EC Functions Across Samples",
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize_number = 8)
dev.off()






# Step 1: Load EC abundance data
ec_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ec_data <- read_tsv(ec_file)

# Step 2: Load EC descriptions
ec_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ec_name.txt.gz"
ec_annot <- read_tsv(ec_annot_file, col_names = FALSE)
colnames(ec_annot) <- c("ec_id", "description")

# Step 3: Annotate ECs
ec_annotated <- ec_data %>%
  rename(ec_id = `function`) %>%
  left_join(ec_annot, by = "ec_id") %>%
  filter(!is.na(description))

# Step 4: Reshape to long format
ec_long <- ec_annotated %>%
  pivot_longer(cols = -c(ec_id, description),
               names_to = "Sample",
               values_to = "Abundance")

# Step 5: Identify top 10 ECs by total abundance
top_ecs <- ec_long %>%
  group_by(description) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 15) %>%
  pull(description)

# Step 6: Filter for top 10
ec_top_long <- ec_long %>%
  filter(description %in% top_ecs)

# Step 7: Assign color palette
ec_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))[1:length(top_ecs)]

# Step 8a: Compute label positions within each bar segment
ec_top_long <- ec_top_long %>%
  arrange(Sample, desc(description)) %>%
  group_by(Sample) %>%
  mutate(midpoint = cumsum(Abundance) - Abundance / 2) %>%
  ungroup()

# Step 8b: Compute total abundance per bar
bar_sums <- ec_top_long %>%
  group_by(Sample) %>%
  summarise(Total = sum(Abundance), .groups = "drop")

# Step 9: Build stacked barplot
p <- ggplot(ec_top_long, aes(x = Sample, y = Abundance, fill = description)) +
  geom_bar(stat = "identity") +
  # Add exact value to each segment
  geom_text(aes(y = midpoint, label = round(Abundance, 1)),
            size = 3, fontface = "bold", color = "black") +
  # Add total per bar
  geom_text(data = bar_sums,
            aes(x = Sample, y = Total, label = round(Total, 1)),
            vjust = -0.5, size = 3, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = ec_colors) +
  theme_minimal() +
  labs(title = "Top 10 Predicted EC Functions Across Samples",
       y = "Predicted Abundance",
       x = "Sample",
       fill = "EC Function") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
p
# Step 10: Save the plot as SVG
ggsave("Top 10 Predicted EC Functions Across Samples.svg",
       plot = p,
       device = "svg",
       dpi = 330,
       width = 12,
       height = 8)



# Step 1: Load EC contribution data
ec_contrib_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz"
ec_contrib <- read_tsv(ec_contrib_file)

# Step 2: Load EC annotations
ec_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ec_name.txt.gz"
ec_annot <- read_tsv(ec_annot_file, col_names = FALSE)
colnames(ec_annot) <- c("ec_id", "description")

# Step 3: Select target EC ID (modify as needed)
target_ec <- "3.6.4.12"  # EC number for DNA helicase

# Step 4: Filter contribution data for the target EC
ec_target <- ec_contrib %>%
  filter(`function` == target_ec)

# Step 5: Get EC description
ec_desc <- ec_annot %>%
  filter(ec_id == target_ec) %>%
  pull(description)

# Step 6: Aggregate contribution by taxon
ec_contrib_summary <- ec_target %>%
  group_by(taxon) %>%
  summarise(Total_Contribution = sum(norm_taxon_function_contrib, na.rm = TRUE)) %>%
  arrange(desc(Total_Contribution)) %>%
  slice_head(n = 15)

# Step 7: Create barplot
p <- ggplot(ec_contrib_summary, aes(x = reorder(taxon, Total_Contribution),
                                    y = Total_Contribution,
                                    fill = taxon)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(Total_Contribution, 2)),
            hjust = -0.2, size = 3, fontface = "bold") +  # Annotate values
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  labs(title = paste("Top Taxonomic Contributors to EC", target_ec, "-", ec_desc),
       x = "Taxon",
       y = "Contribution Score") +
  theme(plot.title = element_text(face = "bold"))
p
# Step 8: Save plot as .svg
ggsave(paste0("Top Taxonomic Contributors to EC ", target_ec, " - ", ec_desc, ".svg"),
       plot = p,
       device = "svg",
       dpi = 330,
       width = 10,
       height = 6)




# ----------- Part A: Total Predicted EC Functions per Sample -----------

# Load EC abundance data
ec_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ec_data <- read_tsv(ec_file)

# Count number of non-zero EC functions per sample
ec_richness <- ec_data %>%
  select(-`function`) %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Num_ECs")

# Barplot of EC richness
p1 <- ggplot(ec_richness, aes(x = Sample, y = Num_ECs, fill = Sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Num_ECs), vjust = -0.3, fontface = "bold", size = 3) +  # Bold annotations
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Total Number of Predicted EC Functions per Sample",
       x = "Sample",
       y = "Number of EC Functions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))
p1
# Save Part A
ggsave("Total Number of Predicted EC Functions per Sample.svg",
       plot = p1, device = "svg", dpi = 330, width = 10, height = 6)

# ----------- Part B: NSTI Scores per Sample -----------

# Load NSTI scores for EC predictions
nsti_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/weighted_nsti.tsv.gz"
ec_nsti <- read_tsv(nsti_file)

# NSTI barplot
p2 <- ggplot(ec_nsti, aes(x = sample, y = weighted_NSTI, fill = sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(weighted_NSTI, 3)), vjust = -0.3, fontface = "bold", size = 3) +
  scale_fill_brewer(palette = "Accent") +
  labs(title = "NSTI Score per Sample (EC Predictions)",
       x = "Sample",
       y = "Mean NSTI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# Save Part B
ggsave("NSTI Score per Sample (EC Predictions).svg",
       plot = p2, device = "svg", dpi = 330, width = 10, height = 6)
p1
p2






# Step 1: Load EC abundance data
ec_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
ec_data <- read_tsv(ec_file)

# Step 2: Create matrix (ECs as rows, samples as columns)
ec_mat <- ec_data %>%
  column_to_rownames("function") %>%
  as.matrix()

# Step 3: Transpose (samples as rows)
ec_mat_t <- t(ec_mat)

# Step 4: Convert to relative abundance
ec_rel <- ec_mat_t / rowSums(ec_mat_t)

# Step 5: Compute Bray-Curtis distance matrix
bray_ec <- vegdist(ec_rel, method = "bray")

# Step 6: Run PCoA
pcoa_ec <- cmdscale(bray_ec, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_ec$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(pcoa_df)

# Step 7: Define consistent color palette for samples
samples <- pcoa_df$Sample
sample_colors <- setNames(brewer.pal(length(samples), "Accent")[1:length(samples)], samples)

# Step 8: Plot with consistent sample colors
p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Sample, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel(size = 4, fontface = "bold") +
  scale_color_manual(values = sample_colors) +
  theme_minimal() +
  labs(title = "PCoA of EC-Based Functional Profiles",
       x = "PCoA Axis 1",
       y = "PCoA Axis 2") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# Step 9: Save as high-resolution SVG
ggsave("PCoA of EC-Based Functional Profiles.svg",
       plot = p, device = "svg", dpi = 330, width = 8, height = 6)




library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Step 1: Load MetaCyc pathway abundance data
pathway_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz"
pathway_data <- read_tsv(pathway_file)

# Step 2: Load MetaCyc pathway descriptions
pathway_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/metacyc-pwy_name.txt.gz"
pathway_annot <- read_tsv(pathway_annot_file, col_names = FALSE)
colnames(pathway_annot) <- c("pathway_id", "description")

# Step 3: Clean and aggregate (do NOT annotate yet)
pathway_data_clean <- pathway_data %>%
  rename(pathway_id = `pathway`) %>%
  group_by(pathway_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Step 4: Convert to matrix (Pathway x Sample)
pathway_mat <- pathway_data_clean %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

# Step 5: Transpose (Samples as rows)
pathway_mat_t <- t(pathway_mat)

# Step 6: Identify Top 10 Pathway IDs by total abundance across all samples
top_pathway_ids <- names(sort(colSums(pathway_mat_t), decreasing = TRUE)[1:10])

# Step 7: Subset matrix to top 10 pathways
top_path_mat <- pathway_mat_t[, top_pathway_ids]

# Step 8: Row-wise z-score normalization
scaled_path_mat <- t(scale(t(top_path_mat)))

# Step 9: Annotate pathways with descriptions
top_path_annot <- pathway_annot %>% filter(pathway_id %in% colnames(scaled_path_mat))
desc_labels <- setNames(top_path_annot$description, top_path_annot$pathway_id)

# Step 10: Rename columns as "Pathway_ID | Description"
colnames(scaled_path_mat) <- paste0(colnames(scaled_path_mat), " | ", desc_labels[colnames(scaled_path_mat)])

# Step 11: Transpose for heatmap display (samples as rows, pathways as columns)
transposed_mat <- t(scaled_path_mat)

# Step 12: Prepare annotation labels
annotation_labels <- matrix(
  sprintf("%.2f", transposed_mat),
  nrow = nrow(transposed_mat),
  dimnames = dimnames(transposed_mat)
)

# Step 13: Plot and save as SVG
svg("Heatmap_Top15_MetaCyc_Pathways.svg", width = 12, height = 8)
pheatmap(transposed_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = annotation_labels,
         number_color = "black",
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
         main = "Top 10 MetaCyc Pathways by Total Abundance",
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize_number = 8)
dev.off()





library(svglite)
# Step 1: Load pathway abundance data
pathway_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz"
pathway_data <- readr::read_tsv(pathway_file)

# Step 2: Load pathway descriptions
pathway_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/metacyc-pwy_name.txt.gz"
pathway_annot <- readr::read_tsv(pathway_annot_file, col_names = FALSE)
colnames(pathway_annot) <- c("pathway_id", "description")

# Step 3: Merge and clean
pathway_annotated <- pathway_data %>%
  dplyr::left_join(pathway_annot, by = c("pathway" = "pathway_id")) %>%
  dplyr::filter(!is.na(description))

# Step 4: Reshape to long format
pathway_long <- pathway_annotated %>%
  tidyr::pivot_longer(cols = -c(pathway, description),
                      names_to = "Sample",
                      values_to = "Abundance")

# Step 5: Identify top 10 abundant pathways
top_paths <- pathway_long %>%
  dplyr::group_by(description) %>%
  dplyr::summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::arrange(desc(total_abundance)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::pull(description)

# Step 6: Filter to top pathways
top_pathway_long <- pathway_long %>%
  dplyr::filter(description %in% top_paths)

# Step 7: Use same color logic as EC script
path_colors <- c(RColorBrewer::brewer.pal(8, "Set1"),
                 RColorBrewer::brewer.pal(8, "Dark2"))[1:length(top_paths)]

# Step 8a: Compute segment midpoint
top_pathway_long <- top_pathway_long %>%
  dplyr::arrange(Sample, desc(description)) %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(midpoint = cumsum(Abundance) - Abundance / 2) %>%
  dplyr::ungroup()

# Step 8b: Total abundance per sample
bar_totals <- top_pathway_long %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(Total = sum(Abundance), .groups = "drop")

# Step 9: Plot and save as SVG
svg("Top 10 Predicted MetaCyc Pathways Across Samples.svg", width = 12, height = 8)
ggplot2::ggplot(top_pathway_long, aes(x = Sample, y = Abundance, fill = description)) +
  geom_bar(stat = "identity") +
  # Segment-level labels
  geom_text(aes(y = midpoint, label = round(Abundance, 1)),
            size = 3, fontface = "bold", color = "black") +
  # Total per bar
  geom_text(data = bar_totals,
            aes(x = Sample, y = Total, label = round(Total, 1)),
            vjust = -0.4, size = 3, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = path_colors) +
  theme_minimal() +
  labs(title = "Top 10 Predicted MetaCyc Pathways Across Samples",
       x = "Sample",
       y = "Predicted Abundance",
       fill = "MetaCyc Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
dev.off()





# Step 1: Load pathway contribution data
pathway_contrib_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_contrib.tsv.gz"
path_contrib <- read_tsv(pathway_contrib_file)

# Step 2: Load pathway descriptions
pathway_annot_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/metacyc-pwy_name.txt.gz"
path_annot <- read_tsv(pathway_annot_file, col_names = FALSE)
colnames(path_annot) <- c("pathway_id", "description")

# Step 3: Define pathway of interest
target_pathway <- "PWY0-1586"  # Replace with your chosen MetaCyc ID

# Step 4: Filter and annotate
path_target <- path_contrib %>%
  filter(`function` == target_pathway)

# Get pathway description
path_desc <- path_annot %>%
  filter(pathway_id == target_pathway) %>%
  pull(description)

# Step 5: Summarize contribution by taxon
path_contrib_summary <- path_target %>%
  group_by(taxon) %>%
  summarise(Total_Contribution = sum(taxon_rel_function_abun, na.rm = TRUE)) %>%
  arrange(desc(Total_Contribution)) %>%
  slice_head(n = 15)

# Step 6: Save barplot with annotations
filename <- paste0("Top Taxa Contributing to ", target_pathway, " - ", path_desc, ".svg")

svg(filename, width = 10, height = 6)
ggplot(path_contrib_summary, aes(x = reorder(taxon, Total_Contribution),
                                 y = Total_Contribution,
                                 fill = taxon)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(Total_Contribution, 2)),
            hjust = -0.2, fontface = "bold", size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))[1:nrow(path_contrib_summary)]) +
  theme_minimal() +
  labs(title = paste("Top Taxa Contributing to", target_pathway, "-", path_desc),
       x = "Taxon",
       y = "Relative Contribution Score") +
  theme(plot.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()







# ----------- Part A: Pathway Richness per Sample -----------

# Step 1: Load unstratified pathway abundance data
pathway_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz"
path_data <- read_tsv(pathway_file)

# Step 2: Count non-zero pathways per sample
pathway_richness <- path_data %>%
  select(-pathway) %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Num_Pathways")

# Step 3: Annotated barplot (SVG export)
svg("Number of Predicted MetaCyc Pathways per Sample.svg", width = 10, height = 6)
ggplot(pathway_richness, aes(x = Sample, y = Num_Pathways, fill = Sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Num_Pathways), vjust = -0.2, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c(brewer.pal(8, "Accent"), brewer.pal(8, "Set2"))[1:nrow(pathway_richness)]) +
  labs(title = "Number of Predicted MetaCyc Pathways per Sample",
       x = "Sample",
       y = "Pathway Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 13, face = "bold"))
dev.off()


# ----------- Part B: NSTI Scores per Sample (Proxy from EC) -----------

# Step 4: Load NSTI scores
nsti_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/weighted_nsti.tsv.gz"
nsti_data <- read_tsv(nsti_file)

# Step 5: Annotated NSTI barplot (SVG export)
svg("NSTI Score per Sample (MetaCyc Inference Confidence).svg", width = 10, height = 6)
ggplot(nsti_data, aes(x = sample, y = weighted_NSTI, fill = sample)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(weighted_NSTI, 3)), vjust = -0.2, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set3"))[1:nrow(nsti_data)]) +
  labs(title = "NSTI Score per Sample (Inference Confidence)",
       x = "Sample",
       y = "Mean NSTI Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 13, face = "bold"))
dev.off()





# Step 1: Load pathway abundance data
pathway_file <- "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz"
pathway_data <- read_tsv(pathway_file)

# Step 2: Convert to matrix (pathways as rows)
path_mat <- pathway_data %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Step 3: Transpose (samples as rows)
path_mat_t <- t(path_mat)

# Step 4: Normalize to relative abundance
path_rel <- path_mat_t / rowSums(path_mat_t)

# Step 5: Compute Bray-Curtis distance
bray_path <- vegdist(path_rel, method = "bray")

# Step 6: Perform PCoA
pcoa_path <- cmdscale(bray_path, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_path$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(pcoa_df)

# Step 7: Save annotated PCoA plot
svg("PCoA of Predicted MetaCyc Pathway Profiles.svg", width = 10, height = 6)
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Sample, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel(size = 4, fontface = "bold") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "PCoA of Predicted MetaCyc Pathway Profiles",
       x = "PCoA Axis 1",
       y = "PCoA Axis 2") +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 10))
dev.off()




# Function to compute absolute + relative abundance per sample for PICRUSt2 results
per_sample_function_abundance <- function(input_file, 
                                          annot_file = NULL, 
                                          output_file, 
                                          function_type = "KO") {
  # Read input data
  data <- readr::read_tsv(input_file)
  colnames(data)[1] <- "function_id"
  
  # Clean function IDs
  data$function_id <- gsub("^.*:", "", data$function_id)
  
  # Add annotations if provided
  if (!is.null(annot_file)) {
    annot <- readr::read_tsv(annot_file, col_names = FALSE)
    colnames(annot) <- c("function_id", "description")
    annot$function_id <- gsub("^.*:", "", annot$function_id)
    data <- dplyr::left_join(data, annot, by = "function_id")
  } else {
    data$description <- data$function_id
  }
  
  # Extract sample columns
  sample_cols <- setdiff(colnames(data), c("function_id", "description"))
  
  # Calculate absolute abundance
  abs_df <- data
  
  # Calculate relative abundance (%)
  rel_df <- data
  rel_df[, sample_cols] <- sweep(rel_df[, sample_cols], 2, 
                                 colSums(rel_df[, sample_cols]), 
                                 FUN = "/") * 100
  
  # Rename relative columns
  rel_cols <- paste0(sample_cols, " (%)")
  colnames(rel_df)[colnames(rel_df) %in% sample_cols] <- rel_cols
  
  # Combine absolute and relative values
  final_df <- dplyr::left_join(abs_df, rel_df, 
                               by = c("function_id", "description"))
  
  # Add total abundance column for sorting
  final_df$Total_absolute <- rowSums(final_df[, sample_cols], na.rm = TRUE)
  
  # Sort by total abundance
  final_df <- dplyr::arrange(final_df, desc(Total_absolute))
  
  # Export to CSV
  readr::write_csv(final_df, output_file)
}

# Apply to KO predictions
per_sample_function_abundance(
  input_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
  annot_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ko_name.txt.gz",
  output_file = "per_sample_KO_abundance.csv",
  function_type = "KO"
)

# Apply to EC predictions
per_sample_function_abundance(
  input_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
  annot_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/ec_name.txt.gz",
  output_file = "per_sample_EC_abundance.csv",
  function_type = "EC"
)

# Apply to Pathway predictions
per_sample_function_abundance(
  input_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/NEW/humphrey/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz",
  annot_file = "~/Desktop/maize_data/TACC/humphrey_picrust2_results/picrust2-master/picrust2/default_files/description_mapfiles/metacyc-pwy_name.txt.gz",
  output_file = "per_sample_Pathway_abundance.csv",
  function_type = "Pathway"
)

