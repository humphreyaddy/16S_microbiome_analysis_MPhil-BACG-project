# This is an R workshop on microbial analysis

# Author: Kwasi

# Loading libraries

library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.

library(plotly) # A package to create interactive web graphics of use in 3D plots

library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.

library(philr) # This package provides functions for the analysis of compositional data 

library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step

library(adespatial) # Tools for the multiscale spatial analysis of multivariate data

library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library(microbiomeSeq) # Data analysis and visualization

library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(ranacapa) # Data analysis 

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(microbiomeutilities) # some utility tools 


#setwd("/home/humphrey/Desktop/maize_data/R_analysis_viz")
setwd("/home/humphrey/Desktop/maize_data/TACC/output/NEW")
# Create a phyloseq object
physeq <- readRDS("../maize_data_23.06.2025.rds")

phy <- qza_to_phyloseq("table.qza","rooted-tree.qza","taxonomy.qza","sample-metadata.tsv")

phy

# Or

# Importing abundance table
ASVs <- read_qza("table.qza")

# Import Metadata
metadata <- read.table("sample-metadata.tsv", sep='\t',header=T, row.names=1, comment="")

metadata <- metadata[-1,] # remove second line that specifies data type

# Import tree
tree <- read_qza("rooted-tree.qza")

# Import taxonomy
taxonomy <- read_qza("taxonomy.qza")
tax_table <- do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), ";"))
colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table) <- taxonomy$data$Feature.ID


# Creating phyloseq object
physeq <- phyloseq(
  otu_table(ASVs$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(tax_table),
  sample_data(metadata)
)

# Summarize the phyloseq  
summarize_phyloseq(physeq)
print(physeq)
summary(sample_sums(physeq))

# Accessing the phyloseq object
ntaxa(physeq)

nsamples(physeq)

sample_names(physeq)[1:6]

rank_names(physeq)

sample_variables(physeq)

otu_table(physeq)[1:5, 1:4]

tax_table(physeq)[1:5, 1:4]


###Saving a phyloseq Object###
#saveRDS(physeq,"maize_data_23.06.2025.rds")

#### Load phyloseq object ####-------------------------------------------------------
#physeq <- readRDS("final_humphrey.rds")
#physeq <- readRDS("maize_data_23.06.2025.rds")


# Plot the distribution of the reads
plot_read_distribution(physeq, groups = "Production_conditions", plot.type = "density") + theme_biome_utils()
plot_read_distribution(physeq, groups = "isolation_source", plot.type = "density") + theme_biome_utils()

# Rarefy phyloseq object
#physeq_rarefy <- rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
#physeq_rarefy <- rarefy_even_depth(physeq, sample.size=67619, rngseed=42,replace=F)

# Taxa prevalence at phylum level
plot_taxa_prevalence(physeq, "Phylum")
plot_taxa_prevalence(physeq, "Family")
plot_taxa_prevalence(physeq, "Genus")

session_info()

# Composition plots

# Barplot

physeq_df <- microbiomeutilities::phy_to_ldf(physeq, 
                                             transform.counts = "compositional")

# An additonal column Sam_rep with sample names is created for reference purpose
colnames(physeq_df)

# Box plot at Family level

ggstripchart(physeq_df, "food.type", "Abundance", 
             facet.by = "Family", color = "food.type",
             palette = "jco") + rremove("x.text")

# Subset to top 10 most abundant families
physeq_fam_10 <- aggregate_top_taxa2(physeq, level = "Family", top = 10)
physeq_fam_10

# Converting to relative abundance
physeq_fam_10_rel <- microbiome::transform(physeq_fam_10, "compositional")

?aggregate_top_taxa2
?microbiome::transform

# plot
plot_composition(physeq_fam_10_rel, sample.sort = "food.type", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

# alternative plot without phyloseq
taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "food.type")
?taxa_barplot

# To make it interactive
ggplotly(taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "food.type"))

# Save the plot
b.plot <- taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "food.type")

ggsave("barplot_family.png", b.plot,  width = 14, height = 10, dpi = 300)

# Heatmap
taxa_heatmap(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "food.type")

library(pheatmap)
# Alternative
p <- plot_taxa_heatmap(physeq,
                       subset.top = 15,
                       VariableA = c("food.type","Production_conditions", "Raw_material"),
                       transformation = "log10",
                       cluster_rows = T,
                       cluster_cols = F,
                       show_colnames = F,
                       heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
)
#the plot is stored here
p$plot

# To display the plot
print(p$plot)

# Optional: save to file
ggsave("Top15_Taxa_Heatmap.svg", p$plot, width = 12, height = 8, dpi = 300)


# table used for plot is here
p$tax_tab[1:3,1:3]





# Agglomerate to Family level
physeq_fam <- tax_glom(physeq, taxrank = "Family")

# Transform to relative abundance
physeq_fam_10_rel <- transform_sample_counts(physeq_fam, function(x) x / sum(x))

# Filter to top 10 families if needed
top10 <- names(sort(taxa_sums(physeq_fam_10_rel), decreasing = TRUE))[1:10]
physeq_fam_10_rel <- prune_taxa(top10, physeq_fam_10_rel)








# Third option
h.map <- plot_heatmap(physeq_fam_10_rel, method="PCoA", distance="bray", taxa.label = "Family", sample.order = unique(sample_names(physeq))) + facet_grid(~food.type, scales = "free_x", drop = TRUE) + theme_bw() + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)) + theme(legend.key = element_blank(),strip.background = element_rect(colour="black", fill="white"))

# Make bacterial names italics
h.map <- h.map + theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'))

# Change the color palette
h.map <- h.map + scale_fill_distiller("Abundance", palette = "RdYlBu")

# clean the x-axis
h.map <- h.map + rremove("x.text")

print(h.map)

# Saving the plot
ggsave("heatmap_family.svg", h.map,  width = 14, height = 10, dpi = 300)

# Boxplot

physeq_df <- microbiomeutilities::phy_to_ldf(physeq_fam_10, 
                                             transform.counts = "compositional")

# An additonal column Sam_rep with sample names is created for reference purpose
colnames(physeq_df)

# Box plot at Family level

ggstripchart(physeq_df, "food.type", "Abundance", 
             facet.by = "Family", color = "food.type",
             palette = "jco") + rremove("x.text")


# Plot relative abundance
mycols <- c("tomato" ,"#f5f3cc","#958c0d", "#dba11c",  "#00c4ad", "#005456", "#277342", "#2f4858", "#673e00", "#8e005e")

plot_taxa_boxplot(physeq,
                  taxonomic.level = "Phylum",
                  top.otu = 6, 
                  group = "food.type",
                  add.violin= FALSE,
                  title = "Top three family", 
                  keep.other = TRUE,
                  group.order = c("S19_Fermented_maize", "S37_Fermented_Maiz", "S38_D_Maize", "S4_Ogi", "S5_Mawe", "S6_Mawe"),
                  group.colors = mycols,
                  dot.size = 2) + theme_biome_utils()

# Plot taxa/OTU specified by user (which ones gave most enrichment in diversity)
physeq.f <- format_to_besthit(physeq) # top 5 responsible

top_taxa(physeq.f, 5)

#select.taxa <- c("d29fe3c70564fc0f69f2c03e0d1e5561:k__Bacteria", "154709e160e8cada6bfb21115acc80f5: ovatus")

select.taxa <- c("4921d046c2494b5875e699d4a6efc4b6:d__Bacteria", "581a55014641e1cd55cdc272c0365a28:Lactobacillus_fermentum", "c667664df1fb6f8705f9a8a861370a72: Lactobacillus_delbrueckii",
                 "bf1dd9519026f72be520b483b4142308:d__Bacteria", "9e9c4270966cf507485654eb33a69880:d__Bacteria" )

p <- plot_listed_taxa(physeq.f, select.taxa, 
                      group= "food.type",
                      group.order = c("S19_Fermented_maize", "S37_Fermented_Maiz", "S38_D_Maize", "S4_Ogi", "S5_Mawe", "S6_Mawe"),
                      group.colors = mycols,
                      add.violin = F,
                      dot.opacity = 0.25,
                      box.opacity = 0.25,
                      panel.arrange= "grid") + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent)


# Adding statistical test with ggpubr::stat_compare_means()

# If more than two variables/pairwise comparison
comps <- make_pairs(sample_data(physeq.f)$food.type)
print(comps)

p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 

p + scale_y_continuous(labels = scales::percent)

# Plot top 4 genera
physeq.genus <- aggregate_taxa(physeq, "Genus")
top_four <- top_taxa(physeq.genus, 4)
top_four

top_genera <- plot_listed_taxa(physeq.genus, top_four, 
                               group= "food.type",
                               group.order = c("S19_Fermented_maize", "S37_Fermented_Maiz", "S38_D_Maize", "S4_Ogi", "S5_Mawe", "S6_Mawe"),
                               group.colors = mycols,
                               add.violin = F,
                               dot.opacity = 0.25,
                               box.opacity = 0.25,
                               panel.arrange= "wrap")
top_genera

top_genera + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")

?stat_compare_means
?plot_listed_taxa

# Dominant taxa - genus level (can be done at any taxonomic level)
dom.tax <- dominant_taxa(physeq,level = "Genus", group="food.type")
head(dom.tax$dominant_overview)

# Taxa summary
taxa_summary(physeq, "Phylum")

# Group specific abundances of taxa
grp_abundances <- get_group_abundances(physeq,
                                       level = "Phylum",
                                       group="food.type",
                                       transform = "compositional")

# Clean taxa names
grp_abundances$OTUID <- gsub("p__", "",grp_abundances$OTUID)
grp_abundances$OTUID <- ifelse(grp_abundances$OTUID == "", 
                               "Unclassified", grp_abundances$OTUID)

# Plot
mean.plot <- grp_abundances %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=food.type)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("food.type", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

mean.plot

# Samples dominated by specific taxa (not working, debug needed)
physeq.gen <- aggregate_taxa(physeq, "Genus")
bact_dom <- find_samples_taxa(physeq.gen, taxa = "g__Bacteroides")
?find_samples_taxa
#get sample dominated by "g__Bacteroides
ps.sub <- prune_samples(sample_names(physeq.gen) %in% bact_dom, physeq.gen)





# Basic plot with Wilcoxon test
p1 <- ggboxplot(
  physeq_df, 
  x = "food.type", 
  y = "Abundance", 
  fill = "food.type"
) +
  stat_compare_means(comparisons = comps, method = "wilcox.test")  # Add p-values

# Plot with approximate p-values (e.g. for ties or large samples)
p2 <- ggboxplot(
  physeq_df, 
  x = "food.type", 
  y = "Abundance", 
  fill = "food.type"
) +
  stat_compare_means(comparisons = comps, method = "wilcox.test", exact = FALSE)

# Save both plots
ggsave("Boxplot_Exact_Wilcoxon.svg", p1, width = 8, height = 6, dpi = 600)
ggsave("Boxplot_Approx_Wilcoxon.svg", p2, width = 8, height = 6, dpi = 600)

# Save only plotted columns (can adjust as needed)
write.csv(physeq_df[, c("food.type", "Abundance")], 
          "Boxplot_Plot_Values.csv", 
          row.names = FALSE)


# Run stats separately and save
library(dplyr)

# Wilcoxon test for all pairs
wilcox_results <- do.call(rbind, lapply(comps, function(comp) {
  group1 <- comp[1]
  group2 <- comp[2]
  group1_vals <- physeq_df %>% filter(food.type == group1) %>% pull(Abundance)
  group2_vals <- physeq_df %>% filter(food.type == group2) %>% pull(Abundance)
  test <- wilcox.test(group1_vals, group2_vals, exact = FALSE)
  data.frame(Group1 = group1, Group2 = group2, P.value = test$p.value)
}))

# Save to CSV
write.csv(wilcox_results, "Boxplot_Wilcoxon_Stats.csv", row.names = FALSE)





# Alpha diversity
ggrare(physeq, step = 50, color="food.type", label = "Sample", se = TRUE) # check if enough has being sequenced to capture enough ASVs per sample
physeq_rarefy <- ggrare(physeq, step = 50, color="food.type", label = "Sample", se = TRUE) # check if enough has being sequenced to capture enough ASVs per sample

# Alpha diversity plot
plot_richness(physeq_rarefy, measures="Shannon")
plot_richness(physeq, measures="Shannon")

# Alpha diversity by groups
a.div <- plot_richness(physeq_rarefy, x="food.type", measures=c("Shannon", "simpson", "Observed"), color = "food.type") + geom_boxplot() + theme_bw()
a.div <- plot_richness(physeq, x="food.type", measures=c("Shannon", "simpson", "Observed"), color = "food.type") + geom_boxplot() + theme_bw()

a.div

# Adding statistical support
a.div + stat_compare_means(
  comparisons = comps,
  label = "p.signif",
  tip.length = 0.05,
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
    symbols = c("xxxx", "***", "**", "*", "ns")
  ),
  method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90))

# Generate spreadsheet with richness estimates
richness <- estimate_richness(physeq_rarefy,measures=c("Shannon", "simpson", "Observed"))
write.csv(richness, file = "alpha_div.csv")

# Another option for alpha diversity plot
plot_diversity_stats(physeq, group = "food.type", 
                     index = "diversity_shannon", 
                     group.order = c("S19_Fermented_maize", "S37_Fermented_Maiz", "S38_D_Maize", "S4_Ogi", "S5_Mawe", "S6_Mawe"),                      
                     group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE) + ylab("Shannon Diversity") + xlab("")

# Beta diversity

# Non-phylogenetic beta diversity
physeq.ord <- ordinate(physeq_rarefy, "PCoA", "bray")
b.div.bray <- plot_ordination(physeq_rarefy, physeq.ord, type= "samples", color= "food.type") + geom_point(size=3)
b.div.bray <- b.div.bray + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.bray)

# Phylogenetic beta diversity (considers phylogenetic distances)
# convert to relative abundance
physeq_rel <- microbiome::transform(physeq, "compositional")
physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted=T)
b.div.wuni <- plot_ordination(physeq_rel, physeq.ord.wuni, type= "samples", color= "food.type") + geom_point(size=3)
b.div.wuni <- b.div.wuni + stat_ellipse() + ggtitle("Weighted Unifrac")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.wuni)

# Permanova (statistical analysis for beta analysis)

otu <- abundances(physeq_rel)
meta <- meta(physeq_rel)

# Statistics - Bdiv
permanova <- adonis2(t(otu) ~ food.type, data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["food.type", "Pr(>F)"])

# adonis2
permanova <- adonis2(t(otu) ~ food.type, data = meta, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova.txt")


# Assign plot to object
shannon_plot <- plot_diversity_stats(
  physeq, group = "food.type", 
  index = "diversity_shannon", 
  group.order = c("Fante_Kenkey", "Ga_Kenkey"),
  group.colors = mycols,
  label.format = "p.format",
  stats = TRUE
) + ylab("Shannon Diversity") + xlab("") + 
  theme(axis.text.x = element_text(angle = 0))

# Save plot
ggsave("Shannon_Diversity_by_FoodType.svg", shannon_plot, width = 8, height = 6, dpi = 600)


# Extract Shannon values per sample
shannon_vals <- estimate_richness(physeq, measures = "Shannon")

# Add sample metadata (e.g. food type)
shannon_vals$SampleID <- rownames(shannon_vals)
shannon_vals <- merge(shannon_vals, sample_data(physeq), by.x = "SampleID", by.y = "row.names")

# Save to CSV
write.csv(shannon_vals, "Shannon_diversity_values.csv", row.names = FALSE)


ggsave("BetaDiversity_BrayCurtis.svg", b.div.bray, width = 8, height = 6, dpi = 600)



# Extract Bray PCoA coordinates
bray_coords <- as.data.frame(physeq.ord$vectors)
bray_coords$SampleID <- rownames(bray_coords)
bray_coords <- merge(bray_coords, sample_data(physeq_rarefy), by.x = "SampleID", by.y = "row.names")

# Save to CSV
write.csv(bray_coords, "BetaDiversity_BrayCurtis_PCoA.csv", row.names = FALSE)

ggsave("BetaDiversity_WeightedUnifrac.svg", b.div.wuni, width = 8, height = 6, dpi = 600)

# Extract Weighted UniFrac PCoA coordinates
wuni_coords <- as.data.frame(physeq.ord.wuni$vectors)
wuni_coords$SampleID <- rownames(wuni_coords)
wuni_coords <- merge(wuni_coords, sample_data(physeq_rel), by.x = "SampleID", by.y = "row.names")

# Save to CSV
write.csv(wuni_coords, "BetaDiversity_WeightedUnifrac_PCoA.csv", row.names = FALSE)









####--------------------------SAVING CSVs---------------------------------####

# Directory to save CSVs
dir.create("csv_outputs", showWarnings = FALSE)

# 1. Rarefaction curve (14_Rarefraction_on_ggrare.png)
# No direct object assigned, re-run and extract
rarefaction_df <- ggrare(physeq, step = 50, color = "food.type", label = "Sample", se = TRUE)$data
write.csv(rarefaction_df, "csv_outputs/Rarefaction_Curve_Values.csv", row.names = FALSE)

# 2. Alpha diversity measures (16_Alpha_diversity_measures.png)
# Already handled in script
# File: alpha_div.csv and Shannon_diversity_values.csv

# 3. Taxa prevalence at phylum level (2_Taxa_prevalence_at_phylum.png)
# No object saved in script; re-extract:
phylum_prevalence <- plot_taxa_prevalence(physeq, "Phylum")$data
write.csv(phylum_prevalence, "csv_outputs/Taxa_Prevalence_Phylum.csv", row.names = FALSE)

# 4. Taxa abundance at genus level (3_Taxa_abundance_at_genus.png)
genus_abundance <- Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Genus
write.csv(genus_abundance, "csv_outputs/Taxa_Abundance_Genus.csv")

# 5. Top ten family-level abundances (4_Top_Ten_Taxa_Abundance_at_family.png)
top10_fam_abund <- phyloseq::psmelt(physeq_fam_10_rel)
write.csv(top10_fam_abund, "csv_outputs/Top10_Family_Abundance.csv", row.names = FALSE)

# 6. Relative abundance top 15 (5_relative_abundance_to_top_15.png)
# Extract from heatmap plot object
write.csv(p$tax_tab, "csv_outputs/Top15_Taxa_Heatmap_Abundance.csv")

# 7. Relative abundance top 5 (5.1_relative_abundance_to_top_5.png)
# Recreate top 5 table
top5_phylum <- Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Phylum
write.csv(top5_phylum, "csv_outputs/Top5_Taxa_Abundance.csv")

# 8. Family heatmap (9_heatmap_family.png and heatmap_family.png)
# Already handled via: write.csv(p$tax_tab...)
# Confirmed again here just in case:
write.csv(p$tax_tab, "csv_outputs/Heatmap_Family_Abundance.csv")

# 9. Top10 taxa heatmap (Top10_Taxa_Heatmap.png)
# Already same object (p$tax_tab), or if different subset:
# If different, re-define p with top 10 instead of 15
# (Optional re-run if needed)

# 10. Boxplot (Boxplot_Exact_Wilcoxon.png)
write.csv(physeq_df[, c("food.type", "Abundance")], 
          "csv_outputs/Boxplot_Abundance_by_FoodType.csv", 
          row.names = FALSE)

# 11. Beta diversity: Bray-Curtis (BetaDiversity_BrayCurtis.png)
write.csv(bray_coords, "csv_outputs/BetaDiversity_BrayCurtis_PCoA.csv", row.names = FALSE)

# 12. Beta diversity: Weighted UniFrac (BetaDiversity_WeightedUnifrac.png)
write.csv(wuni_coords, "csv_outputs/BetaDiversity_WeightedUnifrac_PCoA.csv", row.names = FALSE)

# 13. Shannon diversity boxplot by food type (Shannon_Diversity_by_FoodType.svg)
write.csv(shannon_vals, "csv_outputs/Shannon_Diversity_BySample.csv", row.names = FALSE)

# 14. Sequencing depth (1_sequencing_depth.png)
# Not assigned but inferred from: sample_sums(physeq)
seq_depth <- data.frame(SampleID = names(sample_sums(physeq)), Reads = sample_sums(physeq))
write.csv(seq_depth, "csv_outputs/Sequencing_Depth_Per_Sample.csv", row.names = FALSE)






# Load saved phyloseq object
physeq <- readRDS("maize_data_23.06.2025.rds")

# Restore ASVs and taxonomy
ASVs <- otu_table(physeq)
tax_table <- tax_table(physeq)

# Re-run missing data exports
genus_abundance <- Summarize.Taxa(ASVs, as.data.frame(tax_table))$Genus
write.csv(genus_abundance, "csv_outputs/Taxa_Abundance_Genus.csv")

top5_phylum <- Summarize.Taxa(ASVs, as.data.frame(tax_table))$Phylum
write.csv(top5_phylum, "csv_outputs/Top5_Taxa_Abundance.csv")

# Shannon diversity
shannon_vals <- estimate_richness(physeq, measures = "Shannon")
shannon_vals$SampleID <- rownames(shannon_vals)
shannon_vals <- merge(shannon_vals, as.data.frame(sample_data(physeq)), by.x = "SampleID", by.y = "row.names")
write.csv(shannon_vals, "csv_outputs/Shannon_Diversity_BySample.csv", row.names = FALSE)

# Rarefy if not yet done
physeq_rarefy <- rarefy_even_depth(physeq, rngseed=42, sample.size=67619, replace=FALSE)

# Beta diversity: Bray-Curtis
physeq.ord <- ordinate(physeq_rarefy, "PCoA", "bray")
bray_coords <- as.data.frame(physeq.ord$vectors)
bray_coords$SampleID <- rownames(bray_coords)
bray_coords <- merge(bray_coords, as.data.frame(sample_data(physeq_rarefy)), by.x = "SampleID", by.y = "row.names")
write.csv(bray_coords, "csv_outputs/BetaDiversity_BrayCurtis_PCoA.csv", row.names = FALSE)

# Beta diversity: Weighted UniFrac
physeq_rel <- microbiome::transform(physeq, "compositional")
physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted = TRUE)
wuni_coords <- as.data.frame(physeq.ord.wuni$vectors)
wuni_coords$SampleID <- rownames(wuni_coords)
wuni_coords <- merge(wuni_coords, as.data.frame(sample_data(physeq_rel)), by.x = "SampleID", by.y = "row.names")
write.csv(wuni_coords, "csv_outputs/BetaDiversity_WeightedUnifrac_PCoA.csv", row.names = FALSE)

