# ==============================================================================
# PAIRWISE SNP DISTANCE & DIVERGENCE HEATMAPS
# Generates Heatmaps for Clonal Lineage Comparisons
# 
# Description: This script processes PLINK distance matrices (absolute SNP 
# differences and percentage divergence/1-IBS). It standardizes sample names, 
# enforces a strict manual ordering to group clonal lineages visually, and 
# plots color-blind friendly heatmaps without hierarchical clustering to 
# preserve the biological grouping.
# ==============================================================================

# Install required packages if not already installed
# install.packages(c("pheatmap", "RColorBrewer"))

library(pheatmap)
library(RColorBrewer)

# Set working directory to the location of your distance matrix files
# setwd("path/to/your/working/directory") 

# ==============================================================================
# 1. LOAD AND FORMAT DATA
# ==============================================================================
# Read the sample IDs. 
# PLINK writes '.id' files as 'FID IID'; we extract the second column (IID).
sample_ids <- read.table("ljul_clones_haploid_mapped_absolute_diffs.dist.id", stringsAsFactors = FALSE)[,2]

# Load the Absolute SNP Differences matrix
dist_mat <- as.matrix(read.table("ljul_clones_haploid_mapped_absolute_diffs.dist"))
rownames(dist_mat) <- sample_ids
colnames(dist_mat) <- sample_ids

# Load the Percentage Divergence matrix (1 - Identity By State)
div_mat <- as.matrix(read.table("ljul_clones_haploid_mapped_percent_divergence.mdist"))
rownames(div_mat) <- sample_ids
colnames(div_mat) <- sample_ids

# ==============================================================================
# 2. RENAME AND REORDER SAMPLES FOR PUBLICATION
# ==============================================================================
# Create a dictionary to map the raw sequencing/bioinformatics IDs to the 
# formal nomenclature used in the publication.
rename_map <- c(
  "C01" = "UNAH1",
  "C02" = "UNAH2",
  "C03" = "UNAH3",
  "C04" = "Suyapa_mother",
  "C05" = "Suyapa_offspring1",
  "C06" = "Suyapa_offspring2",
  "F01" = "Zoo_mother",
  "F03" = "Zoo_offspring"
)

# Apply the formal publication names to the row and column names of both matrices
rownames(dist_mat) <- rename_map[rownames(dist_mat)]
colnames(dist_mat) <- rename_map[colnames(dist_mat)]
rownames(div_mat) <- rename_map[rownames(div_mat)]
colnames(div_mat) <- rename_map[colnames(div_mat)]

# Define the exact visual sequence for the axes.
# This manually groups the lineages (UNAH, then Suyapa, then Zoo) together.
desired_order <- c("UNAH1", "UNAH2", "UNAH3", 
                   "Suyapa1", "Suyapa_offspring1", "Suyapa_offspring2", 
                   "Zoo_mother", "Zoo_offspring")

# Subset and reorder both matrices to match the desired biological order
dist_mat <- dist_mat[desired_order, desired_order]
div_mat <- div_mat[desired_order, desired_order]

# ==============================================================================
# 3. AESTHETIC SETTINGS FOR PUBLICATION
# ==============================================================================
# Create a color-blind friendly, high-contrast palette. 
# Dark Blue = 0 distance (Identical), Yellow = Medium Divergence, Dark Red = High Divergence
pub_palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)

# Global heatmap graphical parameters
cell_border <- "grey80"  # Subtle grid lines to define individual matrix cells
font_base_size <- 12
number_size <- 10

# ==============================================================================
# 4. GENERATE FIGURE 2A: ABSOLUTE SNP DIFFERENCES
# ==============================================================================
# Open a high-resolution PDF graphics device
pdf("Fig2A_Absolute_SNP_Differences.pdf", width = 8, height = 7)

pheatmap(
  dist_mat,
  color = pub_palette,
  border_color = cell_border,
  display_numbers = TRUE,
  number_format = "%.0f",            # Uses whole numbers (integers) for exact SNP counts
  fontsize_number = number_size,
  fontsize = font_base_size,
  angle_col = 90,                    # Rotates X-axis labels 90 degrees to prevent overlapping
  cluster_rows = FALSE,              # CRITICAL: Disables auto-clustering to enforce desired_order
  cluster_cols = FALSE,              # CRITICAL: Disables auto-clustering to enforce desired_order
  main = "Pairwise Absolute SNP Differences"
)

dev.off() # Close the graphics device and save the PDF

# ==============================================================================
# 5. GENERATE FIGURE 2B: PERCENTAGE DIVERGENCE
# ==============================================================================
pdf("Fig2B_Percentage_Divergence.pdf", width = 8, height = 7)

pheatmap(
  div_mat,
  color = pub_palette,
  border_color = cell_border,
  display_numbers = TRUE,
  number_format = "%.3f",            # Uses 3 decimal places (e.g., 0.001) for fractional divergence
  fontsize_number = number_size,
  fontsize = font_base_size,
  angle_col = 90,
  cluster_rows = FALSE,              # CRITICAL: Disables auto-clustering
  cluster_cols = FALSE,              # CRITICAL: Disables auto-clustering
  main = "Pairwise Sequence Divergence (1 - IBS)"
)

dev.off()
