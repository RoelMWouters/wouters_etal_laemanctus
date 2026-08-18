# ==============================================================================
# Plotting Genomic Concordance Factors (gCF) using ggtree
# ==============================================================================
# Description:
# This script reads IQ-TREE concordance factor output files (.tree and .stat),
# applies customizable taxonomic label translations, parses embedded gCF values 
# and decisive gene counts (n), and generates a publication-ready cladogram 
# with uniform branch lengths.
#
# Dependencies:
#   install.packages(c("tidyverse", "ape", "ggplot2"))
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("ggtree", "treeio"))
# ==============================================================================

library(ggtree)
library(treeio)
library(tidyverse)
library(ape)
library(ggplot2)

# ==============================================================================
# 1. READ INPUT FILES AND RENAME TAXA
# ==============================================================================
tree_file <- "subgenome_concordance1000.cf.tree"
stat_file <- "subgenome_concordance1000.cf.stat"

# Load the tree and concordance factor statistics table
tree <- read.tree(tree_file)
cf_stats <- read.table(stat_file, header = TRUE, comment.char = "#", na.strings = "NA", fill = TRUE)

# Preserve the tree structure (assuming root has already been set natively by IQ-TREE)
rooted_tree <- tree 

# Taxonomic translation / renaming 
# (Modify these regex patterns or names to match your specific sample identifiers)
rooted_tree$tip.label <- gsub("C01", "Laemanctus_UNAH.1", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("C02", "Laemanctus_UNAH.2", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("C03", "Laemanctus_UNAH.3", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("C04", "Laemanctus_Suyapa.1", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("C05", "Laemanctus_Suyapa.2", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("C06", "Laemanctus_Suyapa.3", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("F01", "Laemanctus_Zoo.1", rooted_tree$tip.label)
rooted_tree$tip.label <- gsub("F03", "Laemanctus_Zoo.3", rooted_tree$tip.label)

# ==============================================================================
# 2. PARSE gCF STATS AND MERGE WITH TREE DATA
# ==============================================================================
# Initialize ggtree object with uniform branch lengths and a ladderized layout
p <- ggtree(rooted_tree, branch.length = "none", ladderize = TRUE)

tree_data <- p$data

# Extract the gCF value natively embedded in the Newick node label (e.g., "bootstrap/gCF")
tree_data <- tree_data %>%
  mutate(
    gCF_extracted = ifelse(grepl("/", label), as.numeric(sub(".*/", "", label)), NA)
  )

# Match the extracted gCF value with the closest corresponding gN (decisive gene count)
# from the .stat file to reliably bypass internal node ID re-indexing issues.
tree_data$gN <- sapply(tree_data$gCF_extracted, function(val) {
  if (is.na(val)) return(NA)
  closest_row <- which.min(abs(cf_stats$gCF - val))
  return(cf_stats$gN[closest_row])
})

# Format text labels for the plot (e.g., "64.3% (n=297)")
merged_data <- tree_data %>%
  mutate(
    cf_text = ifelse(!is.na(gCF_extracted), paste0(gCF_extracted, "% (n=", gN, ")"), "")
  )

p$data <- merged_data

# ==============================================================================
# 3. GENERATE PUBLICATION-READY CLADOGRAM
# ==============================================================================
p_final <- p + 
  # Tip labels (styled in bold italics with an offset for spacing)
  geom_tiplab(size = 3.8, offset = 0.5, fontface = "bold.italic") + 
  
  # Node points (strictly filtered to omit root and tip nodes lacking gCF metrics)
  geom_nodepoint(data = function(d) filter(d, !is.na(gCF_extracted)), color = "#0072B2", size = 3.2, alpha = 0.85) +
  
  # Concordance factor text labels positioned adjacent to nodes
  geom_text(aes(x = x, y = y, label = cf_text), color = "#0072B2", size = 3, fontface = "bold", vjust = -1.2, hjust = 1.1) +
  
  theme_tree() +
  hexpand(0.6) + # Expand canvas horizontally on the right side for tip labels
  vexpand(0.15) + # Add vertical padding
  coord_cartesian(clip = "off") + # Prevent labels from being clipped outside drawing bounds
  
  labs(
    title = "Genomic Concordance Factor (gCF) Cladogram"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 120) # Generous left margin to accommodate deep node text
  )

# Preview plot in R session
print(p_final)

# Export high-resolution vector graphic for publication or repository archive
ggsave("GenomeConcordance.svg", plot = p_final, width = 10, height = 8, dpi = 300)
