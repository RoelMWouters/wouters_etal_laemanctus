# ==============================================================================
# PHYLOGENOMIC TREE VISUALIZATION PIPELINE
# Generates Maximum Likelihood & eGCF Concordance Trees
# 
# Description: This script processes IQ-TREE output, performs exact binomial 
# tests on gene concordance factors (gCF), and plots the topologies. It includes
# a custom engine to override ggtree's default ladderization.
# ==============================================================================

library(ggtree)
library(treeio)
library(tidyverse)
library(ape)

# Set working directory to the location of your tree and stat files
# setwd("path/to/your/working/directory") 

# Shared styling constants for the plots
node_blue <- "#1f78b4"
brlen_red <- "#8B0000"

# ==============================================================================
# Taxonomic Translation
# Standardizes tip labels from raw alignment IDs to formal publication names.
# ==============================================================================
translate_taxa <- function(labels) {
  labels <- gsub("^Outgroup$", "B.vittatus_GCA_046270045.1", labels)
  labels <- gsub("^Parent1_Ref$", "L.serratus_GCA_041295975.1", labels)
  labels <- gsub("^C10$", "L.serratus_zoo", labels)
  labels <- gsub("^C09$", "L.longipes_zoo", labels)
  
  # Standardize clonal lineage IDs
  labels <- gsub("^C01_", "Laemanctus_UNAH.1_", labels)
  labels <- gsub("^C02_", "Laemanctus_UNAH.2_", labels)
  labels <- gsub("^C03_", "Laemanctus_UNAH.3_", labels)
  labels <- gsub("^C04_", "Laemanctus_Suyapa.1_", labels)
  labels <- gsub("^C05_", "Laemanctus_Suyapa.2_", labels)
  labels <- gsub("^C06_", "Laemanctus_Suyapa.3_", labels)
  labels <- gsub("^F01_", "Laemanctus_Zoo.1_", labels)
  labels <- gsub("^F03_", "Laemanctus_Zoo.3_", labels)
  
  # Subgenome translation
  labels <- gsub("_Mat$", "_lon", labels)
  labels <- gsub("_Pat$", "_ser", labels)
  
  return(labels)
}

# ==============================================================================
# Deterministic FigTree Layout Engine
# ggtree automatically reorders/ladderizes branches based on clade size, which
# can scramble nested relationships. This vector explicitly defines the exact visual top-to-bottom
# order required.
# ==============================================================================
figtree_top_to_bottom <- c(
  # PATERNAL CLADE (Top to Bottom)
  "Parent1_Ref",
  "C10",
  "C04_Pat",
  "C06_Pat",
  "C05_Pat",
  "C01_Pat",
  "C02_Pat",
  "C03_Pat",
  "F03_Pat",
  "F01_Pat",
  
  # MATERNAL CLADE (Top to Bottom)
  "C09",
  "C03_Mat",
  "C02_Mat",
  "C01_Mat",
  "C04_Mat",
  "C06_Mat",
  "C05_Mat",
  "F01_Mat",
  "F03_Mat",
  
  # OUTGROUP (Absolute Bottom)
  "Outgroup"
)

# ggtree draws from y=1 (bottom) upwards, so we reverse the vector for internal mapping
ggtree_bottom_to_top <- rev(figtree_top_to_bottom)

# This function mathematically locks the tree$edge matrix so ggtree's internal
# C-level postorder sorting cannot override the custom visual sequence above.
lock_figtree_order <- function(tree, bottom_to_top_vec) {
  Ntip <- length(tree$tip.label)
  
  present_target <- bottom_to_top_vec[bottom_to_top_vec %in% tree$tip.label]
  missing_tips   <- setdiff(tree$tip.label, present_target)
  full_target    <- c(present_target, missing_tips)
  
  old_to_new_tip <- match(tree$tip.label, full_target)
  tree$tip.label <- full_target
  
  is_tip_col2 <- tree$edge[, 2] <= Ntip
  tree$edge[is_tip_col2, 2] <- old_to_new_tip[tree$edge[is_tip_col2, 2]]
  
  get_min_rank <- function(node) {
    if (node <= Ntip) {
      return(node)
    } else {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      return(min(sapply(children, get_min_rank)))
    }
  }
  
  new_edges <- integer(0)
  traverse <- function(node) {
    edge_indices <- which(tree$edge[, 1] == node)
    if (length(edge_indices) > 0) {
      children <- tree$edge[edge_indices, 2]
      ranks <- sapply(children, get_min_rank)
      sorted_idx <- edge_indices[order(ranks)]
      for (idx in sorted_idx) {
        new_edges <<- c(new_edges, idx)
        child <- tree$edge[idx, 2]
        if (child > Ntip) {
          traverse(child)
        }
      }
    }
  }
  
  root_node <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))[1]
  traverse(root_node)
  
  tree$edge <- tree$edge[new_edges, , drop = FALSE]
  if (!is.null(tree$edge.length)) {
    tree$edge.length <- tree$edge.length[new_edges]
  }
  attr(tree, "order") <- "cladewise"
  return(tree)
}

# ==============================================================================
# 1. LOAD DATA, COMPUTE BINOMIAL CONCORDANCE & EXPORT SUMMARY CSV
# ==============================================================================
tree_file <- "clean_19tip_concordance.cf.tree"
stat_file <- "clean_19tip_concordance.cf.stat"

my_tree  <- read.tree(tree_file)
cf_stats <- read.table(stat_file, header = TRUE, comment.char = "#", na.strings = "NA", fill = TRUE)

# Calculate Resolved N, Effective gCF (eGCF), and Exact Binomial p-value (H0: p <= 1/3)
cf_stats <- cf_stats %>%
  rowwise() %>%
  mutate(
    resolved_N = gCF_N + gDF1_N + gDF2_N,
    against_N = gDF1_N + gDF2_N,
    eGCF = ifelse(resolved_N > 0, round((gCF_N / resolved_N) * 100, 1), NA_real_),
    binom_p = ifelse(resolved_N > 0, 
                     binom.test(gCF_N, resolved_N, p = 1/3, alternative = "greater")$p.value, 
                     1.0)
  ) %>%
  ungroup()

rooted_tree <- root(my_tree, outgroup = "Outgroup", resolve.root = TRUE)
rooted_tree <- lock_figtree_order(rooted_tree, ggtree_bottom_to_top)

Ntip <- length(rooted_tree$tip.label)
node_map <- rep(NA_integer_, length(rooted_tree$node.label))
clade_desc_vec <- rep("Collapsed Intra-Clonal / Polytomy Noise (Binomial p >= 0.001)", nrow(cf_stats))

# Map nodes to stats iteratively using branch length to prevent collisions
for (i in seq_along(rooted_tree$node.label)) {
  lbl <- rooted_tree$node.label[i]
  node_num <- Ntip + i
  
  if (grepl("/", lbl)) {
    parts <- as.numeric(unlist(strsplit(lbl, "/")))
    gcf_val <- parts[2]
    scf_val <- parts[3]
    
    edge_idx <- which(rooted_tree$edge[, 2] == node_num)
    br_len <- rooted_tree$edge.length[edge_idx]
    
    row_idx <- which.min(abs(cf_stats$gCF - gcf_val) + abs(cf_stats$sCF - scf_val) + abs(cf_stats$Length - br_len)*100)
    node_map[i] <- row_idx
    
    desc_tips <- extract.clade(rooted_tree, node_num)$tip.label
    
    # Auto-generate biological descriptions based on clade tip composition
    desc_str <- case_when(
      setequal(desc_tips, c("Parent1_Ref", "C10", "C01_Pat", "C02_Pat", "C03_Pat", "C04_Pat", "C05_Pat", "C06_Pat", "F01_Pat", "F03_Pat")) ~ "Paternal Subgenome Root (L. serratus Ref + C10 + all ser clones)",
      setequal(desc_tips, c("C10", "C01_Pat", "C02_Pat", "C03_Pat", "C04_Pat", "C05_Pat", "C06_Pat", "F01_Pat", "F03_Pat")) ~ "Paternal Zoo + Hybrids Clade (C10 + all ser clones)",
      setequal(desc_tips, c("C01_Pat", "C02_Pat", "C03_Pat", "C04_Pat", "C05_Pat", "C06_Pat", "F01_Pat", "F03_Pat")) ~ "SINGLE-ORIGIN PATERNAL CLADE: All Hybrid ser Clones",
      setequal(desc_tips, c("C01_Pat", "C02_Pat", "C03_Pat")) ~ "Paternal UNAH Clonal Lineage (UNAH 1-3 ser)",
      setequal(desc_tips, c("C04_Pat", "C05_Pat", "C06_Pat", "F01_Pat", "F03_Pat")) ~ "Paternal Suyapa/Zoo Clonal Lineage (Suyapa/Zoo ser)",
      setequal(desc_tips, c("C09", "C01_Mat", "C02_Mat", "C03_Mat", "C04_Mat", "C05_Mat", "C06_Mat", "F01_Mat", "F03_Mat")) ~ "Maternal Subgenome Root (C09 L. longipes + all lon clones)",
      setequal(desc_tips, c("C01_Mat", "C02_Mat", "C03_Mat", "C04_Mat", "C05_Mat", "C06_Mat", "F01_Mat", "F03_Mat")) ~ "SINGLE-ORIGIN MATERNAL CLADE: All Hybrid lon Clones",
      setequal(desc_tips, c("C09", "C01_Mat", "C02_Mat", "C03_Mat")) ~ "DUAL-ORIGIN SPLIT: Sexual L. longipes (C09) + UNAH_lon Clade",
      setequal(desc_tips, c("C01_Mat", "C02_Mat", "C03_Mat")) ~ "Maternal UNAH Clonal Lineage (UNAH 1-3 lon)",
      setequal(desc_tips, c("C04_Mat", "C05_Mat", "C06_Mat", "F01_Mat", "F03_Mat")) ~ "Maternal Suyapa/Zoo Clonal Lineage (Suyapa/Zoo lon)",
      TRUE ~ paste0("Sub-clade (", paste(translate_taxa(desc_tips), collapse = " + "), ")")
    )
    clade_desc_vec[row_idx] <- desc_str
  }
}

cf_stats <- cf_stats %>%
  mutate(
    Clade_Description = clade_desc_vec,
    Status = ifelse(binom_p < 0.001, "SIGNIFICANT (Kept)", "NOISE (Collapsed)")
  ) %>%
  select(ID, Clade_Description, Status, gCF_N, gDF1_N, gDF2_N, gDFP_N, against_N, resolved_N, eGCF, binom_p, sCF, Length, everything())

write.csv(cf_stats, "Laemanctus_clean_concordance_statistics_summary.csv", row.names = FALSE)

# ==============================================================================
# FIGURE I: MAXIMUM LIKELIHOOD TREE (ALL LABELS SHOWN)
# ==============================================================================
ml_tree <- rooted_tree

for (i in seq_along(ml_tree$node.label)) {
  lbl <- ml_tree$node.label[i]
  node_num <- Ntip + i
  edge_idx <- which(ml_tree$edge[, 2] == node_num)
  
  if (grepl("/", lbl) && length(edge_idx) == 1 && !is.na(node_map[i])) {
    row_data <- cf_stats[node_map[i], ]
    # Format: eGCF% (Votes For / Votes Against)
    ml_tree$node.label[i] <- paste0(row_data$eGCF, "% (", row_data$gCF_N, "/", row_data$against_N, ")")
  } else {
    ml_tree$node.label[i] <- ""
  }
}

ml_tree$tip.label <- translate_taxa(ml_tree$tip.label)

outgroup_edge_ml <- which(ml_tree$edge[, 2] == which(ml_tree$tip.label == "B.vittatus_GCA_046270045.1"))
true_outgroup_len_ml <- round(ml_tree$edge.length[outgroup_edge_ml], 4)
ml_tree$edge.length[outgroup_edge_ml] <- 0.020

# CRITICAL: ladderize = FALSE prevents ggtree from silently reordering branches
p_ml <- ggtree(ml_tree, linewidth = 0.7, ladderize = FALSE) 

p_ml$data <- p_ml$data %>%
  mutate(
    support_text = ifelse(!isTip & label != "", label, NA),
    brlen_text = case_when(
      label == "B.vittatus_GCA_046270045.1" ~ as.character(true_outgroup_len_ml),
      TRUE ~ NA_character_ 
    ),
    slash_mark = ifelse(label == "B.vittatus_GCA_046270045.1", "//", NA)
  )

p_fig1_ml <- p_ml +
  geom_tiplab(size = 4.2, fontface = "italic", offset = 0.0015) +
  geom_nodepoint(data = function(d) filter(d, !is.na(support_text)), color = node_blue, size = 2.8) +
  geom_text(aes(x = branch, y = y, label = support_text), color = node_blue, size = 3.2, fontface = "bold", vjust = -1.2, na.rm = TRUE) +
  geom_text(aes(x = branch, y = y, label = brlen_text), color = brlen_red, size = 3.2, vjust = 1.8, na.rm = TRUE) +
  geom_text(aes(x = branch, y = y, label = slash_mark), color = "black", size = 6, fontface = "bold", vjust = 0.4, na.rm = TRUE) +
  theme_tree2() + hexpand(0.35) + vexpand(0.05) + coord_cartesian(clip = "off") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0), plot.margin = margin(t = 20, r = 40, b = 20, l = 15))

print(p_fig1_ml)
ggsave("Laemanctus_Fig1_Clean_ML_Tree.jpg", plot = p_fig1_ml, width = 16, height = 12, dpi = 300)
ggsave("Laemanctus_Fig1_Clean_ML_Tree.pdf", plot = p_fig1_ml, width = 16, height = 12)


# ==============================================================================
# FIGURE ii: COLLAPSED eGCF CONCORDANCE TREE (eGCF > 50%)
# ==============================================================================
egcf_tree <- rooted_tree

for (i in seq_along(egcf_tree$node.label)) {
  lbl <- egcf_tree$node.label[i]
  node_num <- Ntip + i
  edge_idx <- which(egcf_tree$edge[, 2] == node_num)
  
  if (grepl("/", lbl) && length(edge_idx) == 1 && !is.na(node_map[i])) {
    row_data <- cf_stats[node_map[i], ]
    
    # Collapse if eGCF is missing or strictly below 50
    if (is.na(row_data$eGCF) || row_data$eGCF < 51) {
      egcf_tree$edge.length[edge_idx] <- 0
      egcf_tree$node.label[i] <- "BLANK" 
    } else {
      egcf_tree$node.label[i] <- paste0("ROW_", node_map[i])
    }
  } else {
    egcf_tree$node.label[i] <- "BLANK" 
  }
}

collapsed_tree <- di2multi(egcf_tree, tol = 1e-6)
collapsed_tree <- read.tree(text = write.tree(collapsed_tree))
collapsed_tree <- lock_figtree_order(collapsed_tree, ggtree_bottom_to_top)

for (i in seq_along(collapsed_tree$node.label)) {
  lbl <- collapsed_tree$node.label[i]
  if (grepl("^ROW_", lbl)) {
    r_idx <- as.numeric(sub("^ROW_", "", lbl))
    row_data <- cf_stats[r_idx, ]
    collapsed_tree$node.label[i] <- paste0(row_data$eGCF, "% (", row_data$gCF_N, "/", row_data$against_N, ")")
  } else {
    collapsed_tree$node.label[i] <- ""
  }
}

collapsed_tree$tip.label <- translate_taxa(collapsed_tree$tip.label)

outgroup_edge <- which(collapsed_tree$edge[, 2] == which(collapsed_tree$tip.label == "B.vittatus_GCA_046270045.1"))
true_outgroup_len <- round(collapsed_tree$edge.length[outgroup_edge], 4)
collapsed_tree$edge.length[outgroup_edge] <- 0.020

p_egcf <- ggtree(collapsed_tree, linewidth = 0.7, ladderize = FALSE) 

p_egcf$data <- p_egcf$data %>%
  mutate(
    support_text = ifelse(!isTip & label != "", label, NA),
    brlen_text = case_when(
      label == "B.vittatus_GCA_046270045.1" ~ as.character(true_outgroup_len),
      TRUE ~ NA_character_ 
    ),
    slash_mark = ifelse(label == "B.vittatus_GCA_046270045.1", "//", NA)
  )

p_fig2_egcf <- p_egcf +
  geom_tiplab(size = 4.2, fontface = "italic", offset = 0.0015) +
  geom_nodepoint(data = function(d) filter(d, !is.na(support_text)), color = node_blue, size = 2.8) +
  # Lifted blue support text higher (vjust = -1.2)
  geom_text(aes(x = branch, y = y, label = support_text), color = node_blue, size = 3.2, fontface = "bold", vjust = -1.2, na.rm = TRUE) +
  # Dropped red branch length text lower (vjust = 1.8)
  geom_text(aes(x = branch, y = y, label = brlen_text), color = brlen_red, size = 3.2, vjust = 1.8, na.rm = TRUE) +
  geom_text(aes(x = branch, y = y, label = slash_mark), color = "black", size = 6, fontface = "bold", vjust = 0.4, na.rm = TRUE) +
  theme_tree2() + hexpand(0.35) + vexpand(0.05) + coord_cartesian(clip = "off") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0), plot.margin = margin(t = 20, r = 40, b = 20, l = 15))

print(p_fig2_egcf)
# Stretched canvas size to 16x12
ggsave("Laemanctus_Fig2_Clean_Collapsed_eGCF_50.jpg", plot = p_fig2_egcf, width = 16, height = 12, dpi = 300)
ggsave("Laemanctus_Fig2_Clean_Collapsed_eGCF_50.pdf", plot = p_fig2_egcf, width = 16, height = 12)

