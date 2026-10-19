# ==============================================================================
# RAxML MITOCHONDRIAL TOPOLOGY PLOTTING
# ==============================================================================
# Load libraries
library(ggtree)
library(treeio)
library(ape)
library(ggplot2)

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")

INPUT_TREE <- "RAxML_bipartitions.mt_markers"
OUTPUT_IMAGE <- "RAxML_Mito_Topology_Final.png"
PLOT_TITLE <- "Mitochondrial Topology with Bootstrap Support"
OUTGROUP_ID <- "basvit"

# Minimum branch length required to draw a length label
MIN_BRANCH_LEN_LABEL <- 0.00035

# Dictionary mapping raw tip IDs to final publication-ready names
TAXA_MAP <- c(
  "C10" = "L.serratus_zoo",
  "C09" = "L.longipes_zoo",
  "C01" = "Laemanctus_UNAH.1",
  "C02" = "Laemanctus_UNAH.2",
  "C03" = "Laemanctus_UNAH.3",
  "C04" = "Laemanctus_Suyapa.1",
  "C05" = "Laemanctus_Suyapa.2",
  "C06" = "Laemanctus_Suyapa.3",
  "F01" = "Laemanctus_Zoo.1",
  "F03" = "Laemanctus_Zoo.3",
  "basvit" = "B.vittatus_GCA_046270045.1",
  "lser" = "L.serratus_GCA_041295975.1"
)

# ==============================================================================
# 2. READ, ROOT, AND TRANSLATE TAXA
# ==============================================================================
my_tree <- read.tree(INPUT_TREE)
rooted_tree <- root(my_tree, outgroup = OUTGROUP_ID, resolve.root = TRUE)

for (raw_name in names(TAXA_MAP)) {
  rooted_tree$tip.label <- gsub(raw_name, TAXA_MAP[[raw_name]], rooted_tree$tip.label)
}

FINAL_OUTGROUP <- TAXA_MAP[[OUTGROUP_ID]]
if(is.null(FINAL_OUTGROUP)) FINAL_OUTGROUP <- OUTGROUP_ID

# ==============================================================================
# 3. SMART BRANCH SHRINKING
# ==============================================================================
root_node <- Ntip(rooted_tree) + 1
basal_edges <- which(rooted_tree$edge[,1] == root_node)
basal_children <- rooted_tree$edge[basal_edges, 2] 

outgroup_idx <- which(rooted_tree$tip.label == FINAL_OUTGROUP)
outgroup_edge <- which(rooted_tree$edge[,2] == outgroup_idx)
true_outgroup_length <- round(rooted_tree$edge.length[outgroup_edge], 6)

rooted_tree$edge.length[basal_edges] <- 0.05 
nodes_to_hide <- c(basal_children, outgroup_idx)

# ==============================================================================
# 4. PLOT AND SAVE
# ==============================================================================
p <- ggtree(rooted_tree) + 
  geom_tiplab(size=3.5, offset=0.0002, fontface="italic") + 
  
  # Display Bootstrap values (Accessing the node labels directly from the tree object)
  geom_nodelab(aes(label=ifelse(label == "Root", "", label)), 
               size=3, vjust=-0.5, color="blue") +
  
  # Bulletproof branch lengths
  geom_text(aes(x = branch, 
                label = ifelse(!is.na(branch.length) & branch.length >= MIN_BRANCH_LEN_LABEL & !(node %in% nodes_to_hide), 
                               round(branch.length, 6), "")), 
            vjust = -0.5, size = 2.5, color = "darkred") + 
  
  # Outgroup annotation
  geom_text2(aes(subset = (node == outgroup_idx), x = branch), 
             label = "//", size = 5, vjust = 0.3, color = "black", fontface = "bold") +
  geom_text2(aes(subset = (node == outgroup_idx), x = branch), 
             label = true_outgroup_length, size = 2.5, vjust = -1.5, color = "darkred") +
  
  theme_tree2() +
  hexpand(0.4) +
  ggtitle(PLOT_TITLE)

ggsave(OUTPUT_IMAGE, plot = p, width = 10, height = 8, dpi = 300)
