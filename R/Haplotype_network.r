# Load libraries
library(ape)
library(pegas)

# Set generalized working directory
setwd("/path/to/local/working/directory")

# ==============================================================================
# 1. LOAD DATA, FILTER, & EXTRACT HAPLOTYPES
# ==============================================================================
mito_dna <- read.dna("Mitochondrial_Supermatrix.fasta", format = "fasta")

# Safely drop duplicate/unwanted outgroups to keep only the primary representative
mito_dna <- mito_dna[rownames(mito_dna) != "outgroup_drop", ]

haps <- haplotype(mito_dna)

# Calculate true pairwise SNP distances
d_haps <- dist.dna(haps, model = "N")

# Generate the base network
net_mito <- haploNet(haps, d = d_haps)

# Clean up the network by removing the infinite alternative dotted links
attr(net_mito, "alter.links") <- NULL

# ==============================================================================
# 2. METADATA MAPPING (Restoring Sample Names)
# ==============================================================================
tip_names <- rownames(mito_dna)
group_labels <- rep("Unknown", length(tip_names))

# Map specific individuals to their taxonomic/geographic groups
group_labels[grep("S10", tip_names)] <- "L. serratus" 
group_labels[grep("S09", tip_names)] <- "L. longipes"
group_labels[grep("S0[1-3]", tip_names)] <- "L. julioi (Pop B)"
group_labels[grep("S0[4-6]", tip_names)] <- "L. julioi (Pop A)"
group_labels[grep("S0[7-8]", tip_names)] <- "L. julioi (Captive)"

# Convert to a factor
group_factor <- as.factor(group_labels)

# Build the contingency table for the pie charts
hap_frequencies <- haploFreq(mito_dna, fac = group_factor, haplo = haps)

# ==============================================================================
# 3. THE LONG BRANCH CAP 
# ==============================================================================
# Cap massive mutational steps at 15/50. The plot will naturally draw hatch 
# marks instead of trying to squeeze hundreds of lines into the space.
if("step" %in% colnames(net_mito)) {
  net_mito[, "step"] <- ifelse(net_mito[, "step"] > 50, 50, net_mito[, "step"])
} else {
  net_mito[, 3] <- ifelse(net_mito[, 3] > 50, 50, net_mito[, 3])
}

# ==============================================================================
# 4. PLOT AND SAVE THE FIGURE
# ==============================================================================
# Colors are applied alphabetically based on the factor levels
my_colors <- c("skyblue", "gold", "mediumseagreen", "darkorange", "purple")

png("Mito_Network_Final.png", width = 10, height = 8, units = "in", res = 300)

par(mar = c(1, 1, 1, 6), xpd = TRUE)

plot(net_mito, 
     size = attr(net_mito, "freq"),   # Node size scales with sample size
     pie = hap_frequencies,           # Colors defined by taxonomic groups
     bg = my_colors,                  
     scale.ratio = 0.2,               # Spaces out the central network
     cex = 1.2,                       # Makes the hatch marks larger and clearer
     labels = FALSE,                  
     show.mutation = 1)               # 1 = Perpendicular Hatch Marks

# Optional: Float text next to the long branch (adjust x and y as needed)
# text(x = 0, y = 5, labels = "// ~400 SNPs", col = "red", font = 2)

legend("right", inset = c(-0.15, 0), 
       legend = levels(group_factor), 
       fill = my_colors, 
       bty = "n", 
       cex = 1.1,
       title = "Samples")

dev.off()
