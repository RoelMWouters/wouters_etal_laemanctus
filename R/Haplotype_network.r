# ==============================================================================
# MITOCHONDRIAL HAPLOTYPE NETWORK
# ==============================================================================
# Load libraries
library(ape)
library(pegas)

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")

INPUT_FASTA <- "Mitochondrial_Supermatrix.fasta"
OUTPUT_IMAGE <- "Mito_Network_Final.png"

# Specific raw samples to drop entirely (e.g., secondary outgroups or failed samples)
SAMPLES_TO_DROP <- c("lser_outgroup")

# Map regex patterns (found in FASTA headers) to their final grouping name
GROUP_MAP <- list(
  "C10"     = "serratus",
  "C09"     = "longipes",
  "C0[1-3]" = "UNAH",
  "C0[4-6]" = "Suyapa",
  "F0[1,3]" = "Zoo"
)

# Lock the order of the factor levels for consistent color mapping
GROUP_ORDER <- c("UNAH", "Suyapa", "Zoo", "longipes", "serratus")

# Colors apply strictly in the order of GROUP_ORDER
GROUP_COLORS <- c("#3366CC", "#FF9933", "#E63946", "#78D134", "#12E0AC")

# Maximum mutational steps to draw before capping (prevents long-branch clutter)
MAX_STEPS <- 50 

# ==============================================================================
# 2. LOAD DATA, FILTER, & EXTRACT HAPLOTYPES
# ==============================================================================
mito_dna <- read.dna(INPUT_FASTA, format = "fasta")

# Safely drop unwanted samples
mito_dna <- mito_dna[!(rownames(mito_dna) %in% SAMPLES_TO_DROP), ]

haps <- haplotype(mito_dna)
d_haps <- dist.dna(haps, model = "N")
net_mito <- haploNet(haps, d = d_haps)

# Clean up the network by removing the infinite alternative dotted links
attr(net_mito, "alter.links") <- NULL

# ==============================================================================
# 3. METADATA MAPPING (Restoring Sample Names & Forcing Order)
# ==============================================================================
tip_names <- rownames(mito_dna)
group_labels <- rep("Unknown", length(tip_names))

for (pattern in names(GROUP_MAP)) {
  group_labels[grep(pattern, tip_names)] <- GROUP_MAP[[pattern]]
}

group_factor <- factor(group_labels, levels = GROUP_ORDER)
hap_frequencies <- haploFreq(mito_dna, fac = group_factor, haplo = haps)

# ==============================================================================
# 4. THE LONG BRANCH CAP 
# ==============================================================================
if("step" %in% colnames(net_mito)) {
  net_mito[, "step"] <- ifelse(net_mito[, "step"] > MAX_STEPS, MAX_STEPS, net_mito[, "step"])
} else {
  net_mito[, 3] <- ifelse(net_mito[, 3] > MAX_STEPS, MAX_STEPS, net_mito[, 3])
}

# ==============================================================================
# 5. PLOT AND SAVE THE FIGURE
# ==============================================================================
png(OUTPUT_IMAGE, width = 10, height = 8, units = "in", res = 300)
par(mar = c(1, 1, 1, 6), xpd = TRUE)

plot(net_mito, 
     size = attr(net_mito, "freq"),   
     pie = hap_frequencies,           
     bg = GROUP_COLORS,                  
     scale.ratio = 0.2,               
     cex = 1.2,                       
     labels = FALSE,                  
     show.mutation = 1)               

legend("right", inset = c(-0.15, 0), 
       legend = levels(group_factor), 
       fill = GROUP_COLORS, 
       bty = "n", 
       cex = 1.1,
       title = "Samples")

dev.off()
