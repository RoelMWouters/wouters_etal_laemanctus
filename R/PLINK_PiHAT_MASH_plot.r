# ==============================================================================
# COMBINED MASH & PI_HAT HEATMAP
# ==============================================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(scales)

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")

# Input files
MASH_FILE <- "kmer_matrix.txt"
PIHAT_FILE <- "IBD.genome"

# Define the exact samples to include (drops any unlisted samples like outgroups)
TARGET_SAMPLES <- c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08")

# Dictionary mapping raw IDs to final publication-ready names
NAME_MAP <- c(
  "S01" = "PopA.1",
  "S02" = "PopA.2",
  "S03" = "PopA.3",
  "S04" = "PopB.1",
  "S05" = "PopB.2",
  "S06" = "PopB.3",
  "S07" = "PopC.1",
  "S08" = "PopC.2"
)

# ==========================================
# 2. LOAD AND PARSE MASH MATRIX (COLORS)
# ==========================================
num_samples <- as.numeric(readLines(MASH_FILE, n = 1))
col_names <- paste0("V", 1:(num_samples + 1))
raw_mash <- read.table(MASH_FILE, skip = 1, fill = TRUE, 
                       stringsAsFactors = FALSE, header = FALSE, 
                       col.names = col_names)

long_paths <- raw_mash$V1
samples <- sapply(strsplit(long_paths, "/"), function(x) x[length(x)-1])

mash_dist <- matrix(0, nrow = num_samples, ncol = num_samples, dimnames = list(samples, samples))

for (i in 2:num_samples) {
  distances <- as.numeric(na.omit(as.numeric(raw_mash[i, 2:(i)])))
  mash_dist[i, 1:length(distances)] <- distances
}

mash_dist[upper.tri(mash_dist)] <- t(mash_dist)[upper.tri(mash_dist)]
mash_long <- reshape2::melt(mash_dist)
colnames(mash_long) <- c("Sample1", "Sample2", "Mash_Dist")

# ==========================================
# 3. LOAD AND PARSE PI_HAT MATRIX (TEXT)
# ==========================================
genome_data <- read.table(PIHAT_FILE, header=TRUE, stringsAsFactors=FALSE)
pairs <- genome_data %>% select(IID1, IID2, PI_HAT)

pairs_mirrored <- data.frame(IID1 = pairs$IID2, IID2 = pairs$IID1, PI_HAT = pairs$PI_HAT)
pihat_long <- rbind(pairs, pairs_mirrored)

self_pairs <- data.frame(IID1 = samples, IID2 = samples, PI_HAT = 1.0)
pihat_long <- rbind(pihat_long, self_pairs)

# ==========================================
# 4. MERGE DATA, FILTER, & APPLY NAMING
# ==========================================
combined_data <- mash_long %>%
  left_join(pihat_long, by = c("Sample1" = "IID1", "Sample2" = "IID2")) %>%
  mutate(PI_HAT_Label = ifelse(is.na(PI_HAT), "-", sprintf("%.2f", PI_HAT))) %>%
  filter(Sample1 %in% TARGET_SAMPLES & Sample2 %in% TARGET_SAMPLES)

combined_data$Sample1 <- NAME_MAP[as.character(combined_data$Sample1)]
combined_data$Sample2 <- NAME_MAP[as.character(combined_data$Sample2)]

long_order <- NAME_MAP[TARGET_SAMPLES]

combined_data$Sample2 <- factor(combined_data$Sample2, levels = long_order)
combined_data$Sample1 <- factor(combined_data$Sample1, levels = rev(long_order))

# ==========================================
# 5. GENERATE COMBINED FIGURE
# ==========================================
plot_combined <- ggplot(combined_data, aes(x = Sample2, y = Sample1, fill = Mash_Dist)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = PI_HAT_Label), color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "navy", 
                       midpoint = 0.0035, limits = c(0, 0.007), 
                       oob = scales::squish, name = "Mash Distance") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 11, face = "bold"), 
        axis.text.y = element_text(size = 11, face = "bold"),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        legend.position = "right") +
  coord_fixed()

ggsave("Combined_Mash_PiHat_Heatmap.png", plot = plot_combined, width = 11, height = 9, dpi = 300, bg = "white")
