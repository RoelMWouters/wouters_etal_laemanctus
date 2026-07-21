# ==============================================================================
# COMBINED MASH & PI_HAT HEATMAP
# ==============================================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(scales)

# Set your working directory
setwd("C:/Users/roelm/OneDrive/Documenten/Charles University - PhD/Genomics/Laemanctus/strict")

# ==============================================================================
# 1. LOAD AND PARSE MASH MATRIX (COLORS)
# ==============================================================================
print("Parsing MASH distance matrix...")

num_samples <- as.numeric(readLines("ljul_kmer_matrix.txt", n = 1))
col_names <- paste0("V", 1:(num_samples + 1))
raw_mash <- read.table("ljul_kmer_matrix.txt", skip = 1, fill = TRUE, 
                       stringsAsFactors = FALSE, header = FALSE, 
                       col.names = col_names)

long_paths <- raw_mash$V1
samples <- sapply(strsplit(long_paths, "/"), function(x) x[length(x)-1])

mash_dist <- matrix(0, nrow = num_samples, ncol = num_samples, dimnames = list(samples, samples))

for (i in 2:num_samples) {
  distances <- as.numeric(na.omit(as.numeric(raw_mash[i, 2:(i)])))
  mash_dist[i, 1:length(distances)] <- distances
}

# Mirror the matrix
mash_dist[upper.tri(mash_dist)] <- t(mash_dist)[upper.tri(mash_dist)]

# Melt into long format
mash_long <- reshape2::melt(mash_dist)
colnames(mash_long) <- c("Sample1", "Sample2", "Mash_Dist")

# ==============================================================================
# 2. LOAD AND PARSE PI_HAT MATRIX (TEXT VALUES)
# ==============================================================================
print("Parsing PI_HAT matrix...")

genome_data <- read.table("ljul_IBD.genome", header=TRUE, stringsAsFactors=FALSE)
pairs <- genome_data %>% select(IID1, IID2, PI_HAT)

# Mirror the data
pairs_mirrored <- data.frame(IID1 = pairs$IID2, IID2 = pairs$IID1, PI_HAT = pairs$PI_HAT)
pihat_long <- rbind(pairs, pairs_mirrored)

# Add self-comparisons (1.0 identity)
self_pairs <- data.frame(IID1 = samples, IID2 = samples, PI_HAT = 1.0)
pihat_long <- rbind(pihat_long, self_pairs)

# ==============================================================================
# 3. MERGE DATA, FILTER, AND APPLY CUSTOM NAMING
# ==============================================================================
print("Merging and formatting datasets...")

# Join the PI_HAT values onto the MASH framework
combined_data <- mash_long %>%
  left_join(pihat_long, by = c("Sample1" = "IID1", "Sample2" = "IID2")) %>%
  mutate(PI_HAT_Label = ifelse(is.na(PI_HAT), "-", sprintf("%.2f", PI_HAT)))

# DEFINE EXACT CUSTOM ORDER (C09 and C10 omitted):
custom_order <- c("C01", "C02", "C03", "C04", "C05", "C06", "F01", "F03")

# Filter out any samples not in our target list (drops C09 and C10)
combined_data <- combined_data %>%
  filter(Sample1 %in% custom_order & Sample2 %in% custom_order)

# Dictionary mapping old IDs to new descriptive names
name_map <- c(
  "C01" = "UNAH.1",
  "C02" = "UNAH.2",
  "C03" = "UNAH.3",
  "C04" = "Suyapa.1",
  "C05" = "Suyapa.2",
  "C06" = "Suyapa.3",
  "F01" = "Zoo.1",
  "F03" = "Zoo.2"
)

# Apply the new descriptive names to the dataframe
combined_data$Sample1 <- name_map[as.character(combined_data$Sample1)]
combined_data$Sample2 <- name_map[as.character(combined_data$Sample2)]

# Map the custom order to the new long names
long_order <- name_map[custom_order]

# STRICT ORDERING: X-axis normal (L to R), Y-axis reversed (Top to Bottom)
combined_data$Sample2 <- factor(combined_data$Sample2, levels = long_order)
combined_data$Sample1 <- factor(combined_data$Sample1, levels = rev(long_order))

# ==============================================================================
# 4. GENERATE COMBINED FIGURE
# ==============================================================================
print("Generating Combined Heatmap...")

plot_combined <- ggplot(combined_data, aes(x = Sample2, y = Sample1, fill = Mash_Dist)) +
  geom_tile(color = "white", size = 1) +
  
  # Text maps to the PI_HAT values
  geom_text(aes(label = PI_HAT_Label), color = "black", size = 4, fontface = "bold") +
  
  # Color maps to the MASH distances
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "navy", 
                       midpoint = 0.0035,           
                       limits = c(0, 0.007),        
                       oob = scales::squish,       
                       name = "Mash Distance") +
  
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(title = "", 
       subtitle = "",
       x = "", y = "") +
  
  # Adjusted angle=45 and hjust=0 so the long names tilt upwards and don't overlap
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 11, face = "bold"), 
        axis.text.y = element_text(size = 11, face = "bold"),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        legend.position = "right") +
  coord_fixed()

# Save the figure (Increased dimensions slightly to fit the long labels)
ggsave("Figure_Combined_Mash_PiHat.png", plot = plot_combined, width = 11, height = 9, dpi = 300, bg = "white")

print("Combined figure successfully generated!")
