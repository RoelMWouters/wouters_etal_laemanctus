# ==============================================================================
# WHOLE-GENOME PAIRWISE BARCODE PLOT
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")

# Number of chromosomes to plot
TOP_N_CHROMOSOMES <- 36
BIN_SIZE <- 500000 

# Define Raw Samples and Automatically Generate Combinations
RAW_SAMPLES <- c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08")
pair_matrix <- combn(RAW_SAMPLES, 2)

# Define the renaming dictionary
RENAME_MAP <- c(
  "S01" = "PopA.1", "S02" = "PopA.2", "S03" = "PopA.3",
  "S04" = "PopB.1", "S05" = "PopB.2", "S06" = "PopB.3",
  "S07" = "PopC.1", "S08" = "PopC.2"
)

# Function to categorize comparisons (Modify based on your biology)
classify_group <- function(sample1, sample2) {
  case_when(
    sample1 %in% c("PopA.1", "PopA.2", "PopA.3") & sample2 %in% c("PopA.1", "PopA.2", "PopA.3") ~ "Intra-Lineage A",
    sample1 %in% c("PopB.1", "PopB.2", "PopB.3", "PopC.1", "PopC.2") & sample2 %in% c("PopB.1", "PopB.2", "PopB.3", "PopC.1", "PopC.2") ~ "Intra-Lineage B/C",
    TRUE ~ "Inter-Lineage"
  )
}

# ==========================================
# 2. SMART DATA LOADER
# ==========================================
load_pair <- function(s1, s2) {
  name1 <- paste0(s1, "_vs_", s2)
  name2 <- paste0(s2, "_vs_", s1)
  
  if (file.exists(paste0(name1, "_coords.txt"))) {
    df <- read.table(paste0(name1, "_coords.txt"), col.names = c("CHROM", "POS"))
    df$Sample1 <- s1; df$Sample2 <- s2
    return(df)
  } else if (file.exists(paste0(name2, "_coords.txt"))) {
    df <- read.table(paste0(name2, "_coords.txt"), col.names = c("CHROM", "POS"))
    df$Sample1 <- s2; df$Sample2 <- s1
    return(df)
  } else {
    warning(paste("Neither", name1, "nor", name2, "coords.txt found. Skipping."))
    return(NULL)
  }
}

data_list <- list()
for (i in 1:ncol(pair_matrix)) {
  data_list[[i]] <- load_pair(pair_matrix[1, i], pair_matrix[2, i])
}
all_data <- bind_rows(data_list)

# Apply formatting
all_data <- all_data %>%
  mutate(
    Sample1 = unname(RENAME_MAP[Sample1]),
    Sample2 = unname(RENAME_MAP[Sample2]),
    Comparison = paste(Sample1, Sample2, sep = "_vs_"),
    Group = classify_group(Sample1, Sample2)
  )

# ==========================================
# 3. GENOMIC MATH & BINNING
# ==========================================
contig_map <- all_data %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(POS)) %>%
  arrange(desc(chr_len)) %>%
  slice_head(n = TOP_N_CHROMOSOMES) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len)

binned_data <- all_data %>%
  filter(CHROM %in% contig_map$CHROM) %>%
  left_join(contig_map, by = "CHROM") %>%
  mutate(
    BP_cum = POS + tot,
    bin = floor(BP_cum / BIN_SIZE)
  ) %>%
  group_by(Comparison, Group, CHROM, bin) %>%
  summarise(snp_count = n(), BP_cum_center = mean(BP_cum), .groups = 'drop') %>%
  complete(nesting(Comparison, Group), nesting(CHROM, bin), fill = list(snp_count = 0)) %>%
  left_join(contig_map %>% select(CHROM, tot), by = "CHROM") %>%
  mutate(BP_cum_center = (bin * BIN_SIZE) + (BIN_SIZE / 2))

vlines <- c(contig_map$tot, max(contig_map$tot + contig_map$chr_len)) / 1e6

# ==========================================
# 4. THE CHROMOSOMAL BARCODE PLOT
# ==========================================
barcode_plot <- ggplot(binned_data, aes(x = BP_cum_center / 1e6, y = Comparison, fill = snp_count)) +
  geom_tile(width = BIN_SIZE / 1e6, height = 0.8) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = vlines, color = "black", alpha = 0.3, linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("white", "#FFD700", "#B22222", "#000080"), 
    values = rescale(c(0, 50, 500, max(binned_data$snp_count, na.rm = TRUE))),
    name = sprintf("SNPs per\n%skb", BIN_SIZE / 1000),
    na.value = "white"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15), expand = c(0,0)) +
  theme_minimal() +
  labs(x = "Cumulative Genomic Position (Mb)", y = "") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 270),
    strip.background = element_rect(fill = "gray85", color = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "gray95", color = NA),
    legend.position = "right"
  )

ggsave("Lineage_Barcode_Master.png", plot = barcode_plot, width = 18, height = 12, dpi = 300, bg = "white")
