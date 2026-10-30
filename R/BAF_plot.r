# ==============================================================================
# B-ALLELE FREQUENCY (BAF) MANHATTAN PLOT
# ==============================================================================
library(tidyverse)
library(data.table)
library(ggpointdensity)

# ==========================================
# 1. USER CONTROLS
# ==========================================
setwd("/path/to/working/directory")

INPUT_TSV <- "BAF_Data.tsv"
INPUT_FAI <- "Genome.fasta.fai"

# Depth and Plotting parameters
MIN_DP <- 15   
MAX_DP <- 30 
TARGET_DENSITY_PER_MB <- 60 
TOP_N_CHROMOSOMES <- 36

# Filters
CONTIGS_TO_REMOVE <- c("contig_x", "contig_y")

# Optional: Truncate a specific contig to remove artifactual chimeras for example
# Set TRUNCATE_CONTIG to NULL to disable
TRUNCATE_CONTIG <- "contig_name" 
TRUNCATE_CUTOFF <- 278731350 #coordinates

# Sample Mapping (Update column headers to match your TSV)
SAMPLE_COLS <- c("PopA" = "S01_AD", "PopB" = "S04_AD", "PopC" = "S07_AD")
# ==========================================

# 2. Load Data & Index
df <- fread(INPUT_TSV, fill = TRUE, colClasses = "character")
fai <- fread(INPUT_FAI)
colnames(fai) <- c("CHROM", "Full_Length", "Offset_fai", "Linebases", "Linewidth")

if (!is.null(TRUNCATE_CONTIG)) {
  fai <- fai %>%
    mutate(Full_Length = ifelse(CHROM == TRUNCATE_CONTIG, TRUNCATE_CUTOFF, as.numeric(Full_Length)))
}

# 3. Synchronize, Filter, & Sort
active_chroms <- intersect(df$CHROM, fai$CHROM)
fai_active <- fai %>% 
  filter(CHROM %in% active_chroms) %>%
  filter(!CHROM %in% CONTIGS_TO_REMOVE) %>%   
  arrange(desc(as.numeric(Full_Length))) %>% 
  head(TOP_N_CHROMOSOMES) %>%
  mutate(tot = cumsum(as.numeric(Full_Length)) - as.numeric(Full_Length)) %>%
  select(CHROM, tot, Full_Length)

# 4. Strict Depth BAF Calculation
calculate_baf_vector <- function(ad_col, min_dp, max_dp) {
  temp <- strsplit(as.character(ad_col), ",", fixed = TRUE)
  sapply(temp, function(x) {
    if(length(x) < 2) return(NA)
    r <- as.numeric(x[1]); a <- as.numeric(x[2])
    total_dp <- r + a
    if(is.na(total_dp) | total_dp < min_dp | total_dp > max_dp) return(NA)
    return(a / total_dp)
  })
}

# 5. Transform and Stitch
plot_data_wide <- df %>%
  filter(CHROM %in% fai_active$CHROM)

if (!is.null(TRUNCATE_CONTIG)) {
  plot_data_wide <- plot_data_wide %>%
    filter(!(CHROM == TRUNCATE_CONTIG & as.numeric(POS) > TRUNCATE_CUTOFF))
}

plot_data_wide <- plot_data_wide %>%
  inner_join(fai_active, by = "CHROM") %>%
  mutate(POS = as.numeric(POS), BP_cum = POS + tot)

# Compute BAF for each requested sample dynamically
for (sample_name in names(SAMPLE_COLS)) {
  col_name <- SAMPLE_COLS[[sample_name]]
  plot_data_wide[[sample_name]] <- calculate_baf_vector(plot_data_wide[[col_name]], MIN_DP, MAX_DP)
}

# 6. Pivot to Long Format
plot_data <- plot_data_wide %>%
  select(CHROM, Full_Length, BP_cum, all_of(names(SAMPLE_COLS))) %>%
  pivot_longer(cols = all_of(names(SAMPLE_COLS)), names_to = "Sample", values_to = "BAF") %>%
  drop_na(BAF) %>%
  mutate(
    Sample = factor(Sample, levels = names(SAMPLE_COLS)),
    color_group = ifelse(as.numeric(factor(CHROM, levels = fai_active$CHROM)) %% 2 == 0, "A", "B")
  )

# 7. Proportional Downsampling & Pre-Jittering
set.seed(42)
plot_data_thinned <- plot_data %>%
  group_by(Sample, CHROM) %>%
  mutate(target_n = round((as.numeric(Full_Length) / 1e6) * TARGET_DENSITY_PER_MB)) %>%
  sample_n(size = min(n(), first(target_n))) %>% 
  ungroup() %>%
  mutate(BAF_jittered = BAF + runif(n(), min = -0.015, max = 0.015))

# 8. Generate Plot
ggplot(plot_data_thinned, aes(x = BP_cum / 1e6, y = BAF_jittered)) +
  geom_pointdensity(size = 0.4, alpha = 0.8) +
  facet_grid(Sample ~ .) +
  scale_color_gradientn(
    colors = c("#e0f3db", "#a8ddb5", "#43a2ca", "#0868ac", "#02213b"),
    name = "Local\nDensity"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#d35400", linewidth = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.02, 1.02)) +
  theme_bw() + 
  labs(x = "Cumulative Genomic Position (Mb)", y = "B-Allele Frequency") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    strip.text.y = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("BAF_Manhattan_DensityScatter.png", width = 24, height = 10, dpi = 400)
