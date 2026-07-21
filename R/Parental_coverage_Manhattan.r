# ==============================================================================
# SUBGENOME PHASING DIAGNOSTIC (BIPOLAR MANHATTAN PLOT)
# ==============================================================================
library(tidyverse)

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")
INPUT_FILE <- "sliding_window_coverage.tsv"

# Plotting parameters
TOP_N_CHROMOSOMES <- 36
PARENT_A_NAME <- "Parent A"
PARENT_B_NAME <- "Parent B"

# Colors: -1.0 = Parent B, 0 = Shared, +1.0 = Parent A
BIAS_PALETTE <- c("#12E0AC", "grey80", "#78D134") 

# ==========================================
# 2. LOAD & CALCULATE BIPOLAR BIAS
# ==========================================
cov_data <- read_tsv(INPUT_FILE, show_col_types = FALSE)

plot_data <- cov_data %>%
  filter(Reads_Parent_A + Reads_Parent_B > 0) %>%
  mutate(
    Total_Reads = Reads_Parent_A + Reads_Parent_B,
    Bipolar_Bias = (Reads_Parent_A - Reads_Parent_B) / Total_Reads,
    Pos_Mb = (Start + End) / 2 / 1e6
  )

# ==========================================
# 3. CHROMOSOME LAYOUT & OFFSETS
# ==========================================
contig_lengths <- plot_data %>%
  group_by(Contig) %>%
  summarize(Length_Mb = max(Pos_Mb)) %>%
  arrange(desc(Length_Mb)) %>%
  head(TOP_N_CHROMOSOMES)

fai_active <- contig_lengths %>%
  mutate(
    tot = cumsum(Length_Mb) - Length_Mb,
    center = tot + (Length_Mb / 2)
  )

plot_data_manhattan <- plot_data %>%
  filter(Contig %in% fai_active$Contig) %>%
  left_join(fai_active %>% select(Contig, tot), by = "Contig") %>%
  mutate(
    BP_cum = Pos_Mb + tot,
    Color_Group = ifelse(as.numeric(factor(Contig, levels = fai_active$Contig)) %% 2 == 0, "Even", "Odd")
  )

# ==========================================
# 4. GENERATE THE MANHATTAN PLOT
# ==========================================
manhattan_plot <- ggplot(plot_data_manhattan, aes(x = BP_cum, y = Bipolar_Bias)) +
  geom_point(aes(color = Bipolar_Bias), size = 0.8, alpha = 0.7) +
  scale_color_gradientn(
    colors = BIAS_PALETTE, 
    limits = c(-1, 1),
    name = sprintf("Subgenome Dominance\n(-1 = %s | +1 = %s)", PARENT_B_NAME, PARENT_A_NAME)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.6, 0.6), linetype = "dotted", color = "darkred", linewidth = 0.5) +
  scale_x_continuous(
    label = fai_active$Contig, 
    breaks = fai_active$center, 
    expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(breaks = c(-1.0, -0.6, 0.0, 0.6, 1.0), limits = c(-1.05, 1.05)) +
  theme_bw() +
  labs(
    title = sprintf("Subgenome Phasing Landscape (Top %s Contigs)", TOP_N_CHROMOSOMES),
    subtitle = "Diagnostic Manhattan Plot for Chimeras and Conserved Chromosomes",
    x = "Cumulative Genomic Position (Mb)",
    y = "Divergence Score (Bipolar Bias)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("Phasing_Diagnostic_Manhattan_Bipolar.png", plot = manhattan_plot, width = 20, height = 6, dpi = 300)
