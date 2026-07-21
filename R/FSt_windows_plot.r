# ==============================================================================
# FST DIVERGENCE DENSITY PLOT
# ==============================================================================
library(ggplot2)
library(dplyr)
library(readr)
library(scales) 

# ==========================================
# 1. USER CONTROLS & I/O
# ==========================================
setwd("/path/to/working/directory")

# Map your input files to their respective display labels
comparisons <- list(
  "fst_A_vs_B.windowed.weir.fst" = "Pop A vs Pop B",
  "fst_A_vs_C.windowed.weir.fst" = "Pop A vs Pop C",
  "fst_B_vs_C.windowed.weir.fst" = "Pop B vs Pop C"
)

# Colors for each comparison
comp_colors <- c(
  "Pop A vs Pop B" = "navy", 
  "Pop A vs Pop C" = "darkorange",
  "Pop B vs Pop C" = "lightseagreen"
)

# ==========================================
# 2. LOAD AND CLEAN DATA
# ==========================================
data_list <- list()

for (file in names(comparisons)) {
  if (file.exists(file)) {
    df <- read_tsv(file, show_col_types = FALSE) %>%
      mutate(Comparison = comparisons[[file]])
    data_list[[file]] <- df
  }
}

fst_combined <- bind_rows(data_list) %>%
  rename(Chrom = CHROM, Start = BIN_START, End = BIN_END, Fst = WEIGHTED_FST) %>%
  mutate(
    Fst = ifelse(Fst < 0, 0, Fst),
    Midpoint = (Start + End) / 2
  )

fst_combined$Comparison <- factor(fst_combined$Comparison, levels = unname(unlist(comparisons)))

# ==========================================
# 3. GENERATE DENSITY PLOT
# ==========================================
density_plot <- ggplot(fst_combined, aes(x = Fst, fill = Comparison, color = Comparison)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values = comp_colors) +
  scale_color_manual(values = comp_colors) +
  scale_y_continuous(labels = label_comma()) +
  theme_minimal() +
  labs(
    x = expression(F[ST]), 
    y = "Number of Genomic Windows"
  ) + 
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black")
  ) + 
  guides(fill = guide_legend(nrow = 3), color = "none")

ggsave("Fst_Density_Distribution.png", plot = density_plot, width = 8, height = 6, dpi = 300)
