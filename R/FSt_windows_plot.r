library(ggplot2)
library(dplyr)
library(readr)
setwd("/path/to/local/working/directory")

fst_AB <- read_tsv("fst_A_vs_B.windowed.weir.fst") %>% mutate(Comparison = "Pop A vs Pop B", Color = "navy")
fst_AC <- read_tsv("fst_A_vs_C.windowed.weir.fst") %>% mutate(Comparison = "Pop A vs Pop C", Color = "darkorange")

fst_combined <- bind_rows(fst_AB, fst_AC) %>%
  rename(Chrom = CHROM, Start = BIN_START, End = BIN_END, Fst = WEIGHTED_FST) %>%
  mutate(Fst = ifelse(Fst < 0, 0, Fst))

fst_combined$Comparison <- factor(fst_combined$Comparison, levels = c("Pop A vs Pop B", "Pop A vs Pop C"))

density_plot <- ggplot(fst_combined, aes(x = Fst, fill = Comparison, color = Comparison)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values = c("Pop A vs Pop B" = "navy", "Pop A vs Pop C" = "darkorange")) +
  scale_color_manual(values = c("Pop A vs Pop B" = "navy", "Pop A vs Pop C" = "darkorange")) +
  theme_minimal() +
  labs(title = "Distribution of Genomic Divergence", x = expression(F[ST]), y = "Number of 100kb Windows") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("Figure_5_Fst_Density.png", plot = density_plot, width = 8, height = 6, dpi = 300)
