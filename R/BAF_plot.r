library(tidyverse)
library(data.table)
library(ggpointdensity)
setwd("/path/to/local/working/directory")

MIN_DP <- 15   
MAX_DP <- 30 
TARGET_DENSITY_PER_MB <- 60 

df <- fread("BAF_Data.tsv", fill = TRUE, colClasses = "character")
fai <- fread("genome.fasta.fai")
colnames(fai) <- c("CHROM", "Full_Length", "Offset_fai", "Linebases", "Linewidth")

active_chroms <- intersect(df$CHROM, fai$CHROM)
fai_active <- fai %>% filter(CHROM %in% active_chroms) %>%
  arrange(desc(as.numeric(Full_Length))) %>% head(36) %>%
  mutate(tot = cumsum(as.numeric(Full_Length)) - as.numeric(Full_Length))

axis_df <- fai_active %>% mutate(center = tot + (as.numeric(Full_Length) / 2))

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

plot_data_wide <- df %>% filter(CHROM %in% fai_active$CHROM) %>%
  inner_join(fai_active, by = "CHROM") %>%
  mutate(POS = as.numeric(POS), BP_cum = POS + tot) %>%
  mutate(SampleA = calculate_baf_vector(S01_AD, MIN_DP, MAX_DP))

plot_data <- plot_data_wide %>%
  select(CHROM, Full_Length, BP_cum, SampleA) %>%
  pivot_longer(cols = c("SampleA"), names_to = "Sample", values_to = "BAF") %>%
  drop_na(BAF)

set.seed(42)
plot_data_thinned <- plot_data %>% group_by(Sample, CHROM) %>%
  mutate(target_n = round((as.numeric(Full_Length) / 1e6) * TARGET_DENSITY_PER_MB)) %>%
  sample_n(size = min(n(), first(target_n))) %>% ungroup() %>%
  mutate(BAF_jittered = BAF + runif(n(), min = -0.015, max = 0.015))

ggplot(plot_data_thinned, aes(x = BP_cum / 1e6, y = BAF_jittered)) +
  geom_pointdensity(size = 0.4, alpha = 0.8) +
  facet_grid(Sample ~ .) +
  scale_color_gradientn(colors = c("#e0f3db", "#a8ddb5", "#43a2ca", "#0868ac", "#02213b")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#d35400", linewidth = 0.5) +
  scale_x_continuous(label = fai_active$CHROM, breaks = axis_df$center / 1e6) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.02, 1.02)) +
  theme_bw() + 
  labs(title = "Whole-Genome BAF Landscape", x = "Cumulative Genomic Position (Mb)", y = "B-Allele Frequency") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, hjust = 1))

ggsave("BAF_Manhattan_DensityScatter.png", width = 24, height = 10, dpi = 400)
