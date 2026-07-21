setwd("/path/to/local/working/directory")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

genome_data <- read.table("IBD.genome", header=TRUE, stringsAsFactors=FALSE)
pairs <- genome_data %>% select(IID1, IID2, PI_HAT)
pairs_mirrored <- data.frame(IID1 = pairs$IID2, IID2 = pairs$IID1, PI_HAT = pairs$PI_HAT)
all_pairs <- rbind(pairs, pairs_mirrored)

samples <- sort(unique(c(pairs$IID1, pairs$IID2)))
self_pairs <- data.frame(IID1 = samples, IID2 = samples, PI_HAT = 1.0)
all_pairs <- rbind(all_pairs, self_pairs)

all_pairs$IID1 <- factor(all_pairs$IID1, levels = samples)
all_pairs$IID2 <- factor(all_pairs$IID2, levels = rev(samples))

pihat_plot <- ggplot(all_pairs, aes(x = IID1, y = IID2, fill = PI_HAT)) +
  geom_tile(color = "black") + 
  geom_text(aes(label = round(PI_HAT, 2)), color = "black", size = 5) + 
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "PI_HAT\n(Identity)") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(title = "Clonal Lineages", subtitle = "Identity-by-Descent (PI_HAT) Matrix", x = "", y = "") +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        plot.title = element_text(size = 16, face="bold"))
ggsave("Figure_1_PI_HAT_Heatmap.png", plot = pihat_plot, width = 8, height = 6, dpi = 300)

ibs_dist <- read.table("IBS.mdist", header=FALSE)
ibs_ids <- read.table("IBS.mdist.id", header=FALSE, stringsAsFactors=FALSE)
sample_names <- sort(ibs_ids$V2)
colnames(ibs_dist) <- sample_names
ibs_dist$IID1 <- sample_names

ibs_long <- ibs_dist %>% pivot_longer(cols = -IID1, names_to = "IID2", values_to = "Distance")
ibs_long$IID1 <- factor(ibs_long$IID1, levels = sample_names)
ibs_long$IID2 <- factor(ibs_long$IID2, levels = rev(sample_names))

ibs_plot <- ggplot(ibs_long, aes(x = IID1, y = IID2, fill = Distance)) +
  geom_tile(color = "black") + 
  geom_text(aes(label = round(Distance, 4)), color = "black", size = 4) + 
  scale_fill_gradient(low = "firebrick", high = "white", name = "IBS Distance\n(0 = Identical)") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(title = "Overall Genetic Distance", subtitle = "Identity-by-State (IBS) Distance Matrix", x = "", y = "") +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        plot.title = element_text(size = 16, face="bold"))
ggsave("Figure_2_IBS_Distance_Heatmap.png", plot = ibs_plot, width = 8, height = 6, dpi = 300)

pca_data <- read.table("PCA.eigenvec", header=FALSE, stringsAsFactors=FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_data)-2)))
eigenvals <- read.table("PCA.eigenval", header=FALSE)$V1
pve <- round((eigenvals / sum(eigenvals)) * 100, 2)

pca_data$Lineage <- ifelse(pca_data$IID %in% c("S01", "S02", "S03"), "Lineage 1", "Lineage 2")

upgraded_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Lineage, label = IID)) +
  geom_point(size = 6, alpha = 0.7) +
  geom_text_repel(size = 5, fontface = "bold", box.padding = 0.8, point.padding = 0.5, show.legend = FALSE) + 
  scale_color_manual(values = c("Lineage 1" = "navy", "Lineage 2" = "firebrick")) +
  theme_minimal() +
  labs(title = "PCA Scatterplot",
       x = paste0("Principal Component 1 (", pve[1], "%)"), 
       y = paste0("Principal Component 2 (", pve[2], "%)"), color = "") +
  theme(axis.title = element_text(size = 12, face="bold"),
        axis.text = element_text(size = 11), plot.title = element_text(size = 16, face="bold"),
        legend.position = "bottom", legend.title = element_text(face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))
ggsave("Figure_3_Upgraded_PCA.png", plot = upgraded_pca, width = 8, height = 6, dpi = 300)
