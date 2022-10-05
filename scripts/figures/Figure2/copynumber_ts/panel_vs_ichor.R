library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(data.table)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/panel_cnv"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

### Read in files
samples <- read.delim(samples)
cnmops <- read.delim(list.files(path, "CHARM_LFS_panelCNmops.txt", full.names = TRUE, recursive = TRUE))
viscap <- read.delim(list.files(path, "CHARM_LFS_VisCap.txt", full.names = TRUE, recursive = TRUE))
ichor <- read.delim(list.files(path, "CHARM_LFS_corrected_depths.txt", full.names = TRUE, recursive = TRUE))

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
samples <- samples[!(samples$TS) %in% exclude, ]

### Format cnmops data
cnmops <- cnmops[cnmops$Chromosome == "chr17" & cnmops$Gene == "TP53", ]
row.names(cnmops) <- paste0(cnmops$Chromosome, "_", cnmops$Start)
cnmops <- cnmops[, colnames(cnmops) %in% samples$TS]
samples_cnmops <- samples[samples$TS %in% colnames(cnmops), ]

### Format viscap data
viscap <- viscap[viscap$Chromosome == "chr17" & viscap$Start < 8000000, ]
row.names(viscap) <- paste0(viscap$Chromosome, "_", viscap$Start - 1)
viscap <- viscap[, colnames(viscap) %in% samples$TS]
samples_viscap <- samples[samples$TS %in% colnames(viscap), ]

### Format ichorCNA data
ichor <- ichor[ichor$chr == "17" &
                 ichor$start == 7000001, colnames(ichor) %in% samples_cnmops$sWGS]
ichor <- t(ichor)

### Make cnmops and viscap comparison
cnmops_melt <- cnmops
cnmops_melt$location <- row.names(cnmops_melt)
cnmops_melt <- melt(cnmops_melt, id = "location")

viscap_melt <- viscap
viscap_melt$location <- row.names(viscap_melt)
viscap_melt <- melt(viscap_melt, id = "location")

data_comp <- merge(cnmops_melt, viscap_melt, by = c("location", "variable"))
colnames(data_comp) <- c("location", "sample", "cnmops", "viscap")

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               strip.background = element_rect(fill = NA),
               strip.text = element_text(size = 10),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.background = element_blank(),
               legend.position = "none",
               legend.direction="vertical",
               axis.text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.title = element_text(size = 10))

### Plot comparison
comp_plot <- ggplot(data_comp, aes(cnmops, viscap)) + 
  geom_point(stroke = 0, pch = 21, alpha = 0.5, size = 2, fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("") +
  xlab("PanelCNmops (CopyNumber)") +
  ylab("VisCap (CopyNumber)") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.4, aes(label = ..rr.label..)) +
  theme + 
  scale_fill_manual(values = c("grey", "#fa4e0c")) +
  scale_y_continuous(limits = c(-2.1, 1)) + 
  scale_x_continuous(limits = c(-2.1, 1)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
comp_plot

### Make ichor comparison table
cnmops_mean <- colMeans(cnmops)
viscap_mean <- colMeans(viscap)

means <- cbind(cnmops_mean, viscap_mean)
means <- merge(means, samples, by.x = "row.names", by.y = "TS")
means <- means[, c("Row.names", "cnmops_mean", "viscap_mean", "sWGS")]
means <- merge(means, ichor, by.x = "sWGS", by.y = "row.names")
colnames(means) <- c("sWGS", "TS", "cnmops", "viscap", "ichor")

### Plot comparisons
cnmop_plot <- ggplot(means, aes(cnmops, ichor)) + 
  geom_point(stroke = 0, pch = 21, alpha = 0.5, size = 2, fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("") +
  xlab("PanelCNmops (Mean)") +
  ylab("IchorCNA (CopyNumber)") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.35, aes(label = ..rr.label..)) +
  theme + 
  scale_fill_manual(values = c("grey", "#fa4e0c")) +
  scale_y_continuous(limits = c(-0.5, 0.5)) + 
  scale_x_continuous(limits = c(-1.1, 0.5)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
cnmop_plot

vis_plot <- ggplot(means, aes(viscap, ichor)) + 
  geom_point(stroke = 0, pch = 21, alpha = 0.5, size = 2, fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("") +
  xlab("VisCap (Mean)") +
  ylab("IchorCNA (CopyNumber)") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.35, aes(label = ..rr.label..)) +
  theme + 
  scale_fill_manual(values = c("grey", "#fa4e0c")) +
  scale_y_continuous(limits = c(-0.5, 0.5)) + 
  scale_x_continuous(limits = c(-1.1, 0.5)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
vis_plot

figure <- ggarrange(comp_plot, cnmop_plot, vis_plot, align = "v", nrow = 3)
figure
ggsave(file.path(outdir, "panel_cnv_comparisons.pdf"), figure, height = 7, width = 3)
