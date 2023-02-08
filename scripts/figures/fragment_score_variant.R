library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/fragment_score"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_score"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Import data
germline <- read.delim(list.files(path, "germline", full.names = TRUE))
somatic <- read.delim(list.files(path, "somatic", full.names = TRUE))
oicr <- read.delim(list.files(path, "oicr", full.names = TRUE))
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
data_samples <- data_samples[!(data_samples$TS %in% exclude), ]
data_samples <- data_samples[data_samples$TS %in% germline$sample, ]
germline <- germline[germline$sample %in% data_samples$TS, ]
somatic <- somatic[somatic$sample %in% data_samples$TS, ]
oicr <- oicr[oicr$sample %in% data_samples$TS, ]

### Prune germline variants
germline$vaf <- germline$COUNT_VAR/(germline$COUNT_WT + germline$COUNT_VAR)
germline <- germline[germline$COUNT_WT > 100, ]
germline <- germline[germline$vaf > 0.45 & germline$vaf < 0.55, ]
germline <- germline[germline$LENGTH_WT >= 150 & germline$LENGTH_WT <= 185, ]

### Make germline difference histogram
germline$difference <- germline$LENGTH_WT - germline$LENGTH_VAR
germline_difference_median <- median(germline$difference)
germline_difference_mad <- mad(germline$difference)

### Find the center point and covariance matrix
germline_center <- colMeans(germline[,c("WFS", "VFS")])
germline_cov <- cov(germline[,c("WFS", "VFS")])

### Find ellipses limits and coordinates
germline_rad <- qchisq(p = 0.95, df = ncol(germline))
germline_rad <- sqrt(germline_rad)
germline_ellipse <- car::ellipse(center = germline_center , shape = germline_cov , radius = germline_rad ,
                        segments = 150 , draw = FALSE)
germline_ellipse <- as.data.frame(germline_ellipse)

### Calculate the Mahalanobis Distances
cutoff <- qchisq(p = 0.95, df = ncol(germline))
germline$score <- mahalanobis(germline[, c("WFS", "VFS")], germline_center, germline_cov)
germline$colour <- ifelse(germline$score > cutoff, "red", "black")

somatic$score <- mahalanobis(somatic[, c("WFS", "VFS")], germline_center, germline_cov)
somatic$colour <- ifelse(somatic$score > cutoff, "red", "black")

### Prune somatic mutations
pre_prune <- nrow(somatic)
somatic$vaf <- somatic$COUNT_VAR/(somatic$COUNT_WT + somatic$COUNT_VAR)
somatic <- somatic[somatic$COUNT_WT > 100, ]
somatic <- somatic[!((somatic$vaf > 0.45 & somatic$vaf < 0.55) |
                       somatic$vaf > 0.9), ]
somatic <- somatic[somatic$LENGTH_WT >= 150 & somatic$LENGTH_WT <= 185, ]
somatic$difference <- somatic$LENGTH_WT - somatic$LENGTH_VAR

somatic_filtered <- somatic[#somatic$difference > germline_difference_median + 3*germline_difference_mad | 
                            #  somatic$difference < germline_difference_median - 3*germline_difference_mad |
                            #  somatic$t_test < 0.05 |
                            #  somatic$ks_test < 0.05 |
                              somatic$score > cutoff, ]

### Calculate difference in somatic lengths
somatic_difference_median <- median(somatic$difference)
somatic_difference_mad <- mad(somatic$difference)

### Combine dataframe and calculate ks test
data <- bind_rows(germline, somatic)

ks <- ks.test(data$difference[data$type == "germline"],
              data$difference[data$type == "somatic"])$p.value
ks <- format(ks, nsmall = 2)

data_stats <- data %>%
  group_by(type) %>%
  dplyr::summarise(mean = mean(VFS),
                   sd = sd(VFS),
                   N = n())
data_stats$pvalue <- c(0, t.test(data$VFS[data$type == "germline"], data$VFS[data$type == "somatic"])$p.value)
data_stats$annot <- c("", "*")

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot the germline variant scores w/ 95% CI
fig_germ <- ggplot() +
  geom_point(data = germline, aes(WFS, VFS, color = colour), size = 1, alpha = 0.5) +
  geom_polygon(data = germline_ellipse, aes(x,y), fill = NA , color = "red" , alpha = 0.5)+
  geom_point(aes(germline_center[1] , germline_center[2]) , size = 3 , color = "black") +
  geom_text(aes(x = -1.5, y = 1, label = paste0("n = ", nrow(germline)))) +
  stat_regline_equation(data = germline, label.x = -1.9, label.y = 1.5, aes(WFS, VFS, label = ..rr.label..)) +
  geom_abline(slope = 1, intercept = c(0,0)) +
  ylab("Mutant Fragment Score") + 
  xlab("Wildtype Fragment Score") +
  ggtitle("Germline Variants") + 
  theme +
  scale_color_manual(values = c(black = "grey", red = "red")) +
  scale_y_continuous(limits = c(-2.55, 2.55)) +
  scale_x_continuous(limits = c(-2.1, 1)) 
fig_germ

fig_somatic <- ggplot() +
  geom_point(data = somatic, aes(WFS, VFS, color = colour), size = 1, alpha = 0.5) +
  geom_polygon(data = germline_ellipse, aes(x,y), fill = NA , color = "red" , alpha = 0.5)+
  geom_point(aes(germline_center[1] , germline_center[2]) , size = 3 , color = "black") +
  geom_text(aes(x = -1.5, y = 1, label = paste0("n = ", nrow(somatic)))) +
  stat_regline_equation(data = somatic, label.x = -1.9, label.y = 1.5, aes(WFS, VFS, label = ..rr.label..)) +
  geom_abline(slope = 1, intercept = c(0,0)) +
  ylab("Mutant Fragment Score") + 
  xlab("Wildtype Fragment Score") +
  ggtitle("Somatic Variants") + 
  theme +
  scale_color_manual(values = c(black = "grey", red = "red")) +
  scale_y_continuous(limits = c(-2.55, 2.55)) +
  scale_x_continuous(limits = c(-2.1, 1)) 
fig_somatic

fig_scores <- ggplot(data, aes(type, VFS, fill = type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5,
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5, alpha = 0.5) +
  geom_text(data = data_stats, aes(x = type, y = 2.5, label = annot), size = 7) +
  ylab("Variant Score") + 
  xlab("Variant Type") +
  ggtitle("Score") +
  scale_fill_manual(values = c(germline = "black", somatic = "red"),) +
  theme +
  theme(axis.text.x = element_blank())
fig_scores

Figure <- ggarrange(fig_scores, fig_germ, fig_somatic, nrow = 1, widths = c(1.5, 2, 2), align = "hv")
Figure

ggsave(file.path(outdir, "fragment_score_variant_scores.pdf"), Figure, width = 9, height = 3.5)


fig_hist <- ggplot() +
  geom_density(data = data, aes(y = stat(density), x = difference, fill = type), alpha = 0.25, color = NA) +
  geom_text(aes(x = -15, y = 0.15, label = paste0("KS Test\np-value = ", ks))) +
  ylab("Frequency") + 
  xlab("Fragment Size Difference") +
  labs(fill = "Variant Type") +
  ggtitle("Mean Fragment Length") +
  scale_fill_manual(values = c(germline = "black", somatic = "red"),
                    labels = c(germline = "Germline", somatic = "Somatic")) +
  theme + 
  theme(legend.position = c(0.8, 0.75)) +
  scale_x_continuous(limits = c(-25, 25), expand = c(0,0))
fig_hist

ggsave(file.path(outdir, "fragment_score_variant_difference.pdf"), fig_hist, width = 3.5, height = 3.5)
