library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpmisc)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/fragment_score"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_score"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/fragment_score"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Import data
score <- read.delim(list.files(path, "panel", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "panel", full.names = TRUE))
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
data_samples <- data_samples[!(data_samples$TS %in% exclude), ]
data_samples <- data_samples[data_samples$TS %in% score$sample, ]
score <- score[score$sample %in% data_samples$TS, ]

## Merge scores
scores <- merge(data_samples, score, by.x = "TS", by.y = "sample", all = TRUE)
scores <- bind_rows(scores, normal)
scores[is.na(scores)] <- "Healthy"

scores$cancer_status <- factor(scores$cancer_status, levels = c("Healthy", "negative", "positive"),
                                  labels = c("Healthy", "LFS Cancer\nNegative", "LFS Cancer\nPositive"))

## Calculate statistics (cancer status)
data_stats_cancer <- scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(score, na.rm = TRUE),
                   Mean=mean(score, na.rm = TRUE),
                   SD=sd(score, na.rm = TRUE),
                   N=n())
a <- t.test(scores$score[scores$cancer_status=="Healthy"], scores$score[scores$cancer_status=="LFS Cancer\nNegative"])$p.value
b <- t.test(scores$score[scores$cancer_status=="Healthy"], scores$score[scores$cancer_status=="LFS Cancer\nPositive"])$p.value
c <- t.test(scores$score[scores$cancer_status=="LFS Cancer\nNegative"], scores$score[scores$cancer_status=="LFS Cancer\nPositive"])$p.value
t_test <- c(1, a, b)
t_test2 <- c(1, 1, c)
data_stats_cancer$pvalue <- t_test
data_stats_cancer$pvalue2 <- t_test2
data_stats_cancer$annot <- ifelse(data_stats_cancer$pvalue < 0.05 & data_stats_cancer$pvalue > 0.01, "*",
                           ifelse(data_stats_cancer$pvalue < 0.01 & data_stats_cancer$pvalue > 0.001, "**",
                                  ifelse(data_stats_cancer$pvalue < 0.001, "***", "")))
data_stats_cancer$annot2 <- ifelse(data_stats_cancer$pvalue2 < 0.05 & data_stats_cancer$pvalue2 > 0.01, "*",
                                  ifelse(data_stats_cancer$pvalue2 < 0.01 & data_stats_cancer$pvalue2 > 0.001, "**",
                                         ifelse(data_stats_cancer$pvalue2 < 0.001, "***", "")))
rm(a,b,c,t_test,t_test2)

my_comparisons <- list( c("Healthy", "LFS Cancer\nNegative"), 
                        c("Healthy", "LFS Cancer\nPositive"),
                        c("LFS Cancer\nNegative", "LFS Cancer\nPositive"))

### Make healthy medians
healthy_medianB <- as.double(data_stats_cancer[data_stats_cancer$cancer_status == "Healthy", colnames(data_stats_cancer) == "Median"])

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 11),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot figures
fig_charm <- ggplot(scores, aes(cancer_status, score, fill = cancer_status)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_cancer, aes(x = cancer_status, y = -0.95, label = N), size = 4) +
  xlab("") + 
  ylab("TP53 Fragment Score") +
  ggtitle("TP53 TS") + 
  scale_fill_manual(values = c("grey", "#a6cee3", "#fb9a99")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(-1, 0.3), expand = c(0,0)) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     size = 5,
                     tip.length = 0,
                     step.increase = 0.11,
                     position = 0.2)
fig_charm
ggsave(file.path(outdir, "fragment_score_panel.pdf"), width = 2.5, height = 4)

