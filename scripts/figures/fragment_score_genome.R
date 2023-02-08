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
ichor <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA/CHARM_LFS_ichorCNA_summary_reviewed.txt"
mutation <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot/Oncoplot_full.txt"
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Import data
score_wg <- read.delim(list.files(path, "genome", full.names = TRUE))
score_ts <- read.delim(list.files(path, "panel", full.names = TRUE))
normal_wg <- read.delim(list.files(healthy_path, "genome", full.names = TRUE))
data_ichor <- read.delim(ichor)
data_mutation <- read.delim(mutation)
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
failed_TS <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
failed_WG <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% failed_WG), ]
data_samples <- data_samples[!(data_samples$TS %in% failed_TS), ]
data_samples <- data_samples[data_samples$sWGS %in% score_wg$sample, ]
score_wg <- score_wg[score_wg$sample %in% data_samples$sWGS, ]
score_ts <- score_ts[score_ts$sample %in% data_samples$TS, ]

### Make detection cutoff
ts_limit <- quantile(score_ts$score[score_ts$sample %in% data_samples$TS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"]], 0.9)
wg_limit <- quantile(score_wg$score[score_wg$sample %in% data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"]], 0.9)

## Merge scores
colnames(score_wg) <- c("sample", "genome")
score_wg$diag <- "LFS"
score_wg <- merge(score_wg, data_samples, by.x = "sample", by.y = "sWGS")

colnames(score_ts) <- c("sample", "tp53")
scores <- merge(score_wg, score_ts, by.x = "TS", by.y = "sample", all = TRUE)
colnames(normal_wg) <- c("sample", "genome")
scores <- bind_rows(scores, normal_wg)
scores$diag[is.na(scores$diag)] <- "Healthy"
scores$cancer_status[is.na(scores$cancer_status)] <- "Healthy"

scores$cancer_status <- factor(scores$cancer_status, levels = c("Healthy", "negative", "positive"),
                                  labels = c("Healthy", "LFS Cancer\nNegative", "LFS Cancer\nPositive"))

### Make ichor vs fragment table
data_ichor <- merge(data_ichor[, c(1:6)], scores, by = "sample")
data_ichor$TF_all <- ifelse(data_ichor$TF_short > data_ichor$TF, data_ichor$TF_short, data_ichor$TF)
data_ichor$cancer_status <- factor(data_ichor$cancer_status, levels = c("LFS Cancer\nNegative", "LFS Cancer\nPositive"),
                                   labels = c("Negative", "Positive"))

### Make mutation fragment table
data_mutation <- merge(data_mutation, scores, by.x = "sample_ID", by.y = "TS")
data_mutation$cancer_status <- factor(data_mutation$cancer_status, levels = c("LFS Cancer\nNegative", "LFS Cancer\nPositive"),
                                      labels = c("Negative", "Positive"))

### Seperate LFS patients only
data_comp <- scores[scores$diag == "LFS", ]
data_comp <- data_comp[complete.cases(data_comp), ]
data_comp$cancer_status <- factor(data_comp$cancer_status, levels = c("LFS Cancer\nNegative", "LFS Cancer\nPositive"),
                                  labels = c("Cancer Negative", "Cancer Positive"))

## Calculate statistics (cancer status)
data_stats_cancer <- scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(genome, na.rm = TRUE),
                   Mean=mean(genome, na.rm = TRUE),
                   SD=sd(genome, na.rm = TRUE),
                   N=n())
a <- t.test(scores$genome[scores$cancer_status=="Healthy"], scores$genome[scores$cancer_status=="LFS Cancer\nNegative"])$p.value
b <- t.test(scores$genome[scores$cancer_status=="Healthy"], scores$genome[scores$cancer_status=="LFS Cancer\nPositive"])$p.value
c <- t.test(scores$genome[scores$cancer_status=="LFS Cancer\nNegative"], scores$genome[scores$cancer_status=="LFS Cancer\nPositive"])$p.value
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

### Set healthy median
healthy_medianB <- data_stats_cancer$Median[data_stats_cancer$cancer_status == "Healthy"]

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
fig_charm <- ggplot(scores, aes(cancer_status, genome, fill = cancer_status)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_cancer, aes(x = cancer_status, y = -0.95, label = N), size = 4) +
  xlab("") + 
  ylab("Genome Fragment Score") +
  ggtitle("Genome sWGS") + 
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
ggsave(file.path(outdir, "fragment_score_genome.pdf"), width = 2.5, height = 4)

### Compare fragment scores (panel vs swgs)
fig_comp <- ggplot(data_comp, aes(genome, tp53)) + 
  geom_point(aes(fill = cancer_status), alpha = 0.5, stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  stat_regline_equation(label.y = -0.8, label.x = -0.35, aes(label = ..rr.label..), size = 5) +
  geom_vline(xintercept = wg_limit, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = ts_limit, linetype = "dashed", size = 0.5) +
  ggtitle("Fragment Score") +
  xlab("Genome Score") +
  ylab("TP53 Score") +
  labs(fill = "Cancer Status") +
  theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = c(0.25, 0.85),
        legend.title = element_text(size = 11),
        legend.background = element_blank()) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(-1, 0), expand = c(0,0)) + 
  scale_x_continuous(limits = c(-1, 0)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_comp

ggsave(file.path(outdir, "fragment_score_comparison.pdf"), fig_comp, device = "pdf", width = 4, height = 4, units = "in")

### Compare fragment scores to ichorCNA
fig_wg <- ggplot(data_ichor, aes(TF_all, genome)) + 
  geom_point(aes(fill = cancer_status), alpha = 0.5, stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  geom_vline(xintercept = 0.03, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = wg_limit, linetype = "dashed", size = 0.5) +
  stat_regline_equation(label.y = -0.15, label.x = 0.5, aes(label = ..rr.label..), size = 5) +
  ggtitle("Genome Score vs IchorCNA") +
  xlab("Tumor Fraction") +
  ylab("Genome Score") +
  labs(fill = "Cancer Status") +
  theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.background = element_blank()) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(-1, 0), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 0.75)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_wg

fig_ts <- ggplot(data_ichor, aes(TF_all, tp53)) + 
  geom_point(aes(fill = cancer_status), alpha = 0.5, stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  geom_vline(xintercept = 0.03, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = wg_limit, linetype = "dashed", size = 0.5) +
  stat_regline_equation(label.y = -0.15, label.x = 0.5, aes(label = ..rr.label..), size = 5) +
  ggtitle("TP53 Score vs IchorCNA") +
  xlab("Tumor Fraction") +
  ylab("TP53 Score") +
  labs(fill = "Cancer Status") +
  theme + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(-1, 0), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 0.75)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_ts

fig_mut <- ggplot(data_mutation, aes(vaf, tp53)) + 
  geom_point(aes(fill = cancer_status), alpha = 0.5, stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = ts_limit, linetype = "dashed", size = 0.5) +
  stat_regline_equation(label.y = -0.15, label.x = 0.1, aes(label = ..rr.label..), size = 5) +
  ggtitle("TP53 Score vs Mutation VAF") +
  xlab("TP53 Mutation VAF") +
  ylab("TP53 Score") +
  labs(fill = "Cancer Status") +
  theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.background = element_blank()) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(-1, 0), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 0.15)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_mut

Figure <- ggarrange(fig_wg, fig_ts, fig_mut, nrow = 1, widths = c(1, 1, 1.3))
Figure

ggsave(file.path(outdir, "fragment_score_genomic_comparison.pdf"), Figure, device = "pdf", width = 13, height = 4, units = "in")







