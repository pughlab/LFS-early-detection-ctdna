library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpmisc)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_proportion"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
ichor <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA/CHARM_LFS_ichorCNA_summary_reviewed.txt"
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Find files
proportion <- list.files(path, "proportion", full.names = TRUE)
normal_prop <- list.files(healthy_path, "proportion", full.names = TRUE)
normal_freq <- list.files(healthy_path, "freq", full.names = TRUE)

### Import data
data_proportion <- read.delim(proportion)
normal_prop <- read.delim(normal_prop)
data_ichor <- read.delim(ichor)
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_proportion <- data_proportion[!(data_proportion$cancer == "unknown"), ]
data_proportion <- data_proportion[!(data_proportion$sample %in% exclude), ]
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data_proportion$sample, ]

## Seperate LFS and combine with external healthy
data_prop <- data_proportion[data_proportion$type %in% c("LFS"), ]
colnames(normal_prop) <- c("sample", "cancer", "type", "p_20_150")
normal_prop$cancer <- "healthy"
colnames(data_prop) <- colnames(normal_prop)
data_prop <- rbind(data_prop, normal_prop)
data_prop$cancer <- factor(data_prop$cancer, levels = c("healthy", "negative", "positive"),
                           labels = c("Healthy", "LFS Cancer\nNegative", "LFS Cancer\nPositive"))

### Make ichor vs fragment table
data_ichor <- merge(data_ichor, data_prop, by = "sample")
data_ichor$TF_all <- ifelse(data_ichor$TF_short > data_ichor$TF, data_ichor$TF_short, data_ichor$TF)
data_ichor$cancer_status <- factor(data_ichor$cancer_status, levels = c("negative", "positive"),
                                   labels = c("Negative", "Positive"))

## Calculate statistics (cancer status)
data_stats_cancer <- data_prop %>%
  group_by(cancer)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop$p_20_150[data_prop$cancer=="Healthy"], data_prop$p_20_150[data_prop$cancer=="LFS Cancer\nNegative"])$p.value
b <- t.test(data_prop$p_20_150[data_prop$cancer=="Healthy"], data_prop$p_20_150[data_prop$cancer=="LFS Cancer\nPositive"])$p.value
c <- t.test(data_prop$p_20_150[data_prop$cancer=="LFS Cancer\nNegative"], data_prop$p_20_150[data_prop$cancer=="LFS Cancer\nPositive"])$p.value
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

### Compare against cancer history
data_prop_LFS <- data_prop[data_prop$type == "LFS", ]
data_prop_LFS <- merge(data_prop_LFS, data_samples, by.x = "sample", by.y = "sWGS")
data_prop_LFS$previous_cancer <- factor(data_prop_LFS$previous_cancer, levels = c("no", "yes"),
                                        labels = c("No", "Yes"))
data_stats_previous <- data_prop_LFS[!(data_prop_LFS$previous_cancer == ""), ] %>%
  group_by(previous_cancer)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop_LFS$p_20_150[data_prop_LFS$previous_cancer=="No"], data_prop_LFS$p_20_150[data_prop_LFS$previous_cancer=="Yes"])$p.value
t_test <- c(1, a)
data_stats_previous$pvalue <- t_test

### Compare with negative, no history LFS
samples_neg <- data_samples[data_samples$cancer_status == "negative" & 
                              data_samples$previous_cancer == "no", ]
data_prop_hbc <- data_prop[data_prop$type == "healthy", ]
data_prop_neg <- data_prop_LFS[data_prop_LFS$sample %in% samples_neg$sWGS, ]
data_prop_neg <- bind_rows(data_prop_neg, data_prop_hbc)
data_prop_neg[is.na(data_prop_neg)] <- "Healthy"
data_prop_neg$type.x <- ifelse(data_prop_neg$type.x == "LFS", "LFS Healthy", data_prop_neg$type.x)

data_stats_neg <- data_prop_neg %>%
  group_by(type.x)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop_neg$p_20_150[data_prop_neg$type.x =="Healthy"], data_prop_neg$p_20_150[data_prop_neg$type.x =="LFS Healthy"])$p.value
t_test <- c(1, a)
data_stats_neg$pvalue <- t_test
data_stats_neg$annot <- ifelse(data_stats_neg$pvalue < 0.05 & data_stats_neg$pvalue > 0.01, "*",
                                  ifelse(data_stats_neg$pvalue < 0.01 & data_stats_neg$pvalue > 0.001, "**",
                                         ifelse(data_stats_neg$pvalue < 0.001, "***", "")))

### Make healthy medians
healthy_medianB <- as.double(data_stats_cancer[data_stats_cancer$cancer == "Healthy", colnames(data_stats_cancer) == "Median"])

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
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               plot.margin = unit(c(0,1,0,1), "lines"))

### Plot figures
fig_charm <- ggplot(data_prop, aes(cancer, p_20_150, fill = cancer)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_cancer, aes(x = cancer, y = 0.03, label = N), size = 4) +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("Cancer Status") + 
  scale_fill_manual(values = c("grey", "#a6cee3", "#fb9a99")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0)) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     size = 5,
                     tip.length = 0,
                     step.increase = 0.11,
                     position = 0.2)
fig_charm

fig_neg <- ggplot(data_prop_neg, aes(type.x, p_20_150) ) +
  geom_boxplot(aes(fill = type), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = type), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_neg, aes(x = type.x, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_neg, aes(x = type.x, y = 0.4, label = annot), size = 5.5) +
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  labs(fill = "") +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("LFS\nHealthy") + 
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0))
fig_neg  

fig_previous <- ggplot(data_prop_LFS[!(data_prop_LFS$previous_cancer == ""), ], aes(previous_cancer, p_20_150, fill = previous_cancer)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = data_stats_previous, aes(x = previous_cancer, y = 0.03, label = N), size = 4) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("Previous\nCancer") + 
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0))
fig_previous

data_prop_LFS$years <- as.numeric(data_prop_LFS$years)
fig_age <- ggplot(data_prop_LFS, aes(years, p_20_150)) + 
  geom_point(aes(fill = "grey"), stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.45, aes(label = ..rr.label..), size = 5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "red") +
  ggtitle("Age (years)") +
  xlab("") +
  ylab("Proportion under 150bp") +
  labs(fill = "Cancer Status") +
  theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c("grey")) +
  scale_y_continuous(limits = c(0, 0.55), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 75)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_age

Figure <- ggarrange(fig_charm, fig_neg, fig_previous, fig_age, align = "h", widths = c(1.25, 1, 1, 1.75), nrow = 1)
Figure
ggsave(file.path(outdir, "fragment_proportions.pdf"), Figure, device = "pdf", width = 10, height = 5, units = "in")

### Compare fragment size to ichor TF
fig_ichor <- ggplot(data_ichor, aes(TF, p_20_150)) + 
  geom_point(aes(fill = cancer_status), alpha = 0.5, stroke = 0, pch = 21, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", fill = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.35, label.x = 0.2, aes(label = ..rr.label..), size = 5) +
  geom_hline(yintercept = healthy_medianB, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "grey") +
  ggtitle("Fragment size vs ichorCNA") +
  xlab("ichorCNA (TF)") +
  ylab("Proportion under 150bp") +
  labs(fill = "Cancer Status") +
  theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = c(0.75,0.2),
        legend.title = element_text(size = 11),
        legend.background = element_blank()) +
  scale_fill_manual(values = c("black", "red")) +
  scale_y_continuous(limits = c(0.1, 0.4), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 0.4)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_ichor

ggsave(file.path(outdir, "fragment_proportions_ichor.pdf"), fig_ichor, device = "pdf", width = 4, height = 3.5, units = "in")
