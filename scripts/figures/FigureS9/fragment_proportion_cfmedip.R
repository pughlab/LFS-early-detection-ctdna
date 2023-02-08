library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(ggh4x)

### Set paths
path <- ""
outdir <- ""
healthy_path <- ""
normal_samples <- "hbc_sample_list.txt"
samples <- "sample_list.txt"

## Import data
cm <- read.delim(list.files(path, "LFS_fragment_cfmedip.txt", full.names = TRUE))
normal_cm <- read.delim(list.files(healthy_path, "fragment_cfmedip.txt", full.names = TRUE))
wg <- read.delim(list.files(path, "LFS_fragment_freq.txt", full.names = TRUE))
normal_wg <- read.delim(list.files(healthy_path, "fragment_freq.txt", full.names = TRUE))
normal_samples <- read.delim(normal_samples)
samples <- read.delim(samples)

### Restrict sizes (10-500bp)
length_wg <- c(10:500)
length_cm <- c(30:500)

cm <- cm[cm$length %in% length_cm, ]
normal_cm <- normal_cm[normal_cm$length %in% length_cm, ]

wg <- wg[wg$length %in% length_wg, ]
normal_wg <- normal_wg[normal_wg$length %in% length_wg, ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]
samples <- samples[(samples$cancer_status %in% c("positive", "negative")), ]
samples <- samples[samples$cfMeDIP %in% colnames(cm), ]

row.names(cm) <- cm$length
cm <- cm[ , colnames(cm) %in% samples$cfMeDIP]
cm[is.na(cm)] <- 0

row.names(normal_cm) <- normal_cm$length
normal_cm <- normal_cm[, !(colnames(normal_cm) == "length")]

row.names(wg) <- wg$length
wg <- wg[ , colnames(wg) %in% samples$sWGS]
wg[is.na(wg)] <- 0

row.names(normal_wg) <- normal_wg$length
normal_wg <- normal_wg[, !(colnames(normal_wg) == "length")]

### Calculate proportions of size bins
frag_1 <- colSums(as.matrix(cm[row.names(cm) %in% c(30:130), ]))/colSums(cm[row.names(cm) %in% c(30:230), ])
frag_2 <- colSums(as.matrix(cm[row.names(cm) %in% c(131:230), ]))/colSums(cm[row.names(cm) %in% c(30:230), ])
frag_3 <- colSums(as.matrix(cm[row.names(cm) %in% c(231:500), ]))/colSums(cm)

norm_1 <- colSums(as.matrix(normal_cm[row.names(normal_cm) %in% c(30:130), ]))/colSums(normal_cm[row.names(normal_cm) %in% c(30:230), ])
norm_2 <- colSums(as.matrix(normal_cm[row.names(normal_cm) %in% c(131:230), ]))/colSums(normal_cm[row.names(normal_cm) %in% c(30:230), ])
norm_3 <- colSums(as.matrix(normal_cm[row.names(normal_cm) %in% c(231:500), ]))/colSums(normal_cm)

frag_wg_1 <- colSums(as.matrix(wg[row.names(wg) %in% c(10:150), ]))/colSums(wg[row.names(wg) %in% c(10:250), ])
frag_wg_2 <- colSums(as.matrix(wg[row.names(wg) %in% c(151:250), ]))/colSums(wg[row.names(wg) %in% c(10:250), ])
frag_wg_3 <- colSums(as.matrix(wg[row.names(wg) %in% c(251:500), ]))/colSums(wg)

norm_wg_1 <- colSums(as.matrix(normal_wg[row.names(normal_wg) %in% c(10:150), ]))/colSums(normal_wg[row.names(normal_wg) %in% c(10:250), ])
norm_wg_2 <- colSums(as.matrix(normal_wg[row.names(normal_wg) %in% c(151:250), ]))/colSums(normal_wg[row.names(normal_wg) %in% c(10:250), ])
norm_wg_3 <- colSums(as.matrix(normal_wg[row.names(normal_wg) %in% c(251:500), ]))/colSums(normal_wg)

### Make proportions table (cfMeDIP)
data_frag <- data.frame(sample = colnames(cm),
                        frag_1 = frag_1,
                        frag_2 = frag_2,
                        frag_3 = frag_3)
data_frag <- merge(data_frag, samples, by.x = "sample", by.y = "cfMeDIP")
data_norm <- data.frame(sample = colnames(normal_cm),
                        frag_1 = norm_1,
                        frag_2 = norm_2,
                        frag_3 = norm_3)
data <- bind_rows(data_frag, data_norm)
data[is.na(data)] <- "Healthy"
data$diag <- ifelse(data$cancer_status == "positive", "LFS Cancer\nPositive",
                    ifelse(data$cancer_status == "negative", "LFS Cancer\nNegative", "Healthy"))
data <- data[colnames(data) %in% c("sample", "frag_1", "frag_2", "frag_3", "diag")]
data_melt <- reshape2::melt(data, id = c("sample", "diag"))
data_melt$diag <- factor(data_melt$diag, levels = c("Healthy", "LFS Cancer\nNegative", "LFS Cancer\nPositive"))
data_melt$variable <- factor(data_melt$variable, levels = c("frag_1", "frag_2", "frag_3"),
                             labels = c("30-130bp", "131-230bp", "231-500bp"))

### Make proportions table (cfMeDIP)
data_frag_wg <- data.frame(sample = colnames(wg),
                           frag_wg_1 = frag_wg_1,
                           frag_wg_2 = frag_wg_2,
                           frag_wg_3 = frag_wg_3)
data_frag_wg <- merge(data_frag_wg, samples, by.x = "sample", by.y = "sWGS")
data_norm_wg <- data.frame(sample = colnames(normal_wg),
                           frag_wg_1 = norm_wg_1,
                           frag_wg_2 = norm_wg_2,
                           frag_wg_3 = norm_wg_3)
data_norm_wg <- merge(data_norm_wg, normal_samples[, 4:5], by.x = "sample", by.y = "sWGS")
data_wg <- bind_rows(data_frag_wg, data_norm_wg)
data_wg[is.na(data_wg)] <- "Healthy"
data_wg$diag <- ifelse(data_wg$cancer_status == "positive", "LFS Cancer\nPositive",
                       ifelse(data_wg$cancer_status == "negative", "LFS Cancer\nNegative", "Healthy"))
data_wg <- data_wg[colnames(data_wg) %in% c("sample", "frag_wg_1", "frag_wg_2", "frag_wg_3", "diag", "cfMeDIP")]
data_wg_melt <- reshape2::melt(data_wg, id = c("sample", "diag", "cfMeDIP"))
data_wg_melt$diag <- factor(data_wg_melt$diag, levels = c("Healthy", "LFS Cancer\nNegative", "LFS Cancer\nPositive"))
data_wg_melt$variable <- factor(data_wg_melt$variable, levels = c("frag_wg_1", "frag_wg_2", "frag_wg_3"),
                                labels = c("30-130bp", "131-230bp", "231-500bp"))

### Calculate stats
data_stats <- data_melt %>%
  group_by(diag, variable) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value),
                   N=n())

patients <- unique(data_stats$diag)
sizes <- unique(data_stats$variable)
t_test <- c()
for (patient in patients) {
  for (size in sizes) {
    a <- t.test(data_melt$value[data_melt$diag == "Healthy" & data_melt$variable == size], data_melt$value[data_melt$diag == patient & data_melt$variable == size])$p.value
    t_test <- c(t_test, a)
  }
}
data_stats$pvalue <- t_test
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))

### Set limits
data_stats$max <- ifelse(data_stats$variable == "30-130bp", 0.34,
                         ifelse(data_stats$variable == "131-230bp", 0.92, 0.35))
data_stats$min <- ifelse(data_stats$variable == "30-130bp", 0,
                         ifelse(data_stats$variable == "131-230bp", 0, -0.02))

### Plot differences
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.title = element_text(size = 13),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(size = 13),
               axis.text.x = element_blank(),
               legend.position = "none",
               legend.text = element_text(size = 13),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 13))

plot <- ggplot(data_melt) +
  geom_boxplot(aes(diag, value, fill = diag), alpha = 0.5) +
  geom_text(data = data_stats, aes(diag, max, label = annot), size = 6) +
  geom_text(data = data_stats, aes(diag, min, label = N), size = 4) +
  ggtitle("Fragment Size Proportions") +
  xlab("Cancer Status") + 
  ylab("Proportion of Fragments") +
  facet_wrap(.~variable, scales = "free_y") +
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(0, 0.35)),
                               scale_y_continuous(limits = c(0, 0.95)),
                               scale_y_continuous(limits = c(-0.02, 0.35)))) +
  theme +
  scale_fill_manual(name = " ", values =c("black", "#1F78B4", "#E31A1C"))
plot

ggsave(file.path(outdir, "fragment_proportions_cfmedip.pdf"), width = 4.5, height = 3)

### Make comparison table
data_comp <- merge(data_melt[data_melt$variable == "30-130bp", ], 
                   data_wg_melt[data_wg_melt$variable == "30-130bp", -1], 
                   by.x = c("sample", "diag", "variable"), by.y = c("cfMeDIP", "diag", "variable"))

plot_comp <- ggplot(data_comp, aes(value.x, value.y, color = diag)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.35, aes(label = ..rr.label..), color = "black") +
  ggtitle("Proportion of Short (<150bp) Fragments") +
  xlab("cfMeDIP") + 
  ylab("sWGS") +
  facet_grid(.~diag) +
  theme +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(name = " ", values =c("black", "#1F78B4", "#E31A1C")) +
  scale_x_continuous(limits = c(0.1, 0.35), breaks = c(0.1, 0.2, 0.3))
plot_comp

ggsave(file.path(outdir, "fragment_proportions_comp.pdf"), width = 6, height = 3)
