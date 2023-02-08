library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- ""
outdir <- ""
healthy_path <- ""
samples <- "sample_list.txt"

### Find files
frequency <- read.delim(list.files(path, "LFS_fragment_freq_cfmedip", full.names = TRUE))
data_normal <- read.delim(list.files(healthy_path, "freq_cfmedip", full.names = TRUE))
data_samples <- read.delim(samples)

### Restrict size (10-320bp)
frequency <- frequency[frequency$length %in% c(30:230), ]
data_normal <- data_normal[data_normal$length %in% c(10:230), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$cfMeDIP %in% colnames(frequency), ]
row.names(frequency) <- frequency$length
frequency <- frequency[ , colnames(frequency) %in% data_samples$cfMeDIP]
frequency[is.na(frequency)] <- 0

## Set samples to graph and order factors
data_samples$cfMeDIP <- as.factor(data_samples$cfMeDIP)
data_samples$cancer_status <- factor(data_samples$cancer_status,
                                     levels = c("positive", "negative"),
                                     labels = c("Positive", "Negative"))
row.names(data_samples) <- data_samples$cfMeDIP

## Make healthy median
row.names(data_normal) <- data_normal$length
data_normal <- data_normal[, !(colnames(data_normal) == "length")]
normal_median <- rowMedians(as.matrix(data_normal))
normal_sd <- rowSds(as.matrix(data_normal))

### Make previvor median
samples_neg <- data_samples[data_samples$cancer_status == "Negative" &
                              data_samples$previous_cancer == "no", ]
data_neg <- frequency[ , colnames(frequency) %in% samples_neg$cfMeDIP]
previvor_median <- rowMedians(as.matrix(data_neg))
previvor_sd <- rowSds(as.matrix(data_neg))

### Make cancer positive median
samples_pos <- data_samples[data_samples$cancer_status == "Positive", ]
data_pos <- frequency[ , colnames(frequency) %in% samples_pos$cfMeDIP]
positive_median <- rowMedians(as.matrix(data_pos))
positive_sd <- rowSds(as.matrix(data_pos))

### Make comparisons and table
hbc_v_previvor <- (previvor_median - normal_median)/sqrt(((normal_sd^2)/ncol(data_normal)) + ((previvor_sd^2)/ncol(data_neg)))
previvor_v_pos <- (positive_median - previvor_median)/sqrt(((positive_sd^2)/ncol(data_pos)) + ((previvor_sd^2)/ncol(data_neg)))

length <- c(30:230)

data <- cbind(length, normal_median, previvor_median, positive_median, hbc_v_previvor, previvor_v_pos)
data <- as.data.frame(data)

data$fill <- ifelse(data$hbc_v_previvor > 0, "red", "blue")
data$fill2 <- ifelse(data$previvor_v_pos > 0, "red", "blue")

### Set facets
data$size <- ifelse(data$length <= 130, "1", "2")
data$size <- factor(data$size,
                    levels = c("1", "2"),
                    labels = c("30-130bp", "131-230bp"))


### Plot differences
theme <- theme(plot.title = element_text(hjust = 0.5, size = 12), 
               axis.title = element_text(size = 12),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.text = element_text(size = 10),
               legend.position = "none",
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               strip.background = element_rect(fill = "white"),
               strip.text = element_text(size = 10),
               panel.spacing = unit(0.1, "lines"))

plot_freq <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  geom_vline(data = filter(data, size == "131-230bp"), aes(xintercept = 147), linetype = "dashed", color = "black") +
  ggtitle("Fragment Frequency Distribution") +
  xlab("") + 
  ylab("Frequency (%)") +
  facet_grid(.~size, scales = "free", space = "free") +
  theme +
  theme(legend.position = "top",
        legend.key = element_blank()) +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS Healthy", "LFS Active Cancer")) +
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275)) +
  scale_y_continuous(limits=c(0, 3), expand = c(0,0))
plot_freq

plot_inset <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  theme +
  theme(legend.position = "none",
        legend.key = element_blank(),
        plot.margin = margin(-0.5, 0, -0.5, -0.5, "cm"),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS Healthy", "LFS Active Cancer")) +
  scale_x_continuous(limits = c(70, 130), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
plot_inset

plot_previvor <- ggplot(data) +
  geom_bar(aes(length, hbc_v_previvor, fill = fill), stat="identity") +
  ggtitle("LFS Healthy vs Healthy Control") +
  xlab("") + 
  ylab("Z-score") +
  facet_grid(.~size, scales = "free", space = "free") +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99")) +
  theme +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275)) +
  scale_y_continuous(limits=c(-3.5, 3.5), expand = c(0,0))
plot_previvor

plot_positive <- ggplot(data) +
  geom_bar(aes(length, previvor_v_pos, fill = fill2), stat="identity") +
  ggtitle("LFS Active Cancer vs LFS Healthy") +
  xlab("Fragment Size") + 
  ylab("Z-score") +
  facet_grid(.~size, scales = "free", space = "free") +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99")) +
  theme +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275)) +
  scale_y_continuous(limits=c(-3.5, 3.5), expand = c(0,0))
plot_positive

Figure <- ggarrange(plot_freq + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank()),
                    plot_previvor + theme(axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          axis.title.x = element_blank()), 
                    plot_positive, 
                    ncol = 1, align = "v", heights = c(1, 0.75, 0.9))

plot_combined <- ggdraw() +
  draw_plot(Figure) +
  draw_plot(plot_inset, x = 0.13, y = .67, width = .33, height = .18)
plot_combined

ggsave(file.path(outdir, "fragment_distributions_cfmedip.pdf"), plot_combined,  width = 4, height = 8)

### Melt data to plot individual curves
data_melt <- frequency
data_melt$length <- c(30:230)
data_melt <- reshape2::melt(data_melt, id = "length")
data_melt <- merge(data_melt, data_samples, by.x = "variable", by.y = "cfMeDIP")

data_medians <- as.data.frame(cbind(normal_median, previvor_median))
data_medians$length <- c(30:230)

### Plot individual curves
plot_curves1 <- ggplot(data_melt, aes(x = length, y = value)) +
  geom_line(aes(group = variable), size = 0.05, alpha = 0.75) +
  geom_line(data = data_medians, aes(x = length, y = normal_median, color = "red"), size = 0.5) +
  geom_vline(aes(xintercept = 167), linetype = "dashed", color = "black") +
  facet_grid(.~cancer_status) +
  ggtitle("vs Healthy Control") +
  xlab("") + 
  ylab("Frequency (%)") +
  theme +
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275)) +
  scale_y_continuous(limits=c(0, 3.25), expand = c(0,0))
plot_curves1

plot_curves2 <- ggplot(data_melt, aes(x = length, y = value)) +
  geom_line(aes(group = variable), size = 0.05, alpha = 0.75) +
  geom_line(data = data_medians, aes(x = length, y = previvor_median, color = "red"), size = 0.5) +
  geom_vline(aes(xintercept = 167), linetype = "dashed", color = "black") +
  facet_grid(.~cancer_status) +
  ggtitle("vs LFS Healthy") +
  xlab("") + 
  ylab("Frequency (%)") +
  theme +
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275)) +
  scale_y_continuous(limits=c(0, 3.25), expand = c(0,0))
plot_curves2

figure_distributions <- ggarrange(plot_curves1 + theme(axis.text.x = element_blank(),
                                                       axis.ticks.x = element_blank(),
                                                       axis.title.x = element_blank()), 
                                  plot_curves2,
                                  ncol = 1, align = "v", heights = c(0.8, 1))
figure_distributions
ggsave(file.path(outdir, "fragment_curves_cfmedip.pdf"), figure_distributions,  width = 6, height = 5)
