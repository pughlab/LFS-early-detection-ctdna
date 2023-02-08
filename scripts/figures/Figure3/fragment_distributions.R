library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_distributions"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

### Find files
frequency <- list.files(path, "LFS_fragment_freq.txt", full.names = TRUE)
normal_freq <- list.files(healthy_path, "freq.txt", full.names = TRUE)

## Import data
data_frequency <- read.delim(frequency)
data_normal <- read.delim(normal_freq)
data_samples <- read.delim(samples)

### Restrict size (10-320bp)
data_frequency <- data_frequency[data_frequency$length %in% c(10:320), ]
data_normal <- data_normal[data_normal$length %in% c(10:320), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
row.names(data_frequency) <- data_frequency$length
data_frequency <- data_frequency[ , colnames(data_frequency) %in% data_samples$sWGS]
data_frequency[is.na(data_frequency)] <- 0

## Set samples to graph and order factors
data_samples$sWGS <- as.factor(data_samples$sWGS)
data_samples$cancer_status <- factor(data_samples$cancer_status,
                                     levels = c("positive", "negative"),
                                     labels = c("Positive", "Negative"))
row.names(data_samples) <- data_samples$sWGS

## Make healthy median
row.names(data_normal) <- data_normal$length
data_normal <- data_normal[, !(colnames(data_normal) == "length")]
normal_median <- rowMedians(as.matrix(data_normal))
normal_sd <- rowSds(as.matrix(data_normal))

### Make previvor median
samples_neg <- data_samples[data_samples$cancer_status == "Negative" &
                              data_samples$previous_cancer == "no", ]
data_neg <- data_frequency[ , colnames(data_frequency) %in% samples_neg$sWGS]
previvor_median <- rowMedians(as.matrix(data_neg))
previvor_sd <- rowSds(as.matrix(data_neg))

### Make cancer positive median
samples_pos <- data_samples[data_samples$cancer_status == "Positive", ]
data_pos <- data_frequency[ , colnames(data_frequency) %in% samples_pos$sWGS]
positive_median <- rowMedians(as.matrix(data_pos))
positive_sd <- rowSds(as.matrix(data_pos))

### Calculate mean, median, and mode fragment lengths
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

meanfrag_normal <- mean(rep(c(10:320), round(normal_median*1000000, 0)))
meanfrag_previvor <- mean(rep(c(10:320), round(previvor_median*1000000, 0)))
meanfrag_positive <- mean(rep(c(10:320), round(positive_median*1000000, 0)))

medfrag_normal <- median(rep(c(10:320), round(normal_median*1000000, 0)))
medfrag_previvor <- median(rep(c(10:320), round(previvor_median*1000000, 0)))
medfrag_positive <- median(rep(c(10:320), round(positive_median*1000000, 0)))

modefrag_normal <- getmode(rep(c(10:320), round(normal_median*1000000, 0)))
modefrag_previvor <- getmode(rep(c(10:320), round(previvor_median*1000000, 0)))
modefrag_positive <- getmode(rep(c(10:320), round(positive_median*1000000, 0)))

### Make comparisons and table
hbc_v_previvor <- (previvor_median - normal_median)/sqrt(((normal_sd^2)/ncol(data_normal)) + ((previvor_sd^2)/ncol(data_neg)))
previvor_v_pos <- (positive_median - previvor_median)/sqrt(((positive_sd^2)/ncol(data_pos)) + ((previvor_sd^2)/ncol(data_neg)))

length <- c(10:320)

data <- cbind(length, normal_median, previvor_median, positive_median, hbc_v_previvor, previvor_v_pos)
data <- as.data.frame(data)

data$fill <- ifelse(data$hbc_v_previvor > 0, "red", "blue")
data$fill2 <- ifelse(data$previvor_v_pos > 0, "red", "blue")

### Set facets
data$size <- ifelse(data$length <= 89, "1",
                    ifelse(data$length >=90 & data$length <=150, "2",
                           ifelse(data$length >=151 & data$length <=220, "3", "4")))
data$size <- factor(data$size,
                    levels = c("1", "2", "3", "4"),
                    labels = c("10-89bp", "90-150bp", "151-220bp", "221-320bp"))


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
               strip.background = element_blank(),
               strip.text = element_text(size = 9),
               panel.spacing = unit(0.1, "lines"))

plot_freq <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  geom_vline(data = filter(data, size == "151-220bp"), aes(xintercept = 167), linetype = "dashed", color = "black") +
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
  scale_x_continuous(breaks = c(0, 25, 75, 125, 175, 225, 275, 325)) +
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
  scale_x_continuous(limits = c(90, 150), expand = c(0, 0)) +
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
  scale_y_continuous(limits=c(-5, 5), expand = c(0,0))
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
  scale_y_continuous(limits=c(-5, 5), expand = c(0,0))
plot_positive

Figure <- ggarrange(plot_freq + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank()),
                    plot_previvor + theme(axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          axis.title.x = element_blank()), 
                    plot_positive, 
                    ncol = 1, align = "v", heights = c(1, 0.75, 0.85))

plot_combined <- ggdraw() +
  draw_plot(Figure) +
  draw_plot(plot_inset, x = 0.11, y = .66, width = .32, height = .2)
plot_combined

ggsave(file.path(outdir, "fragment_distributions.pdf"), plot_combined,  width = 6, height = 8)

### Melt data to plot individual curves
data_melt <- data_frequency
data_melt$length <- c(10:320)
data_melt <- reshape2::melt(data_melt, id = "length")
data_melt <- merge(data_melt, data_samples, by.x = "variable", by.y = "sWGS")

data_medians <- as.data.frame(cbind(normal_median, previvor_median))
data_medians$length <- c(10:320)

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
ggsave(file.path(outdir, "fragment_curves.pdf"), figure_distributions,  width = 6.25, height = 5)
