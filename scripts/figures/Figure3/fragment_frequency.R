library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(magrittr)
library(ggeffects)
library(sjmisc)
library(lme4)
library(splines)
library(broom)
library(psych)
library(caret)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_frequency"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

## Import data
data_frequency <- read.delim(list.files(path, "LFS_fragment_freq.txt", full.names = TRUE))
data_samples <- read.delim(samples)

### Restrict size (10-320bp)
data_frequency <- data_frequency[data_frequency$length %in% c(10:320), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
row.names(data_frequency) <- data_frequency$length
data_frequency <- data_frequency[ , colnames(data_frequency) %in% data_samples$sWGS]
data_frequency[is.na(data_frequency)] <- 0

## Set samples to graph and order factors
data_samples$sWGS <- as.factor(data_samples$sWGS)
data_samples$cancer_status <- factor(data_samples$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
data_samples$Age <- factor(data_samples$Age,
                           levels = c("adult", "pediatric"),
                           labels = c("Adult", "Pediatric"))
data_samples$previous_cancer <- factor(data_samples$previous_cancer, levels = c("yes", "no", ""),
                                       labels = c("Yes", "No", "Unknown"))
data_samples$stage <- factor(data_samples$stage, levels = c("high", "low", "", "unknown"),
                             labels = c("Stage III/IV", "Stage 0/I/II", "none", "Unknown"))
row.names(data_samples) <- data_samples$sWGS
data_samples <- data_samples[order(data_samples$sample_parent,
                                   data_samples$timepoint), ]
data_frequency <- data_frequency[, data_samples$sWGS]

### Make previvor median
samples_neg <- data_samples[data_samples$cancer_status == "Negative" &
                              data_samples$previous_cancer == "No", ]
data_neg <- data_frequency[ , colnames(data_frequency) %in% samples_neg$sWGS]
previvor_median <- rowMedians(as.matrix(data_neg))
previvor_sd <- rowSds(as.matrix(data_neg))
data_median <- cbind(previvor_median, previvor_sd)
data_median <- as.data.frame(data_median)
row.names(data_median) <- c(10:320)
colnames(data_median) <- c("median", "sd")

## Normalize to fold change over median
data_change <- (data_frequency - data_median$median)/data_median$sd
data_change <- as.matrix(t(data_change))

data_neg <- data_change[row.names(data_change) %in% samples_neg$sWGS, ]

## Calculate summed Z-scores for each compartment
sums_neg0 <- rowSums(data_neg[ , colnames(data_change) %in% c(10:89)], na.rm = TRUE)
sums0 <- rowSums(data_change[ , colnames(data_change) %in% c(10:89)], na.rm = TRUE)
mean <- mean(sums_neg0, na.rm = TRUE)
sd <- sd(sums_neg0)
zscores0 <- (sums0 - mean)/sd

sums_neg1 <- rowSums(data_neg[ , colnames(data_change) %in% c(90:150)], na.rm = TRUE)
sums1 <- rowSums(data_change[ , colnames(data_change) %in% c(90:150)], na.rm = TRUE)
mean <- mean(sums_neg1, na.rm = TRUE)
sd <- sd(sums_neg1)
zscores1 <- (sums1 - mean)/sd

sums_neg2 <- rowSums(data_neg[ , colnames(data_change) %in% c(151:160)], na.rm = TRUE) +
  rowSums(abs(data_neg[ , colnames(data_change) %in% c(161:220)]), na.rm = TRUE)
sums2 <- rowSums(data_change[ , colnames(data_change) %in% c(151:160)], na.rm = TRUE) +
  rowSums(abs(data_change[ , colnames(data_change) %in% c(161:220)]), na.rm = TRUE)
mean <- mean(sums_neg2, na.rm = TRUE)
sd <- sd(sums_neg2)
zscores2 <- (sums2 - mean)/sd

sums_neg3 <- rowSums(data_neg[ , colnames(data_change) %in% c(221:275)], na.rm = TRUE) +
  rowSums(abs(data_neg[ , colnames(data_change) %in% c(276:320)]), na.rm = TRUE)
sums3 <- rowSums(data_change[ , colnames(data_change) %in% c(221:275)], na.rm = TRUE) +
  rowSums(abs(data_change[ , colnames(data_change) %in% c(276:320)]), na.rm = TRUE)
mean <- mean(sums_neg3, na.rm = TRUE)
sd <- sd(sums_neg3)
zscores3 <- (sums3 - mean)/sd

data_scores <- cbind(zscores0, zscores1, zscores2, zscores3)
data_scores <- as.data.frame(data_scores)
data_scores$sample <- row.names(data_change)

### Compare metrics
data_comp <- data_scores
names(data_comp)[names(data_comp) == "zscores0"] <- "10-89bp"
names(data_comp)[names(data_comp) == "zscores1"] <- "90-150bp"
names(data_comp)[names(data_comp) == "zscores2"] <- "151-220bp"
names(data_comp)[names(data_comp) == "zscores3"] <- "221-320bp"

pdf(file.path(outdir, "fragment_scores_comparison.pdf"), width = 7, height = 7)
pairs.panels(data_comp[, c(1:4)])
dev.off()

### Merge scores data with metadata
data_scores <- merge(data_scores, data_samples, by.x = "sample", by.y = "sWGS")
row.names(data_scores) <- data_scores$sample
data_scores <- data_scores[data_samples$sWGS, ]

### Train a logistic regression model
set.seed(123)
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 10, verboseIter = TRUE, classProbs = TRUE, savePredictions = TRUE)
model <- caret::train(x = data_scores[2:5], y = data_scores$cancer_status, trControl = train_control, method = "glm")

### Get training results
predictions <- model[["pred"]][["Positive"]]
index <- model[["pred"]][["rowIndex"]]
predictions <- as.data.frame(cbind(predictions, index))
resamples <- data_scores[index, ]
predictions <- cbind(predictions, resamples$sample)
colnames(predictions) <- c("zscores", "index", "sample")
predictions <- predictions %>%
  group_by(sample) %>%
  dplyr::summarise(Median=median(zscores, na.rm = TRUE),
                   zscores=mean(zscores, na.rm = TRUE),
                   SD=sd(zscores, na.rm = TRUE),
                   N=n())

### Normalize predictions to previvors (z-score)
predictions_neg <- predictions[predictions$sample %in% row.names(data_neg), ]
neg_mean <- mean(predictions_neg$zscores)
neg_sd <- sd(predictions_neg$zscores)
limit <- quantile(predictions_neg$zscores, 0.95)

predictions$zscores <- (predictions$zscores - neg_mean)/neg_sd
limit <- (limit - neg_mean)/neg_sd
data_scores <- merge(data_scores, predictions[, c("sample", "zscores")], by = "sample")

### Calculate limits for each score
limit0 <- quantile(data_scores[data_scores$sample %in% row.names(data_neg), c("zscores0")], 0.95)
limit1 <- quantile(data_scores[data_scores$sample %in% row.names(data_neg), c("zscores1")], 0.95)
limit2 <- quantile(data_scores[data_scores$sample %in% row.names(data_neg), c("zscores2")], 0.95)
limit3 <- quantile(data_scores[data_scores$sample %in% row.names(data_neg), c("zscores3")], 0.95)

### Reformat Cancer status
data_scores$cancer_status <- factor(data_scores$cancer_status, levels = c("Negative", "Positive"),
                                    labels = c("Negative", "Positive"))

### Calculate statistics and n's
data_stats <- data_scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(zscores, na.rm = TRUE),
                   Mean=mean(zscores, na.rm = TRUE),
                   SD=sd(zscores, na.rm = TRUE),
                   N=n())
a <- t.test(data_scores$zscores[data_scores$cancer_status == "Negative"], data_scores$zscores[data_scores$cancer_status == "Positive"])$p.value
data_stats$pvalue <- c(1, a)
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))

data_stats0 <- data_scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(zscores0, na.rm = TRUE),
                   Mean=mean(zscores0, na.rm = TRUE),
                   SD=sd(zscores0, na.rm = TRUE),
                   N=n())
a <- t.test(data_scores$zscores0[data_scores$cancer_status =="Negative"], data_scores$zscores0[data_scores$cancer_status =="Positive"])$p.value
data_stats0$pvalue <- c(1, a)
data_stats0$annot <- ifelse(data_stats0$pvalue < 0.05 & data_stats0$pvalue > 0.01, "*",
                            ifelse(data_stats0$pvalue < 0.01 & data_stats0$pvalue > 0.001, "**",
                                   ifelse(data_stats0$pvalue < 0.001, "***", "")))

data_stats1 <- data_scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(zscores1, na.rm = TRUE),
                   Mean=mean(zscores1, na.rm = TRUE),
                   SD=sd(zscores1, na.rm = TRUE),
                   N=n())
a <- t.test(data_scores$zscores1[data_scores$cancer_status =="Negative"], data_scores$zscores1[data_scores$cancer_status =="Positive"])$p.value
data_stats1$pvalue <- c(1, a)
data_stats1$annot <- ifelse(data_stats1$pvalue < 0.05 & data_stats1$pvalue > 0.01, "*",
                               ifelse(data_stats1$pvalue < 0.01 & data_stats1$pvalue > 0.001, "**",
                                      ifelse(data_stats1$pvalue < 0.001, "***", "")))

data_stats2 <- data_scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(zscores2, na.rm = TRUE),
                   Mean=mean(zscores2, na.rm = TRUE),
                   SD=sd(zscores1, na.rm = TRUE),
                   N=n())
a <- t.test(data_scores$zscores2[data_scores$cancer_status =="Negative"], data_scores$zscores2[data_scores$cancer_status =="Positive"])$p.value
a <- round(a, 4)
data_stats2$pvalue <- c(1, a)
data_stats2$annot <- ifelse(data_stats2$pvalue < 0.05 & data_stats2$pvalue > 0.01, "*",
                               ifelse(data_stats2$pvalue < 0.01 & data_stats2$pvalue > 0.001, "**",
                                      ifelse(data_stats2$pvalue < 0.001, "***", "")))

data_stats3 <- data_scores %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(zscores3, na.rm = TRUE),
                   Mean=mean(zscores3, na.rm = TRUE),
                   SD=sd(zscores3, na.rm = TRUE),
                   N=n())
a <- t.test(data_scores$zscores3[data_scores$cancer_status =="Negative"], data_scores$zscores3[data_scores$cancer_status =="Positive"])$p.value
a <- round(a, 4)
data_stats3$pvalue <- c(1, a)
data_stats3$annot <- ifelse(data_stats3$pvalue < 0.05 & data_stats3$pvalue > 0.01, "*",
                               ifelse(data_stats3$pvalue < 0.01 & data_stats3$pvalue > 0.001, "**",
                                      ifelse(data_stats3$pvalue < 0.001, "***", "")))

### Plot Z-scores (vs HBCs)
theme <- theme(plot.title = element_text(hjust = 0.5, size = 12), 
               axis.title = element_text(size = 12),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.text = element_text(size = 12),
               legend.position = "none",
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank())

col <- c("grey65", "#FB9A99")

outliers0 <- boxplot.stats(data_scores$zscores0[data_scores$cancer_status == "Negative"])$out > 5
outliers0_neg <- length(outliers0[outliers0 == TRUE])
outliers0 <- boxplot.stats(data_scores$zscores0[data_scores$cancer_status == "Positive"])$out > 5
outliers0_pos <- length(outliers0[outliers0 == TRUE])
outliers0 <- as.data.frame(cbind(c(rep("Negative", outliers0_neg), rep("Positive", outliers0_pos)), 
                                 c(seq(4.5, 4.5 + 0.1*(outliers0_neg - 1), 0.1))))
outliers0$V2 <- as.numeric(outliers0$V2)

plot0 <- ggplot(data_scores, aes(cancer_status, zscores0)) +
  geom_boxplot(aes(fill = cancer_status), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = cancer_status), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_point(data = outliers0, aes(V1, V2), pch = 2) +
  geom_text(data = data_stats0, aes(x = cancer_status, y = -2.7, label = N), size = 4) +
  geom_text(data = data_stats0, aes(x = cancer_status, y = 4.7, label = annot), size = 6) +
  scale_fill_manual(values = col) +
  ggtitle("10-89bp") +
  xlab("") + 
  ylab("Fragment Size Z-Score") +
  theme +
  scale_y_continuous(limits=c(-3, 5), expand = c(0,0))
plot0

outliers1 <- boxplot.stats(data_scores$zscores1[data_scores$cancer_status == "Negative"])$out > 5
outliers1_neg <- length(outliers1[outliers1 == TRUE])
outliers1 <- boxplot.stats(data_scores$zscores0[data_scores$cancer_status == "Positive"])$out > 5
outliers1_pos <- length(outliers1[outliers1 == TRUE])
outliers1 <- as.data.frame(cbind(c(rep("Negative", outliers1_neg), rep("Positive", outliers1_pos)), 
                                 c(seq(4.5, 4.5 + 0.1*(outliers1_neg - 1), 0.1))))
outliers1$V2 <- as.numeric(outliers1$V2)

plot1 <- ggplot(data_scores, aes(cancer_status, zscores1)) +
  geom_boxplot(aes(fill = cancer_status), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = cancer_status), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_point(data = outliers1, aes(V1, V2), pch = 2) +
  geom_text(data = data_stats1, aes(x = cancer_status, y = -2.7, label = N), size = 4) +
  geom_text(data = data_stats1, aes(x = cancer_status, y = 4.7, label = annot), size = 6) +
  scale_fill_manual(values = col) +
  ggtitle("90-150bp") +
  xlab("") + 
  ylab("Fragment Size Z-Score") +
  theme +
  scale_y_continuous(limits=c(-3, 5), expand = c(0,0))
plot1

outliers2 <- boxplot.stats(data_scores$zscores1[data_scores$cancer_status == "Negative"])$out > 5
outliers2_neg <- length(outliers2[outliers2 == TRUE])
outliers2 <- boxplot.stats(data_scores$zscores0[data_scores$cancer_status == "Positive"])$out > 5
outliers2_pos <- length(outliers2[outliers2 == TRUE])
outliers2 <- as.data.frame(cbind(c(rep("Negative", outliers2_neg), rep("Positive", outliers2_pos)), 
                                 c(seq(4.5, 4.5 + 0.1*(outliers2_neg - 1), 0.1))))
outliers2$V2 <- as.numeric(outliers2$V2)

plot2 <- ggplot(data_scores, aes(cancer_status, zscores2)) +
  geom_boxplot(aes(fill = cancer_status), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = cancer_status), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_point(data = outliers2, aes(V1, V2), pch = 2) +
  geom_text(data = data_stats2, aes(x = cancer_status, y = -2.7, label = N), size = 4) +
  geom_text(data = data_stats2, aes(x = cancer_status, y = 4.7, label = annot), size = 6) +
  scale_fill_manual(values = col) +
  ggtitle("151-220bp") +
  xlab("") + 
  ylab("Fragment Size Z-Score") +
  theme +
  scale_y_continuous(limits=c(-3, 5), expand = c(0,0))
plot2

outliers3 <- boxplot.stats(data_scores$zscores1[data_scores$cancer_status == "Negative"])$out > 5
outliers3_neg <- length(outliers3[outliers3 == TRUE])
outliers3 <- boxplot.stats(data_scores$zscores0[data_scores$cancer_status == "Positive"])$out > 5
outliers3_pos <- length(outliers3[outliers3 == TRUE])
outliers3 <- as.data.frame(cbind(c(rep("Negative", outliers3_neg), rep("Positive", outliers3_pos)), 
                                 c(seq(4.5, 4.5 + 0.1*(outliers3_neg - 1), 0.1))))
outliers3$V2 <- as.numeric(outliers3$V2)

plot3 <- ggplot(data_scores, aes(cancer_status, zscores3)) +
  geom_boxplot(aes(fill = cancer_status), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = cancer_status), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_point(data = outliers3, aes(V1, V2), pch = 2) +
  geom_text(data = data_stats3, aes(x = cancer_status, y = -2.7, label = N), size = 4) +
  geom_text(data = data_stats3, aes(x = cancer_status, y = 4.7, label = annot), size = 6) +
  scale_fill_manual(values = col) +
  ggtitle("221-320bp") +
  xlab("") + 
  ylab("Fragment Size Z-Score") +
  theme +
  scale_y_continuous(limits=c(-3, 5), expand = c(0,0))
plot3

figure <- ggarrange(plot0,
                    plot1 + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank()),
                    plot2 + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank()),
                    plot3 + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank()),
                    nrow = 1, align = "h", widths = c(1.4, 1, 1, 1))
figure
ggsave(file.path(outdir, "fragment_zscores.pdf"), figure,  width = 5, height = 4.5)

## Set fragment annotations
data_fragment <- as.numeric(as.matrix(row.names(data_frequency)))
data_size <- as.matrix(ifelse(data_fragment <= 89, "1",
                  ifelse(data_fragment >=90 & data_fragment <=150, "2",
                         ifelse(data_fragment >=151 & data_fragment <=220, "3", "4"))))
row.names(data_size) <- colnames(data_change)
data_size <- factor(data_size,
                    levels = c("1", "2", "3", "4"),
                    labels = c("10-89bp", "90-150bp", "151-220bp", "221-320bp"))
rm(data_fragment)

### Order patients by cancer status and score
data_scores <- data_scores[order(data_scores$cancer_status,
                                 data_scores$zscores), ]
data_change <- data_change[data_scores$sample, ]

## Set Age
data_age <- as.matrix(data_scores$Age)
row.names(data_age) <- data_scores$sample

## Set Cancer status
data_cancer <- as.matrix(data_scores$cancer_status)
row.names(data_cancer) <- data_scores$sample

## Set Previous cancers
data_history <- as.matrix(data_scores$previous_cancer)
row.names(data_history) <- data_scores$sample

## Set Tumor Stage
data_stage <- as.matrix(data_scores$stage)
row.names(data_stage) <- data_scores$sample

## Set Z-scores
data_zscores <- as.matrix(data_scores$zscores)
row.names(data_zscores) <- data_scores$sample
data_zscores <- cbind(data_zscores, c(rep(0, nrow(data_scores))), c(rep(limit, nrow(data_scores))))

data_zscores0 <- as.matrix(data_scores$zscores0)
row.names(data_zscores0) <- data_scores$sample
data_zscores0 <- cbind(data_zscores0, c(rep(0, nrow(data_scores))), c(rep(limit0, nrow(data_scores))))

data_zscores1 <- as.matrix(data_scores$zscores1)
row.names(data_zscores1) <- data_scores$sample
data_zscores1 <- cbind(data_zscores1, c(rep(0, nrow(data_scores))), c(rep(limit1, nrow(data_scores))))

data_zscores2 <- as.matrix(data_scores$zscores2)
row.names(data_zscores2) <- data_scores$sample
data_zscores2 <- cbind(data_zscores2, c(rep(0, nrow(data_scores))), c(rep(limit2, nrow(data_scores))))

data_zscores3 <- as.matrix(data_scores$zscores3)
row.names(data_zscores3) <- data_scores$sample
data_zscores3 <- cbind(data_zscores3, c(rep(0, nrow(data_scores))), c(rep(limit3, nrow(data_scores))))

## Set medians
data_median <- as.matrix(data_median$median)

## Set colours
col_fun <- colorRamp2(c(-9, -4.5, 0, 4.5, 9), 
                      c("#1f78b4", "#a6cee3", "white", "#fb9a99", "#e31a1c"))
col_cancer <- c("Positive" = "#fb9a99", "Negative" = "#a6cee3")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_stage <- c("Stage III/IV" = "#fb9a99", "Stage 0/I/II" = "#ebdc78", "none" = "grey95", Unknown = "grey65")

## Set additional annotations
left_annotation <- rowAnnotation("Patient Type" = data_age,
                                 "Cancer History" = data_history,
                                 border = FALSE,
                                 annotation_name_side = "top",
                                 show_annotation_name = FALSE,
                                 col = list("Patient Type" = col_age,
                                            "Cancer History" = col_history,
                                            "Cancer Status" = col_cancer),
                                 show_legend = FALSE,
                                 simple_anno_size = unit(0.3, "cm"))

right_annotation <- rowAnnotation("   " = data_stage,
                                  " " = data_cancer,
                                  "10-89bp" = anno_lines(data_zscores0,
                                                         add_points = TRUE,
                                                         pch = c(16, NA, NA),
                                                         gp = gpar(col = c(NA, "black", "red"),
                                                                   lty = c("solid", "dashed", "dashed")),
                                                         pt_gp = gpar(col = c("black", NA, NA)),
                                                         size = unit(1, "mm"),
                                                         axis_param = list(side = "top",
                                                                           labels_rot = 0,
                                                                           gp = gpar(fontsize = 8)),
                                                         border = FALSE),
                                  "90-150bp" = anno_lines(data_zscores1,
                                                          add_points = TRUE,
                                                          pch = c(16, NA, NA),
                                                          gp = gpar(col = c(NA, "black", "red"),
                                                                    lty = c("solid", "dashed", "dashed")),
                                                          pt_gp = gpar(col = c("black", NA, NA)),
                                                          size = unit(1, "mm"),
                                                          axis_param = list(side = "top",
                                                                            labels_rot = 0,
                                                                            gp = gpar(fontsize = 8)),
                                                          border = FALSE),
                                  "151-220bp" = anno_lines(data_zscores2,
                                                           add_points = TRUE,
                                                           pch = c(16, NA, NA),
                                                           gp = gpar(col = c(NA, "black", "red"),
                                                                     lty = c("solid", "dashed", "dashed")),
                                                           pt_gp = gpar(col = c("black", NA, NA)),
                                                           size = unit(1, "mm"),
                                                           axis_param = list(side = "top",
                                                                             labels_rot = 0,
                                                                             gp = gpar(fontsize = 8)),
                                                           border = FALSE),
                                  "221-320bp" = anno_lines(data_zscores3,
                                                           add_points = TRUE,
                                                           pch = c(16, NA, NA),
                                                           gp = gpar(col = c(NA, "black", "red"),
                                                                     lty = c("solid", "dashed", "dashed")),
                                                           pt_gp = gpar(col = c("black", NA, NA)),
                                                           size = unit(1, "mm"),
                                                           axis_param = list(side = "top",
                                                                             labels_rot = 0,
                                                                             gp = gpar(fontsize = 8)),
                                                           border = FALSE),
                                  "Integrated" = anno_lines(data_zscores,
                                                            add_points = TRUE,
                                                            pch = c(16, NA, NA),
                                                            gp = gpar(col = c(NA, "black", "red"),
                                                                      lty = c("solid", "dashed", "dashed")),
                                                            pt_gp = gpar(col = c("black", NA, NA)),
                                                            size = unit(1, "mm"),
                                                            axis_param = list(side = "top",
                                                                              labels_rot = 0,
                                                                              gp = gpar(fontsize = 8)),
                                                            border = FALSE),
                                  gap = unit(2, "mm"),
                                  annotation_name_side = "top",
                                  annotation_name_gp = gpar(fontsize = 8),
                                  annotation_name_rot = 0,
                                  col = list(" " = col_cancer,
                                             "   " = col_stage),
                                  show_legend = FALSE,
                                  simple_anno_size = unit(0.3, "cm"),
                                  width = unit(9, "cm"))

freq_annotation <- HeatmapAnnotation("LFS Healthy\nMedian" = anno_lines(data_median,
                                                                     gp = gpar(col = c("black", "red"),
                                                                               fontsize = 10),
                                                                     axis = FALSE),
                                     border = FALSE,
                                     annotation_name_side = "left",
                                     annotation_name_gp = gpar(fontsize = 8),
                                     height = unit(2, "cm"),
                                     annotation_name_rot = 0)

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Patient Type", 
                                                  at = c("Adult", "Pediatric"),
                                                  legend_gp = gpar(fill = c("#6A3D9A", "#CAB2D6")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer History",
                                                  at = c("Yes"),
                                                  legend_gp = gpar(fill = c("#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Tumor Stage",
                                                  at = c("Stage 0/I/II", "Stage III/IV", "Unknown/None"),
                                                  legend_gp = gpar(fill = c("#ebdc78", "#fb9a99", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Status",
                                                  at = c("Positive", "Negative"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "#a6cee3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "SD from LFS Healthy",
                                                  at = c(-9, 0, 9),
                                                  col_fun = col_fun,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm")),
                                           Legend(title = "Z-score",
                                                  type = "lines",
                                                  at = c("LFS Healthy Median", "95th Percentile"),
                                                  legend_gp = gpar(col = c("black", "red"), 
                                                                   lty = c("dashed", "dashed")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white")),
                               direction = "horizontal")

## Set orders
row_order <- row.names(data_change)
col_order <- colnames(data_change)
#row_split <- data_scores$sample_parent
#order <- unique(row_split)
#row_split <- factor(row_split, levels = order)
col_split <- data_size

## Generate heatmap
pdf(file.path(outdir, "fragment_frequency.pdf"), width = 7, height = 10)
#pdf(file.path(outdir, "fragment_frequency_short.pdf"), width = 7, height = 6)
Heatmap <- Heatmap(data_change,
                   col = col_fun,
                   show_heatmap_legend = FALSE,
                   row_order = row_order,
                   column_order = col_order,
                   left_annotation = left_annotation,
                   right_annotation = right_annotation,
                   top_annotation = freq_annotation,
                   #row_split = row_split,
                   column_split = col_split,
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 8),
                   show_row_names = FALSE,
                   column_title_gp = gpar(fontsize = 8),
                   show_column_names = FALSE,
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom") 
dev.off()

data_scores$limit <- limit
data_scores$limit0 <- limit0
data_scores$limit1 <- limit1
data_scores$limit2 <- limit2
data_scores$limit3 <- limit3
write.table(data_scores, file.path(outdir, "fragment_scores.txt"), sep = "\t", row.names = FALSE)

### Calculate stats for regression
regression <- summary(model)[["coefficients"]]
regression <- regression[-1, ]
regression <- as.data.frame(regression)
colnames(regression) <- c("estimate", "se", "zvalue", "pvalue")
regression$term <- row.names(regression)
regression$term <- factor(regression$term, levels = c("zscores3", "zscores2", "zscores1", "zscores0"),
                          labels = c("221-320bp", "151-220bp", "90-150bp", "10-89bp"))
regression$pvalue <- round(regression$pvalue, 3)
regression$upper <- regression$estimate + 1.96*regression$se
regression$lower <- regression$estimate - 1.96*regression$se

### Plot regression results
plot_score <- ggplot(data_scores, aes(x = cancer_status, y = zscores)) +
  geom_boxplot(aes(fill = cancer_status), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = cancer_status), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = data_stats, aes(x = cancer_status, y = -3, label = N), size = 4) +
  geom_text(data = data_stats, aes(x = cancer_status, y = 9, label = annot), size = 5) +
  xlab("Cancer Status") + 
  ylab("Integrated Z-score") +
  ggtitle("Combined Metric") + 
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  theme(plot.title = element_text(hjust = 0.5, size = 13), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(limits=c(-4, 10), expand = c(0,0))
plot_score

plot_estimates <- ggplot(regression, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_text(aes(x = term, y = 2, label = pvalue), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Compartment") + 
  ylab("Estimate") +
  ggtitle("Regression Estimates") + 
  theme(plot.title = element_text(hjust = 0.5, size = 13), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 0)) +
  scale_y_continuous(limits=c(-1.5, 2.5), expand = c(0,0)) +
  coord_flip()
plot_estimates

figure_reg <- ggarrange(plot_score, plot_estimates, align = "h", nrow = 1, widths = c(1, 1.75))
figure_reg

ggsave(file.path(outdir, "fragment_regression.pdf"), figure_reg, device = "pdf", width = 5, height = 4.5, units = "in")

