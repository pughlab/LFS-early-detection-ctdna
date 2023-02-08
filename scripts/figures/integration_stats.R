library(dplyr)
library(ComplexHeatmap)
library(UpSetR)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
path2 <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/integration"

### Find files
ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", recursive = TRUE, full.names = TRUE))
mutation <- read.delim(file.path(path2, "oncoplot/Oncoplot_full.txt"))
tp53 <- read.delim(file.path(path, "fragment_score/CHARM_LFS_panel_score.txt"))
fragment <- read.delim(file.path(path, "fragment_score/CHARM_LFS_genome_score.txt"))
breast <- read.delim(file.path(path2, "cfMeDIP/breast_matrix.tsv"))
all <- read.delim(file.path(path2, "cfMeDIP/pancancer_score_Vrba.tsv"))
samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Calculate cumulative breast score
breast <- data.frame(One = colSums(breast),
                     cfMeDIP = colnames(breast))

### Calculate methylation limits
breast_limit <- quantile(breast$One[!(breast$cfMeDIP %in% samples$cfMeDIP)], 0.99)
all_limit <- quantile(all$score[!(all$cfMeDIP %in% samples$cfMeDIP)], 0.99)
tp53_limit <- quantile(tp53$score[tp53$sample %in% samples$TS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)
fragment_limit <- quantile(fragment$score[fragment$sample %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)

### Remove failed and unknown samples and format 
exclude_wg <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
exclude_ts <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
exclude_cm <- c("TGL49_0010_Cf_U_PE_319_CM", "TGL49_0035_Cf_U_PE_313_CM", "TGL49_0041_Cf_U_PE_311_CM", "TGL49_0209_Cf_U_PE_355_CM")

samples <- samples[!(samples$sWGS %in% exclude_wg &
                       samples$TS %in% exclude_ts &
                       samples$cfMeDIP %in% exclude_cm), ]
ichorCNA <- ichorCNA[!(ichorCNA$sample %in% exclude_wg), ]
fragment <- fragment[!(fragment$sample %in% exclude_wg), ]
mutation <- mutation[!(mutation$sample_ID %in% exclude_ts), ]
tp53 <- tp53[!(tp53$sample %in% exclude_ts), ]
breast <- breast[!(breast$cfMeDIP %in% exclude_cm), ]
all <- all[!(all$cfMeDIP %in% exclude_cm), ]

### Merge data together
colnames(fragment) <- c("sample", "genome")
colnames(tp53) <- c("sample", "tp53")
ichorCNA$ichor <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA <- ichorCNA[, c("sample", "ichor")]
mutation <- mutation[, c("sample_ID", "TP53_somatic")]
breast <- breast[, c("cfMeDIP", "One")]
colnames(breast) <- c("cfMeDIP", "breast")
all <- all[, c("cfMeDIP", "score")]
colnames(all) <- c("cfMeDIP", "pancancer")

samples <- samples[!(samples$sWGS == "" & samples$TS == "" & samples$cfMeDIP == ""), ]
samples <- merge(samples, fragment, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, ichorCNA, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, mutation, by.x = "TS", by.y = "sample_ID", all = TRUE)
samples <- merge(samples, tp53, by.x = "TS", by.y = "sample", all = TRUE)
samples <- merge(samples, breast, by = "cfMeDIP", all = TRUE)
samples <- merge(samples, all, by = "cfMeDIP", all = TRUE)
samples <- samples[!(is.na(samples$sWGS)), ]
samples <- samples[rowSums(is.na(samples[, c(26:31)])) != 6, ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

### Calculate number of samples with all 3
samples_pos <- samples[samples$cancer_status == "positive", ]
samples_pos_ts <- nrow(samples_pos[!(samples_pos$TS == ""), ])
samples_pos_wg <- nrow(samples_pos[!(samples_pos$sWGS == ""), ])
samples_pos_cm <- nrow(samples_pos[!(samples_pos$cfMeDIP == ""), ])
samples_pos_all <- nrow(samples_pos[!(samples_pos$TS == "" |
                                       samples_pos$sWGS == "" |
                                       samples_pos$cfMeDIP == ""), ])
samples_pos_any <- nrow(samples_pos)

samples_neg <- samples[samples$cancer_status == "negative", ]
samples_neg_ts <- nrow(samples_neg[!(samples_neg$TS == ""), ])
samples_neg_wg <- nrow(samples_neg[!(samples_neg$sWGS == ""), ])
samples_neg_cm <- nrow(samples_neg[!(samples_neg$cfMeDIP == ""), ])
samples_neg_all <- nrow(samples_neg[!(samples_neg$TS == "" |
                                        samples_neg$sWGS == "" |
                                        samples_neg$cfMeDIP == ""), ])
samples_neg_any <- nrow(samples_neg)

### Format data
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

### Format detection
samples$TP53_somatic <- ifelse(samples$TP53_somatic == "", "neg",
                               ifelse(is.na(samples$TP53_somatic), NA, "pos"))
samples$ichor <- ifelse(samples$ichor >= 0.03, "pos", "neg")
samples$genome <- ifelse(samples$genome > fragment_limit, "pos", "neg")
samples$tp53 <- ifelse(samples$tp53 > tp53_limit, "pos", "neg")
samples$breast <- ifelse(samples$breast > breast_limit, "pos", "neg")
samples$pancancer <- ifelse(samples$pancancer > all_limit, "pos", "neg")
m <- c("genome", "tp53", "ichor", "TP53_somatic", "breast", "pancancer")
samples <- samples[rowSums(is.na(samples[, colnames(samples) %in% m])) != length(m), ]
samples[is.na(samples)] <- "NS"

### Stats for cancer positive patients
data_pos <- samples[samples$cancer_status == "positive", ]
pos_samples <- nrow(data_pos)

pos_mut <- nrow(data_pos[data_pos$TP53_somatic == "pos", ])

pos_ichor <- nrow(data_pos[data_pos$ichor == "pos", ])
pos_fragment <- nrow(data_pos[data_pos$genome == "pos", ])
pos_tp53 <- nrow(data_pos[data_pos$tp53 == "pos", ])

pos_breast <- nrow(data_pos[data_pos$breast == "pos",])
pos_all <- nrow(data_pos[data_pos$pancancer == "pos",])

pos_any <- nrow(data_pos[data_pos$TP53_somatic == "pos" |
                           data_pos$ichor == "pos" |
                           data_pos$genome == "pos" |
                           data_pos$breast == "pos" |
                           data_pos$pancancer == "pos", ])

### Stats for high stage patients
data_high <- data_pos[data_pos$stage == "high", ]
high_samples <- nrow(data_high)

high_ts <- nrow(data_high[!(data_high$TS %in% c("NS", NA)), ])
high_mut <- nrow(data_high[data_high$TP53_somatic == "pos", ])
high_tp53 <- nrow(data_high[data_high$tp53 == "pos", ])

high_wg <- nrow(data_high[!(data_high$sWGS %in% c("NS", NA)), ])
high_ichor <- nrow(data_high[data_high$ichor == "pos", ])
high_fragment <- nrow(data_high[data_high$genome == "pos", ])

high_cm <- nrow(data_high[!(data_high$cfMeDIP %in% c("NS", NA)), ])
high_breast <- nrow(data_high[data_high$breast == "pos", ])
high_all <- nrow(data_high[data_high$pancancer == "pos", ])

high_any <- nrow(data_high[data_high$TP53_somatic == "pos" |
                             data_high$ichor == "pos" |
                             data_high$genome == "pos" |
                             data_high$breast == "pos" |
                             data_high$pancancer == "pos", ])

### Stats for low stage patients
data_low <- data_pos[data_pos$stage == "low", ]
low_samples <- nrow(data_low)

low_ts <- nrow(data_low[!(data_low$TS %in% c("NS", NA)), ])
low_mut <- nrow(data_low[data_low$TP53_somatic == "pos", ])
low_tp53 <- nrow(data_low[data_low$tp53 == "pos", ])

low_wg <- nrow(data_low[!(data_low$sWGS %in% c("NS", NA)), ])
low_ichor <- nrow(data_low[data_low$ichor == "pos", ])
low_fragment <- nrow(data_low[data_low$genome == "pos", ])

low_cm <- nrow(data_low[!(data_low$cfMeDIP %in% c("NS", NA))])
low_breast <- nrow(data_low[data_low$breast == "pos", ])
low_all <- nrow(data_low[data_low$pancancer == "pos", ])

low_any <- nrow(data_low[data_low$TP53_somatic == "pos" |
                             data_low$ichor == "pos" |
                             data_low$genome == "pos" |
                             data_low$breast == "pos" |
                             data_low$pancancer == "pos", ])

### Stats for breast
data_breast <- data_pos[data_pos$cancer_type == "breast", ]
breast_samples <- nrow(data_breast)

breast_ts <- nrow(data_breast[!(data_breast$TS %in% c("NS", NA)), ])
breast_mut <- nrow(data_breast[data_breast$TP53_somatic == "pos", ])
breast_tp53 <- nrow(data_breast[data_breast$tp53 == "pos", ])

breast_wg <- nrow(data_breast[!(data_breast$sWGS %in% c("NS", NA)), ])
breast_ichor <- nrow(data_breast[data_breast$ichor == "pos", ])
breast_fragment <- nrow(data_breast[data_breast$genome == "pos", ])

breast_cm <- nrow(data_breast[!(data_breast$cfMeDIP %in% c("NS", NA)), ])
breast_breast <- nrow(data_breast[data_breast$breast == "pos", ])
breast_all <- nrow(data_breast[data_breast$pancancer == "pos", ])

breast_cm_only <- nrow(data_breast[data_breast$breast == "pos" &
                                     data_breast$TP53_somatic == "neg" &
                                     data_breast$ichor == "neg" &
                                     data_breast$genome == "neg", ])

breast_any <- nrow(data_breast[data_breast$TP53_somatic == "pos" |
                               data_breast$ichor == "pos" |
                               data_breast$genome == "pos" |
                               data_breast$breast == "pos" |
                               data_breast$pancancer == "pos", ])

### Stats for false positives
data_neg <- samples[samples$cancer_status == "negative", ]

neg_samples <- nrow(data_neg)
neg_patients <- length(unique(data_neg$ext_ID))

neg_ts <- nrow(data_neg[!(data_neg$TS %in% c("", NA)), ])
neg_mut <- nrow(data_neg[!(data_neg$TP53_somatic %in% c("neg", "NS")), ])
neg_tp53 <- nrow(data_neg[data_neg$tp53 == "pos", ])

neg_wg <- nrow(data_neg[!(data_neg$sWGS %in% c("", NA)), ])
neg_ichor <- nrow(data_neg[!(data_neg$ichor %in% c("neg", "NS")), ])
neg_fragment <- nrow(data_neg[!(data_neg$genome %in% c("neg", "NS")), ])

neg_cm <- nrow(data_neg[!(data_neg$cfMeDIP %in% c("", NA)), ])
neg_breast <- nrow(data_neg[!(data_neg$breast %in% c("neg", "NS")), ])
neg_all <- nrow(data_neg[!(data_neg$pancancer %in% c("neg", "NS")), ])

false_positives <- data_neg[data_neg$TP53_somatic == "pos" |
                              data_neg$tp53 == "pos" |
                              data_neg$ichor == "pos" |
                              data_neg$genome == "pos" |
                              data_neg$breast == "pos" |
                              data_neg$pancancer == "pos", ]
false_positives <- false_positives[complete.cases(false_positives[,1:3]), ]
false_samples <- nrow(false_positives)
false_patients <- length(unique(false_positives$ext_ID))
write.table(false_positives, file.path(outdir, "false_positives.tsv"), sep = "\t", row.names = FALSE)

all_negative <- data_neg[!(data_neg$sWGS %in% false_positives$sWGS), ]
all_negative <- all_negative[!(all_negative$ext_ID %in% false_positives$ext_ID), ]
all_negative_samples <- nrow(all_negative)
all_negative_patients <- length(unique(all_negative$ext_ID))
write.table(all_negative, file.path(outdir, "all_negatives.tsv"), sep = "\t", row.names = FALSE)

### Make table of stats
descriptors <- c("positive samples", "TS", "WG", "CM", "All",  "mutation positive", "TP53 frag pos", "ichorCNA positive", "fragment positive", "breast positive", "pancan positive", "any positive")
stats <- c(pos_samples, samples_pos_ts, samples_pos_wg, samples_pos_cm, samples_pos_all, pos_mut, pos_tp53, pos_ichor, pos_fragment, pos_breast, pos_all, pos_any)
stats_positive <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, pos_samples, pos_samples, pos_samples, pos_samples, samples_pos_ts, samples_pos_ts,samples_pos_wg, samples_pos_wg, samples_pos_cm, samples_pos_cm, pos_samples)
stats_positive$percent <- as.numeric(stats_positive$stats)/denominators*100

descriptors <- c("late stage samples", "TS", "WG", "CM", "mutation positive", "TP53 frag pos", "ichorCNA positive", "fragment positive", "breast positive", "pancan positive", "any positive")
stats <- c(high_samples, high_ts, high_wg, high_cm, high_mut, high_tp53, high_ichor, high_fragment, high_breast, high_all, high_any)
stats_high <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, high_samples, high_samples, high_samples, high_ts, high_ts, high_wg, high_wg, high_cm, high_cm, high_samples)
stats_high$percent <- as.numeric(stats_high$stats)/denominators*100

descriptors <- c("early stage samples", "TS", "WG", "CM", "mutation positive", "TP53 frag pos", "ichorCNA positive", "fragment positive", "breast positive", "pancan positive", "any positive")
stats <- c(low_samples, low_ts, low_wg, low_cm, low_mut, low_tp53, low_ichor, low_fragment, low_breast, low_all, low_any)
stats_low <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, low_samples, low_samples, low_samples, low_ts, low_ts, low_wg, low_wg, low_cm, low_cm, low_samples)
stats_low$percent <- as.numeric(stats_low$stats)/denominators*100

descriptors <- c("Breast cancers", "TS", "WG", "CM", "mutation positive", "TP53 frag pos", "ichorCNA positive", "fragment positive", "breast positive", "pancan positive", "any positive", "methylation only")
stats <- c(breast_samples, breast_ts, breast_wg, breast_cm, breast_mut, breast_tp53, breast_ichor, breast_fragment, breast_breast, breast_all, breast_any, breast_cm_only)
stats_breast <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, breast_samples, breast_samples, breast_samples, breast_ts, breast_ts, breast_wg, breast_wg, breast_cm, breast_cm, breast_samples, breast_samples)
stats_breast$percent <- as.numeric(stats_breast$stats)/denominators*100

descriptors <- c("neg samples", "neg patients", "TS", "WG", "CM", "All", "mutation positive", "TP53 frag pos", "ichorCNA positive", "fragment positive", "breast positive", "pancan positive", 
                 "false positives", "false positive patients", "all negative samples", "all negative patients")
stats <- c(neg_samples, neg_patients, samples_neg_ts, samples_neg_wg, samples_neg_cm, samples_neg_all, neg_mut, neg_tp53, neg_ichor, neg_fragment, neg_breast, neg_all, 
           false_samples, false_patients, all_negative_samples, all_negative_patients)
stats_neg <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, neg_samples, neg_samples, neg_samples, neg_samples, neg_samples, neg_ts, neg_ts, neg_wg, neg_wg, neg_cm, neg_cm, neg_samples, NA, neg_samples, NA)
stats_neg$percent <- as.numeric(stats_neg$stats)/denominators*100

### Make table for Venn Diagram (all positive)
data_pos$sWGS <- ifelse(data_pos$sWGS == "", data_pos$TS, data_pos$sWGS)
data_pos$sWGS <- ifelse(data_pos$sWGS == "", data_pos$cfMeDIP, data_pos$sWGS)
data_venn <- data_pos[, colnames(data_pos) %in% c("sWGS", "genome", "tp53", "ichor", "TP53_somatic", "breast", "pancancer")]
colnames(data_venn) <- c("sample", "Fragment Score (Genome)", "Fragment Score (TP53)", "Genome-wide CNV", "TP53 Alterations", "Breast Methylation", "Pan-cancer Methylation")
data_venn[data_venn == "pos"] <- 1
data_venn[data_venn == "neg"] <- 0
data_venn[data_venn == "NS"] <- 0
data_venn[is.na(data_venn)] <- 0
data_venn[, c(2:ncol(data_venn))] <- sapply(data_venn[, 2:ncol(data_venn)], as.numeric)
data_venn_high <- data_venn[data_venn$sample %in% data_pos$sWGS[data_pos$stage == "high"], ]
data_venn_low <- data_venn[data_venn$sample %in% data_pos$sWGS[data_pos$stage == "low"], ]
data_venn_breast <- data_venn[data_venn$sample %in% data_pos$sWGS[data_pos$cancer_type == "breast"], ]

data_neg$sWGS <- ifelse(data_neg$sWGS == "", data_neg$TS, data_neg$sWGS)
data_neg$sWGS <- ifelse(data_neg$sWGS == "", data_neg$cfMeDIP, data_neg$sWGS)
data_venn_neg <- data_neg[, colnames(data_neg) %in% c("sWGS", "genome", "tp53", "ichor", "TP53_somatic", "breast", "pancancer")]
colnames(data_venn_neg) <- c("sample", "Fragment Score (Genome)", "Fragment Score (TP53)", "Genome-wide CNV", "TP53 Alterations", "Breast Methylation", "Pan-cancer Methylation")
data_venn_neg[data_venn_neg == "pos"] <- 1
data_venn_neg[data_venn_neg == "neg"] <- 0
data_venn_neg[data_venn_neg == "NS"] <- 0
data_venn_neg[is.na(data_venn_neg)] <- 0
data_venn_neg[, c(2:ncol(data_venn_neg))] <- sapply(data_venn_neg[, 2:ncol(data_venn_neg)], as.numeric)

### Make upset plots
## All Cancers
m1 <- make_comb_mat(data_venn)
cs <- comb_size(m1)
upset <- make_comb_mat(data_venn[,-1])

pdf(file.path(outdir, "integration_positive_upset.pdf"), width = 5.5, height = 2.75)
UpSet <- UpSet(upset, 
               top_annotation = HeatmapAnnotation(
                 "# Samples" = anno_barplot(cs,
                                          ylim = c(0, max(cs)*1.1),
                                          border = FALSE, 
                                          gp = gpar(fill = "black"),
                                          height = unit(3, "cm")), 
                 annotation_name_side = "left", 
                 annotation_name_rot = 90),
               comb_order = order(comb_size(m1), decreasing = TRUE))
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("# Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})
grid.text("All Cancers", x = 0.55, y = 0.95, just = "left",
          gp = gpar(fontsize = 12))
dev.off()

## late stage
m1 <- make_comb_mat(data_venn_high[,-1])
cs <- comb_size(m1)
upset <- make_comb_mat(data_venn_high[,-1])

pdf(file.path(outdir, "integration_positive_high_upset.pdf"), width = 5.5, height = 2.75)
UpSet <- UpSet(upset, 
               top_annotation = HeatmapAnnotation(
                 "# Samples" = anno_barplot(cs,
                                            ylim = c(0, max(cs)*1.1),
                                            border = FALSE, 
                                            gp = gpar(fill = "black"),
                                            height = unit(3, "cm")), 
                 annotation_name_side = "left", 
                 annotation_name_rot = 90),
               comb_order = order(comb_size(m1), decreasing = TRUE))
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("# Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})
grid.text("Late Stage Cancers", x = 0.55, y = 0.95, just = "left",
          gp = gpar(fontsize = 12))
dev.off()

## early stage
m1 <- make_comb_mat(data_venn_low[,-1])
cs <- comb_size(m1)
upset <- make_comb_mat(data_venn_low[,-1])

pdf(file.path(outdir, "integration_positive_low_upset.pdf"), width = 4.5, height = 2.75)
UpSet <- UpSet(upset, 
               top_annotation = HeatmapAnnotation(
                 "# Samples" = anno_barplot(cs,
                                            ylim = c(0, max(cs)*1.1),
                                            border = FALSE, 
                                            gp = gpar(fill = "black"),
                                            height = unit(3, "cm")), 
                 annotation_name_side = "left", 
                 annotation_name_rot = 90),
               comb_order = order(comb_size(m1), decreasing = TRUE))
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("# Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})
grid.text("Early Stage\nCancers", x = 0.8, y = 0.925, just = "center",
          gp = gpar(fontsize = 12))
dev.off()

## Breast Cancers
m1 <- make_comb_mat(data_venn_breast[,-1])
cs <- comb_size(m1)
upset <- make_comb_mat(data_venn_breast[,-1])

pdf(file.path(outdir, "integration_breast_upset.pdf"), width = 5.25, height = 2.75)
UpSet <- UpSet(upset, 
               top_annotation = HeatmapAnnotation(
                 "# Samples" = anno_barplot(cs,
                                            ylim = c(0, max(cs)*1.1),
                                            border = FALSE, 
                                            gp = gpar(fill = "black"),
                                            height = unit(3, "cm")), 
                 annotation_name_side = "left", 
                 annotation_name_rot = 90),
               comb_order = order(comb_size(m1), decreasing = TRUE))
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("# Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})
grid.text("Breast Cancers", x = 0.55, y = 0.95, just = "left",
          gp = gpar(fontsize = 12))
dev.off()

## False positives
m1 <- make_comb_mat(data_venn_neg[,-1])
cs <- comb_size(m1)
upset <- make_comb_mat(data_venn_neg[,-1])

pdf(file.path(outdir, "integration_negatives_upset.pdf"), width = 5.25, height = 2.75)
UpSet <- UpSet(upset, 
               top_annotation = HeatmapAnnotation(
                 "# Samples" = anno_barplot(cs,
                                            ylim = c(0, max(cs)*1.1),
                                            border = FALSE, 
                                            gp = gpar(fill = "black"),
                                            height = unit(3, "cm")), 
                 annotation_name_side = "left", 
                 annotation_name_rot = 90),
               comb_order = order(comb_size(m1), decreasing = TRUE))
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("# Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})
grid.text("Cancer Negative", x = 0.55, y = 0.95, just = "left",
          gp = gpar(fontsize = 12))
dev.off()

