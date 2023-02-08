library(tidyverse)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/fragment_score"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_score"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/fragment_score"
ichor <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA/CHARM_LFS_ichorCNA_summary_reviewed.txt"
mutations <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot/Oncoplot_full.txt"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

### Read in data
score_wg <- read.delim(list.files(path, "genome", full.names = TRUE))
score_ts <- read.delim(list.files(path, "panel", full.names = TRUE))
data_mutation <- read.delim(mutations)
data_ichor <- read.delim(ichor)
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
failed_TS <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
failed_WG <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
exclude <- c("TGL49_0316_Pl_T_PE_301_CM", "TGL49_0060_Cf_n_PE_301_CM")
data_samples <- data_samples[!(data_samples$sWGS %in% failed_WG), ]
data_samples <- data_samples[!(data_samples$TS %in% failed_TS), ]
data_samples <- data_samples[!(data_samples$cfMeDIP %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% score_wg$sample, ]
score_wg <- score_wg[score_wg$sample %in% data_samples$sWGS, ]
score_ts <- score_ts[score_ts$sample %in% data_samples$TS, ]
data_ichor <- data_ichor[data_ichor$sample %in% data_samples$sWGS, ]
data_mutation <- data_mutation[data_mutation$sample %in% data_samples$TS, ]

### Merge all data
data_ichor$ichor <- ifelse(data_ichor$TF > data_ichor$TF_short, data_ichor$TF, data_ichor$TF_short)
data_mutation <- data_mutation[, c("sample_ID", "TP53_somatic", "vaf")]
colnames(score_ts) <- c("sample", "tp53")
colnames(score_wg) <- c("sample", "genome")

data <- merge(data_samples, data_ichor[,c("sample", "ichor")], by.x = "sWGS", by.y = "sample", all = TRUE)
data <- merge(data, data_mutation, by.x = "TS", by.y = "sample_ID", all = TRUE)
data <- merge(data, score_ts, by.x = "TS", "sample", all = TRUE)
data <- merge(data, score_wg, by.x = "sWGS", by.y = "sample", all = TRUE)

### Find score limits
limit_ts <- quantile(data$tp53[data$cancer_status == "negative" & data$previous_cancer == "no"], 0.95, na.rm = TRUE)
limit_wg <- quantile(data$genome[data$cancer_status == "negative" & data$previous_cancer == "no"], 0.95, na.rm = TRUE)
data$tp53 <- ifelse(data$tp53 > limit_ts, "pos", "neg")
data$tp53[is.na(data$tp53)] <- "neg"
data$genome <- ifelse(data$genome > limit_wg, "pos", "neg")
data$ichor <- ifelse(data$ichor > 0.03, "pos", "neg")
data$TP53_somatic <- ifelse(data$TP53_somatic == "", "neg", "pos")
data$TP53_somatic[is.na(data$TP53_somatic)] <- "neg"

### Stats for positive samples
data_pos <- data[data$cancer_status == "positive", ]
pos <- nrow(data_pos)
pos_ts <- nrow(data_pos[data_pos$tp53 == "pos", ])
pos_ts_only <- nrow(data_pos[data_pos$TP53_somatic %in% c("neg", NA) &
                              data_pos$ichor == "neg" &
                              data_pos$tp53 == "pos", ])

pos_wg <- nrow(data_pos[data_pos$genome == "pos", ])
pos_wg_only <- nrow(data_pos[data_pos$TP53_somatic %in% c("neg", NA) &
                               data_pos$ichor == "neg" &
                               data_pos$genome == "pos", ])

pos_frag_all <- nrow(data_pos[data_pos$TP53_somatic %in% c("neg", NA) &
                                data_pos$ichor == "neg" &
                                (data_pos$genome == "pos" | data_pos$tp53 == "pos"), ])

pos_frag_overlap <- nrow(data_pos[data_pos$TP53_somatic %in% c("neg", NA) &
                                    data_pos$ichor == "neg" &
                                    data_pos$genome == "pos" &
                                    data_pos$tp53 == "pos", ])

### Stats for negative samples
data_neg <- data[data$cancer_status == "negative", ]
neg <- nrow(data_neg)
neg_ts <- nrow(data_neg[data_neg$tp53 == "pos", ])
neg_ts_only <- nrow(data_neg[data_neg$TP53_somatic %in% c("neg", NA) &
                               data_neg$ichor == "neg" &
                               data_neg$tp53 == "pos", ])

neg_wg <- nrow(data_neg[data_neg$genome == "pos", ])
neg_wg_only <- nrow(data_neg[data_neg$TP53_somatic %in% c("neg", NA) &
                               data_neg$ichor == "neg" &
                               data_neg$genome == "pos", ])

neg_frag_all <- nrow(data_neg[data_neg$TP53_somatic %in% c("neg", NA) &
                                data_neg$ichor == "neg" &
                                (data_neg$genome == "pos" | data_neg$tp53 == "pos"), ])

neg_frag_overlap <- nrow(data_neg[data_neg$TP53_somatic %in% c("neg", NA) &
                                    data_neg$ichor == "neg" &
                                    data_neg$genome == "pos" &
                                    data_neg$tp53 == "pos", ])

### Summarize stats
descriptors <- c("positive samples", "TP53", "Genome", "TP53 Only", "Genome Only", "Frag Only", "Overlap")
stats <- c(pos, pos_ts, pos_wg, pos_ts_only, pos_wg_only, pos_frag_all, pos_frag_overlap)
stats_pos <- as.data.frame(cbind(descriptors, stats))

descriptors <- c("negative samples", "TP53", "Genome", "TP53 Only", "Genome Only", "Frag Only", "Overlap")
stats <- c(neg, neg_ts, neg_wg, neg_ts_only, neg_wg_only, neg_frag_all, neg_frag_overlap)
stats_neg <- as.data.frame(cbind(descriptors, stats))

