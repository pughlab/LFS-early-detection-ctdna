library(tidyverse)

### Set variables
path <- ""
outdir <- ""

ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", full.names = TRUE))
mutations <- read.delim("Oncoplot_full.txt")
samples <- read.delim("sample_list.txt")

### Remove failed data
exclude <- c("TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG", "TGL49_0316_Pl_T_PE_301_CM", "TGL49_0060_Cf_n_PE_301_CM")
failed <- "TGL49_0025_Cf_U_PE_321_WG"
samples <- samples[!(samples$sWGS %in% c(exclude, failed)), ]
samples <- samples[samples$sWGS %in% ichorCNA$sample, ]
samples <- samples[!(samples$cfMeDIP %in% exclude), ]

### Format ichorCNA data
ichorCNA <- ichorCNA[, c("sample", "TF", "TF_short")]
ichorCNA <- merge(samples, ichorCNA, by.x = "sWGS", by.y = "sample")

### Calculate general statistics
cohort <- nrow(ichorCNA)
patients <- length(unique(ichorCNA$ext_ID))

adult <- ichorCNA[ichorCNA$Age == "adult", ]
adult_samples <- nrow(adult)
adult_patients <- length(unique(adult$ext_ID))

pediatric <- ichorCNA[ichorCNA$Age == "pediatric", ]
pediatric_samples <- nrow(pediatric)
pediatric_patients <- length(unique(pediatric$ext_ID))

full_positive <- nrow(ichorCNA[ichorCNA$TF > 0.03, ])
full_positive_pos <- nrow(ichorCNA[ichorCNA$TF > 0.03 & ichorCNA$cancer_status == "positive", ])
full_positive_neg <- nrow(ichorCNA[ichorCNA$TF > 0.03 & ichorCNA$cancer_status == "negative", ])

### full vs short
concordant <- ichorCNA[(ichorCNA$TF > 0.03 & ichorCNA$TF_short > 0.03) |
                        (ichorCNA$TF < 0.03 & ichorCNA$TF_short < 0.03), ]
concordant_neg <- nrow(concordant[concordant$cancer_status == "negative", ])
concordant_pos <- nrow(concordant[concordant$cancer_status == "positive", ])
concordant <- nrow(concordant)

discordant <- ichorCNA[(ichorCNA$TF > 0.03 & ichorCNA$TF_short < 0.03) |
                         (ichorCNA$TF < 0.03 & ichorCNA$TF_short > 0.03), ]
discordant_neg <- nrow(discordant[discordant$cancer_status == "negative", ])
discordant_pos <- nrow(discordant[discordant$cancer_status == "positive", ])

full_only <- nrow(discordant[discordant$TF > 0.03, ])
full_pos <- nrow(discordant[discordant$TF > 0.03 & discordant$cancer_status == "positive", ])
full_neg <- nrow(discordant[discordant$TF > 0.03 & discordant$cancer_status == "negative", ])

short_only <- nrow(discordant[discordant$TF_short > 0.03, ])
short_pos <- nrow(discordant[discordant$TF_short > 0.03 & discordant$cancer_status == "positive", ])
short_neg <- nrow(discordant[discordant$TF_short > 0.03 & discordant$cancer_status == "negative", ])

discordant <- nrow(discordant)

### Calculate stats of ichorCNA
neg <- ichorCNA[ichorCNA$cancer_status == "negative", ]
ichor_neg_neg <- nrow(neg[neg$TF < 0.03, ])
ichor_neg_pos <- nrow(neg[neg$TF > 0.03, ])
neg_samples <- nrow(neg)
neg_patients <- length(unique(neg$ext_ID))

pos <- ichorCNA[ichorCNA$cancer_status == "positive", ]
ichor_pos_neg <- nrow(pos[pos$TF < 0.03, ])
ichor_pos_pos <- nrow(pos[pos$TF > 0.03, ])
pos_samples <- nrow(pos)
pos_patients <- length(unique(pos$ext_ID))

### Compare late and early stage samples
samples_pos <- ichorCNA[ichorCNA$cancer_status == "positive", ]
samples_pos$ichor <- ifelse(samples_pos$TF > 0.03 | samples_pos$TF_short > 0.03, "pos", "neg")
pos_late <- nrow(samples_pos[samples_pos$stage == "high", ])
pos_early <- nrow(samples_pos[samples_pos$stage == "low", ])
pos_late_pos <- nrow(samples_pos[samples_pos$stage == "high" & samples_pos$ichor == "pos", ])
pos_early_pos <- nrow(samples_pos[samples_pos$stage == "low" & samples_pos$ichor == "pos", ])

### Compare ichorCNA to mutations
neg_mutations <- merge(neg, mutations, by.x = "TS", by.y = "sample_ID")
neg_mut_samples <- nrow(neg_mutations)
neg_mutations_neg_neg <- nrow(neg_mutations[neg_mutations$TF < 0.03 & neg_mutations$TF_short < 0.03 & neg_mutations$TP53_somatic == "", ])
neg_mutations_neg_pos <- nrow(neg_mutations[neg_mutations$TF < 0.03 & neg_mutations$TF_short < 0.03 & !(neg_mutations$TP53_somatic == ""), ])
neg_mutations_pos_neg <- nrow(neg_mutations[(neg_mutations$TF > 0.03 | neg_mutations$TF_short > 0.03) & neg_mutations$TP53_somatic == "", ])
neg_mutations_pos_pos <- nrow(neg_mutations[(neg_mutations$TF > 0.03 | neg_mutations$TF_short > 0.03) & !(neg_mutations$TP53_somatic == ""), ])

pos_mutations <- merge(pos, mutations, by.x = "TS", by.y = "sample_ID")
pos_mut_samples <- nrow(pos_mutations)
pos_mutations_neg_neg <- nrow(pos_mutations[pos_mutations$TF < 0.03 & pos_mutations$TF_short < 0.03 & pos_mutations$TP53_somatic == "", ])
pos_mutations_neg_pos <- nrow(pos_mutations[pos_mutations$TF < 0.03 & pos_mutations$TF_short < 0.03 & !(pos_mutations$TP53_somatic == ""), ])
pos_mutations_pos_neg <- nrow(pos_mutations[(pos_mutations$TF > 0.03 | pos_mutations$TF_short > 0.03) & pos_mutations$TP53_somatic == "", ])
pos_mutations_pos_pos <- nrow(pos_mutations[(pos_mutations$TF > 0.03 | pos_mutations$TF_short > 0.03) & !(pos_mutations$TP53_somatic == ""), ])

### Sum up statistics into summaries
# general stats
descriptors <- c("samples", "patients", "pediatric samples", "pediatric patients", "adult", "adult patients", "full ichorCNA positive", "cancer positive", "cancer negative")
stats <- c(cohort, patients, pediatric_samples, pediatric_patients, adult_samples, adult_patients, full_positive, full_positive_pos, full_positive_neg)
stats_cohort <- as.data.frame(cbind(descriptors, stats))

# comparisons
descriptors <- c("concordant", "cancer negative", "cancer positive", "discordant", "cancer negative", "cancer positive", 
                 "full only", "full only cancer negative", "full only cancer positive", "short only", "short only cancer negative", "short only cancer positive")
stats <- c(concordant, concordant_neg, concordant_pos, discordant, discordant_neg, discordant_pos, full_only, full_neg, full_pos, short_only, short_neg, short_pos)
stats_comparison <- as.data.frame(cbind(descriptors, stats))

#negative
ichor <- c("samples", "patients", "negative", "positive", "overlap", "negative", "negative", "positive", "positive")
muts <- c(NA, NA, NA, NA, NA, "negative", "positive", "negative", "positive")
stats <- c(neg_samples, neg_patients, ichor_neg_neg, ichor_neg_pos, neg_mut_samples, neg_mutations_neg_neg, neg_mutations_neg_pos, neg_mutations_pos_neg, neg_mutations_pos_pos)
stats_negative <- as.data.frame(cbind(ichor, muts, stats))
stats_negative$percent <- as.numeric(stats_negative$stats)/as.numeric(c(NA, NA, rep(neg_samples, 3), rep(neg_mut_samples, 4)))*100

#positive
ichor <- c("samples", "patients", "negative", "positive", "late stage samples", "late stage positive", "early stage samples", "early stage positive", 
           "overlap", "negative", "negative", "positive", "positive")
muts <- c(rep(NA, 9), "negative", "positive", "negative", "positive")
stats <- c(pos_samples, pos_patients, ichor_pos_neg, ichor_pos_pos, pos_late, pos_late_pos, pos_early, pos_early_pos,
           pos_mut_samples, pos_mutations_neg_neg, pos_mutations_neg_pos, pos_mutations_pos_neg, pos_mutations_pos_pos)
stats_positive <- as.data.frame(cbind(ichor, muts, stats))
stats_positive$percent <- as.numeric(stats_positive$stats)/as.numeric(c(NA, NA, rep(pos_samples, 2), NA, pos_late, NA, pos_early, pos_samples, rep(pos_mut_samples, 4)))*100

# total
ichor <- c("overlap", "positive samples", "negative samples", "negative", "negative", "positive", "positive")
muts <- c(NA, NA, NA, "negative", "positive", "negative", "positive")
stats <- c(neg_mut_samples + pos_mut_samples, 
           pos_mut_samples, neg_mut_samples,
           neg_mutations_neg_neg + pos_mutations_neg_neg, 
           neg_mutations_neg_pos + pos_mutations_neg_pos, 
           neg_mutations_pos_neg + pos_mutations_pos_neg, 
           neg_mutations_pos_pos + pos_mutations_pos_pos)
stats_total <- as.data.frame(cbind(ichor, muts, stats))
