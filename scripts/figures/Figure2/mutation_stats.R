library(tidyverse)

### Set variables
path <- ""

mutations <- read.delim("Oncoplot_full.txt")
samples <- read.delim("sample_list.txt")

### Format data
mutations <- merge(mutations, samples, by.x = "sample_ID", by.y = "TS")

### Remove failed samples
exclude <- c("TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0316_Pl_T_PE_301_CM", "TGL49_0060_Cf_n_PE_301_CM")
fail <- c("TGL49_0010_Cf_U_PE_334_TS", "TGL49_0025_Cf_U_PE_334_TS") ### Same sample but reparented
mutations <- mutations[!(mutations$sample_ID %in% exclude), ]
sequenced <- nrow(mutations)
mutations <- mutations[!(mutations$sample_ID %in% fail), ]

samples <- samples[!(samples$cfMeDIP %in% exclude), ]

### Calculate mutation stats (general)
cohort <- nrow(mutations)
germline <- nrow(mutations[!(mutations$TP53_germline %in% c("", "deletion", "duplication")), ])
germline_cnv <- nrow(mutations[mutations$TP53_germline %in% c("deletion", "duplication"), ])
no_germline <- cohort - germline - germline_cnv

noseq <- samples[!(samples$notes == "NS"), ]
noseq <- noseq[noseq$TS == "", ]
noseq_adult <- nrow(noseq[noseq$Age == "adult", ])
noseq_ped <- nrow(noseq[noseq$Age == "pediatric", ])
noseq <- nrow(noseq)

patients <- length(unique(mutations$sample_parent))
serial_patients <- mutations[duplicated(mutations$sample_parent), ]
serial_patients <- length(unique(serial_patients$sample_parent))

adult <- mutations[mutations$Age == "adult", ]
adult_samples <- nrow(adult)
adult <- length(unique(adult$sample_parent))

pediatric <- mutations[mutations$Age == "pediatric", ]
pediatric_samples <- nrow(pediatric)
pediatric <- length(unique(pediatric$sample_parent))

### Calculate mutation stats (cancer status)
neg <- mutations[mutations$cancer_status == "negative", ]
neg_TP53 <- nrow(neg[str_detect(neg$TP53_somatic, ";"), ])
neg_TP53_samples <- nrow(neg[!(neg$TP53_somatic == ""), ])
neg_TP53_patients <- neg[!(neg$TP53_somatic == ""), ]
neg_TP53_patients <- length(unique(neg_TP53_patients$sample_parent))
neg_TP53_total <- neg_TP53 + neg_TP53_samples
neg_BRCA1 <- nrow(neg[!(neg$BRCA1 == ""), ])
neg_BRCA2 <- nrow(neg[!(neg$BRCA2 == ""), ])
neg_PALB2 <- nrow(neg[!(neg$PALB2 == ""), ])
neg_MLH1 <- nrow(neg[!(neg$MLH1 == ""), ])
neg_MSH2 <- nrow(neg[!(neg$MSH2 == ""), ])
neg_MSH6 <- nrow(neg[!(neg$MSH6 == ""), ])
neg_PMS2 <- nrow(neg[!(neg$PMS2 == ""), ])
neg_EPCAM <- nrow(neg[!(neg$EPCAM %in% c("", NA)), ])
neg_APC <- nrow(neg[!(neg$APC == ""), ])
neg_samples <- nrow(neg)
neg_patients <- length(unique(neg$sample_parent))

pos <- mutations[mutations$cancer_status == "positive", ]
pos_TP53 <- nrow(pos[str_detect(pos$TP53_somatic, ";"), ])
pos_TP53_samples <- nrow(pos[!(pos$TP53_somatic == ""), ])
pos_TP53_patients <- pos[!(pos$TP53_somatic == ""), ]
pos_TP53_patients <- length(unique(pos_TP53_patients$sample_parent))
pos_TP53_total <- pos_TP53 + pos_TP53_samples
pos_BRCA1 <- nrow(pos[!(pos$BRCA1 == ""), ])
pos_BRCA2 <- nrow(pos[!(pos$BRCA2 == ""), ])
pos_PALB2 <- nrow(pos[!(pos$PALB2 == ""), ])
pos_MLH1 <- nrow(pos[!(pos$MLH1 == ""), ])
pos_MSH2 <- nrow(pos[!(pos$MSH2 == ""), ])
pos_MSH6 <- nrow(pos[!(pos$MSH6 == ""), ])
pos_PMS2 <- nrow(pos[!(pos$PMS2 == ""), ])
pos_EPCAM <- nrow(pos[!(pos$EPCAM %in% c("", NA)), ])
pos_APC <- nrow(pos[!(pos$APC == ""), ])
pos_samples <- nrow(pos)
pos_patients <- length(unique(pos$sample_parent))

high <- nrow(pos[pos$stage == "high", ])
low <- nrow(pos[pos$stage == "low", ])
none <- pos_samples - high - low

pos_pos <- pos[!(pos$TP53_somatic == ""), ]
pos_high <- nrow(pos_pos[pos_pos$stage == "high", ])
pos_low <- nrow(pos_pos[pos_pos$stage == "low", ])

### Sum up statistics into summaries
#general 
descriptors <- c("samples sequenced", "samples passed", "total patients", "postive samples", "negative samples", "patients with serial samples", "not sequenced", "noseq adult",
                 "noseq ped",  "germline mutation", "germline CNV", "no germline", "pediatric samples", "pediatric patients", "adult samples", "adult patients")
stats <- c(sequenced, cohort, patients, pos_samples, neg_samples, serial_patients, noseq, noseq_adult, noseq_ped, germline, germline_cnv, 
           no_germline, pediatric_samples, pediatric, adult_samples, adult)
stats_cohort <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, NA, NA, cohort, cohort, patients, cohort, noseq, noseq, cohort, cohort, cohort, cohort, patients, cohort, patients)
stats_cohort$percent <- as.numeric(stats_cohort$stats)/denominators*100

#negative
descriptors <- c("negative", "patients", "TP53 mutations", "samples with TP53 mutations", "patients with TP53 mutations", "patients with multiple TP53", 
                 "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC")
stats <- c(neg_samples, neg_patients, neg_TP53_total, neg_TP53_samples, neg_TP53_patients, neg_TP53, 
           neg_BRCA1, neg_BRCA2, neg_PALB2, neg_MLH1, neg_MSH2, neg_MSH6, neg_PMS2, neg_EPCAM, neg_APC)
stats_negative <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, NA, NA, neg_samples, neg_patients, neg_patients, rep(neg_samples, 9))
stats_negative$percent <- as.numeric(stats_negative$stats)/denominators*100

#positive
descriptors <- c("positive", "patients", "late stage samples", "early stage samples", "no stage samples", "TP53 mutations", "samples with TP53 mutations", "patients with TP53 mutations", 
                 "late stage mutation positive", "early stage mutation positive", "patients with multiple TP53", "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC")
stats <- c(pos_samples, pos_patients, high, low, none, pos_TP53_total, pos_TP53_samples, pos_TP53_patients, pos_high, pos_low, pos_TP53, 
           pos_BRCA1, pos_BRCA2, pos_PALB2, pos_MLH1, pos_MSH2, pos_MSH6, pos_PMS2, pos_EPCAM, pos_APC)
stats_positive <- as.data.frame(cbind(descriptors, stats))
denominators <- c(NA, NA, rep(pos_samples, 5), pos_patients, high, low, pos_patients, rep(pos_samples, 9))
stats_positive$percent <- as.numeric(stats_positive$stats)/denominators*100
