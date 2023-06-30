library(tidyverse)

### Set paths
path <- ""

### Find files
breast_scores <- read_tsv(file.path(path, "methylation_score_breast.tsv"))
all_scores <- read_tsv(file.path(path, "pancancer_score_vrba.tsv"))
data_samples <- read.delim("sample_list.txt")

### Calculate thresholds
limit_breast <- quantile(breast_scores$score[breast_scores$cancer_type == "HBC"], 0.95)
limit_all <- quantile(all_scores$score[all_scores$cancer_type == "HBC"], 0.95)

### Format data
exclude <- c("TGL49_0010_Cf_U_PE_319_CM", "TGL49_0209_Cf_U_PE_355_CM")
data_samples <- data_samples[!(data_samples$cfMeDIP %in% exclude), ]
data_samples <- data_samples[!(data_samples$cfMeDIP == ""), ]
breast_scores <- breast_scores[!(breast_scores$cfMeDIP %in% exclude | breast_scores$cfMeDIP == "HBC"), ]
breast_scores <- merge(breast_scores, data_samples, by = "cfMeDIP")
all_scores <- all_scores[!(all_scores$cfMeDIP %in% exclude | all_scores$cfMeDIP == "HBC"), ]
all_scores <- merge(all_scores, data_samples, by = "cfMeDIP")

### Calculate general statistics
cohort <- nrow(data_samples)
good_samples <- nrow(all_scores)
patients <- length(unique(data_samples$sample_parent))

adult <- data_samples[data_samples$Age == "adult", ]
adult_samples <- nrow(adult)
adult_patients <- length(unique(adult$sample_parent))

pediatric <- data_samples[data_samples$Age == "pediatric", ]
pediatric_samples <- nrow(pediatric)
pediatric_patients <- length(unique(pediatric$sample_parent))

cancer_positive <- nrow(data_samples[data_samples$cancer_status == "positive", ])
positive_patients <- length(unique(data_samples$sample_parent[data_samples$cancer_status == "positive"]))
cancer_negative <- nrow(data_samples[data_samples$cancer_status == "negative", ])
negative_patients <- length(unique(data_samples$sample_parent[data_samples$cancer_status == "negative"]))

### Sum up statistics into summaries (general)
desc <- c("samples", "patients", "adult samples", "adult patients", "pediatric samples", "pediatric patients", "cancer positive samples", "cancer positive patients",
          "cancer negative samples", "cancer negative patients")
stats <- c(cohort, patients, adult_samples, adult_patients, pediatric_samples, pediatric_patients, cancer_positive, positive_patients, cancer_negative, negative_patients)
stats_general <- as.data.frame(cbind(desc, stats))

### Calculate breast statistics
positive <- breast_scores[breast_scores$score > limit_breast, ]
positive_pos <- nrow(positive[positive$cancer_status == "positive", ])
positive_pos_breast <- nrow(positive[positive$cancer_type.y == "breast", ])
positive_pos_breast_high <- nrow(positive[positive$cancer_type.y == "breast" & positive$stage == "high", ])
positive_pos_breast_low <- nrow(positive[positive$cancer_type.y == "breast" & positive$stage == "low", ])
positive_pos_other <- positive$cancer_type.y[positive$cancer_status == "positive"]
positive_pos_other <- positive_pos_other[!(grepl("breast", positive_pos_other))]
positive_pos_other <- do.call(paste, c(as.list(positive_pos_other), sep = ", "))

breast <- breast_scores[breast_scores$cancer_type.y == "breast", ]
breast_high <- nrow(breast[breast$stage == "high", ])
breast_low <- nrow(breast[breast$stage == "low", ])
breast_un <- nrow(breast[breast$stage == "unknown", ])
breast <- nrow(breast)

positive_neg <- nrow(positive[positive$cancer_status == "negative", ])
positive_neg_prev <- nrow(positive[positive$cancer_status == "negative" & positive$previous_cancer == "yes", ])
positive_neg_breast <- nrow(positive[positive$cancer_status == "negative" & positive$previous_type %like% "breast", ])
positive_neg_other <- positive$previous_type[positive$cancer_status == "negative"]
positive_neg_other <- positive_neg_other[!(positive_neg_other %like% "breast" | positive_neg_other == "")]
positive_neg_other <- do.call(paste, c(as.list(positive_neg_other), sep = ", "))
positive <- nrow(positive)

breast_only <- breast_scores[breast_scores$cancer_type.y == "breast",]
breast_only <- merge(breast_only[,1:2], all_scores, by = "cfMeDIP")
breast_classifier_only <- nrow(breast_only[breast_only$score.x > limit_breast & breast_only$score.y < limit_all, ])
pancan_classifier_only <- nrow(breast_only[breast_only$score.x < limit_breast & breast_only$score.y > limit_all, ])

### Sum up statistics into summaries (breast)
desc <- c("samples", "passed samples", "patients", "breast samples", "high grade", "low grade", "unknown grade", "score negative", "score positive", "postive breast", "positive high grade",
          "positive low grade", "other cancer types", "score positive negative", "previous cancer", "previous breast", "previous other", "breast classifier only", "pancan classifier only")
stats <- c(cohort, good_samples, patients, breast, breast_high, breast_low, breast_un, good_samples - positive, positive, positive_pos_breast, positive_pos_breast_high, positive_pos_breast_low,
           positive_pos_other, positive_neg, positive_neg_prev, positive_neg_breast, positive_neg_other, breast_classifier_only, pancan_classifier_only)
stats_breast <- as.data.frame(cbind(desc, stats))

### Calculate pancancer statistics
positive <- all_scores[all_scores$score > limit_all, ]
positive_pos <- nrow(positive[positive$cancer_status == "positive", ])
positive_pos_all <- positive$cancer_type.y[positive$cancer_status == "positive" & positive$score > limit_all]
positive_pos_all <- do.call(paste, c(as.list(positive_pos_all), sep = ", "))
positive_pos_all_high <- nrow(positive[positive$stage == "high", ])
positive_pos_all_low <- nrow(positive[positive$stage == "low", ])

all_high <- nrow(all_scores[all_scores$stage == "high", ])
all_low <- nrow(all_scores[all_scores$stage == "low", ])
all_un <- nrow(all_scores[all_scores$stage == "unknown", ])

positive_neg <- nrow(positive[positive$cancer_status == "negative", ])
positive_neg_prev <- nrow(positive[positive$cancer_status == "negative" & positive$previous_cancer == "yes", ])
positive_neg_all <- positive$previous_type[positive$cancer_status == "negative"]
positive_neg_all <- do.call(paste, c(as.list(positive_neg_other), sep = ", "))
positive <- nrow(positive)

### Sum up statistics into summaries (pan cancer)
desc <- c("samples", "passed samples", "patients", "high grade", "low grade", "unknown grade", "score negative", "score positive", "cancer positive", "postive cancer types", "positive high grade",
          "positive low grade", "score positive negative", "previous cancer", "previous all")
stats <- c(cohort, good_samples, patients, all_high, all_low, all_un, good_samples - positive, positive, positive_pos, positive_pos_all, positive_pos_all_high, positive_pos_all_low,
           positive_neg, positive_neg_prev, positive_neg_all)
stats_all <- as.data.frame(cbind(desc, stats))

