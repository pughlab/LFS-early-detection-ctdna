library(tidyverse)

### Set variables
path <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_frequency"

data <- read.delim(list.files(path, "fragment_scores.txt", full.names = TRUE))

### Find positive predictive values
positive <- nrow(data[data$cancer_status == "Positive", ])
negative <- nrow(data[data$cancer_status == "Negative", ])

pos <- nrow(data[data$cancer_status == "Positive" &
                   data$zscores > data$limit, ])
pos0 <- nrow(data[data$cancer_status == "Positive" &
                   data$zscores0 > data$limit0, ])
pos1 <- nrow(data[data$cancer_status == "Positive" &
                   data$zscores1 > data$limit1, ])
pos2 <- nrow(data[data$cancer_status == "Positive" &
                   data$zscores2 > data$limit2, ])
pos3 <- nrow(data[data$cancer_status == "Positive" &
                   data$zscores3 > data$limit3, ])
pos_any <- nrow(data[data$cancer_status == "Positive" &
                   (data$zscores0 > data$limit0 |
                      data$zscores1 > data$limit1 |
                      data$zscores2 > data$limit2 |
                      data$zscores3 > data$limit3 |
                      data$zscores > data$limit), ])

neg <- nrow(data[data$cancer_status == "Negative" &
                   data$zscores > data$limit, ])
neg0 <- nrow(data[data$cancer_status == "Negative" &
                    data$zscores0 > data$limit0, ])
neg1 <- nrow(data[data$cancer_status == "Negative" &
                    data$zscores1 > data$limit1, ])
neg2 <- nrow(data[data$cancer_status == "Negative" &
                    data$zscores2 > data$limit2, ])
neg3 <- nrow(data[data$cancer_status == "Negative" &
                    data$zscores3 > data$limit3, ])
neg_any <- nrow(data[data$cancer_status == "Negative" &
                       (data$zscores0 > data$limit0 |
                          data$zscores1 > data$limit1 |
                          data$zscores2 > data$limit2 |
                          data$zscores3 > data$limit3 |
                          data$zscores > data$limit), ])

### Summarize stats
descriptors <- c("positive samples", "integrated", "comp 1", "comp 2", "comp 3", "comp4", "any",
                 "negative samples", "integrated", "comp 1", "comp 2", "comp 3", "comp4", "any")
stats <- c(positive, pos, pos0, pos1, pos2, pos3, pos_any, negative, neg, neg0, neg1, neg2, neg3, neg_any)
totals <- c(NA, rep(positive, 6), NA, rep(negative, 6))
stats_scores <- as.data.frame(cbind(descriptors, stats))
stats_scores$percent <- as.numeric(stats_scores$stats)/totals
