library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/gatk_coverage/output"
samples <- "//Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/GATK_coverage"
project <- "CHARM_LFS"

### Get file lists and names
coverage_swgs <- list.files(file.path(path, "swgs"), pattern = "sample_summary", full.names = TRUE)
coverage_short <- list.files(file.path(path, "short"), pattern = "sample_summary", full.names = TRUE)
coverage_all_unique <- list.files(file.path(path, "all_unique"), pattern = "sample_summary", full.names = TRUE)
coverage_ts_raw <- list.files(file.path(path, "ts_raw"), pattern = "sample_summary", full.names = TRUE)

names_swgs <- list.files(file.path(path, "swgs"), pattern = "sample_summary", full.names = FALSE)
names_swgs <- gsub("(_WG).*", "\\1", names_swgs)
names_short <- list.files(file.path(path, "short"), pattern = "sample_summary", full.names = FALSE)
names_short <- gsub("(_WG).*", "\\1", names_short)
names_all_unique <- list.files(file.path(path, "all_unique"), pattern = "sample_summary", full.names = FALSE)
names_all_unique <- gsub("(_TS).*", "\\1", names_all_unique)
names_ts_raw <- list.files(file.path(path, "ts_raw"), pattern = "sample_summary", full.names = FALSE)
names_ts_raw <- gsub("(_TS).*", "\\1", names_ts_raw)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in files
samples <- read.delim(samples, check.names = FALSE)

datalist <- lapply(coverage_swgs, function(x){read.delim(file = x, nrows = 1, header = TRUE)})
coverage_swgs <- ldply(datalist)

datalist <- lapply(coverage_short, function(x){read.delim(file = x, nrows = 1, header = TRUE)})
coverage_short <- ldply(datalist)

datalist <- lapply(coverage_all_unique, function(x){read.delim(file = x, nrows = 1, header = TRUE, col.names = c("sample"))})
coverage_all_unique <- ldply(datalist)
coverage_all_unique <- data.frame(do.call('rbind', strsplit(as.character(coverage_all_unique$sample),',',fixed=TRUE)))

datalist <- lapply(coverage_ts_raw, function(x){read.delim(file = x, nrows = 1, header = TRUE, col.names = c("sample"))})
coverage_ts_raw <- ldply(datalist)
coverage_ts_raw <- data.frame(do.call('rbind', strsplit(as.character(coverage_ts_raw$sample),',',fixed=TRUE)))

### Format tables and write data
coverage_swgs$sample_id <- names_swgs
coverage_swgs <- coverage_swgs[ , c("sample_id", "mean")]

coverage_short$sample_id <- names_short
coverage_short <- coverage_short[ , c("sample_id", "mean")]

coverage_all_unique$sample_id <- names_all_unique
coverage_all_unique <- coverage_all_unique[ , c("sample_id", "X3")]

coverage_ts_raw$sample_id <- names_ts_raw
coverage_ts_raw <- coverage_ts_raw[ , c("sample_id", "X3")]

### Merge sequencing coverages
summary <- merge(samples, coverage_swgs, by.x = "sWGS", by.y = "sample_id", all = TRUE)
summary <- merge(summary, coverage_short, by.x = "sWGS", by.y = "sample_id", all = TRUE)
summary <- merge(summary, coverage_all_unique, by.x = "TS", by.y = "sample_id", all = TRUE)
summary <- merge(summary, coverage_ts_raw, by.x = "TS", by.y = "sample_id", all = TRUE)
summary <- summary[ , c("Patient_ID", "LIB_ID", "sample_parent", "timepoint", "family", "Age", "TS", "sWGS", "mean.x", "mean.y", "X3.x", "X3.y")]
colnames(summary) <- c("Patient_ID", "LIB_ID", "sample_parent", "timepoint", "family", "Age", "TS", "sWGS", "coverage_swgs", "coverage_short", "coverage_all_unique", "coverage_ts")

summary <- summary[order(summary$sWGS), ]
write.table(summary, file.path(outdir, paste0(project, "_coverage_summary.txt")), row.names = FALSE, sep = "\t")
