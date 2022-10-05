library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/ichorCNA/output"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_project/LFS/samples/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_project/LFS/ichorCNA"
project <- "CHARM_LFS"

ichorCNA <- file.path(path, "ichorCNA")
ichor_short <- file.path(path, "ichorCNA_short")

files <- list.files(path = ichorCNA, pattern = "*params.txt", full.names = TRUE)
files_short <- list.files(path = ichor_short, pattern = "*params.txt", full.names = TRUE)

names <- list.files(path = ichorCNA, pattern = "*params.txt", full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

names_short <- list.files(path = ichor_short, pattern = "*params.txt", full.names = FALSE)
names_short <- sub("(_WG).*", '\\1', names_short)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)

datalist <- lapply(files, function(x){read.delim(file = x, header = TRUE)[,2]})
data <- do.call(cbind, datalist)

datalist <- lapply(files_short, function(x){read.delim(file = x, header = TRUE)[,2]})
data_short <- do.call(cbind, datalist)

### Format ichorCNA data and write data
data <- as.data.frame(data)
data <- data[c(3:8,10,11,14,15), ]
row.names(data) <- c("Sex", "TF", "Ploidy", "Subclone_Fraction", "Fraction_genome_subclonal", "fraction_cna_subcloncal", 
                     "chrY_fraction", "chrX_median_log_ratio", "gamma_rate_init", "GC_map_correction_mad")
data <- as.data.frame(t(data))
data$sample <- names
data[is.na(data)] <- 0
data[data == "NaN"] <- 0
data <- data %>% dplyr::select(sample, everything())
data <- data[data$sample %in% samples$sWGS, ]

write.table(data, file.path(outdir, paste0(project, "_ichorCNA.txt")), row.names = FALSE, sep = "\t")

data_short <- as.data.frame(data_short)
data_short <- data_short[c(3:8,10,11,14,15), ]
row.names(data_short) <- c("Sex", "TF", "Ploidy", "Subclone_Fraction", "Fraction_genome_subclonal", "fraction_cna_subcloncal", 
                           "chrY_fraction", "chrX_median_log_ratio", "gamma_rate_init", "GC_map_correction_mad")
data_short <- as.data.frame(t(data_short))
data_short$sample <- names_short
data_short[is.na(data_short)] <- 0
data_short[data_short == "NaN"] <- 0
data_short <- data_short %>% dplyr::select(sample, everything())
data_short <- data_short[data_short$sample %in% samples$sWGS, ]

write.table(data_short, file.path(outdir, paste0(project, "_ichorCNA_short.txt")), row.names = FALSE, sep = "\t")

### Merge and summarize ichorCNA runs
data <- data[ , c(1:4)]
data_short <- data_short[ , c(1,3,4)]
data <- merge(data, data_short, by = c("sample"), all = TRUE)
colnames(data) <- c("sample", "sex", "TF", "ploidy", "TF_short", "ploidy_short")
data <- merge(data, samples, by.x = "sample", by.y = "sWGS")
write.table(data, file.path(outdir, paste0(project, "_ichorCNA_summary.txt")), row.names = FALSE, sep = "\t")

### Import IchorCNA copy number depths
files <- list.files(path = ichorCNA, pattern = "*Depth.txt", full.names = TRUE)
files_short <- list.files(path = ichor_short, pattern = "*Depth.txt", full.names = TRUE)

rows <- files[1]
rows <- read.delim(rows, header = TRUE)
rows <- rows[, 1:3]

datalist <- lapply(files, function(x){read.delim(file = x, header = TRUE)[,4]})
data <- do.call(cbind, datalist)

datalist <- lapply(files_short, function(x){read.delim(file = x, header = TRUE)[,4]})
data_short <- do.call(cbind, datalist)

### Format IchorCNA copy number depths and write data
data <- as.data.frame(data)
colnames(data) <- names
data <- data[ , colnames(data) %in% samples$sWGS]
data <- cbind(rows, data)

write.table(data, file.path(outdir, paste0(project, "_corrected_depths.txt")), row.names = FALSE, sep = "\t")

data_short <- as.data.frame(data_short)
colnames(data_short) <- names_short
data_short <- data_short[ , colnames(data_short) %in% samples$sWGS]
data_short <- cbind(rows, data_short)

write.table(data_short, file.path(outdir, paste0(project, "_corrected_depths_short.txt")), row.names = FALSE, sep = "\t")
