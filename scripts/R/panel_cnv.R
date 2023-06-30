library(tidyverse)
library(dplyr)
library(readxl)

### Set working variables
path <- ""
outdir <- ""
project <- "CHARM_LFS"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Get Files
panelcnmops <- read.delim(list.files(path = file.path(path, "/pipeline_output/panelCNmops"), pattern = "mops_ratio_matrix.tsv", full.names = TRUE))
panel_buffy <- read.delim(list.files(path = file.path(path, "/pipeline_output/buffy_coat/panelCNmops"), pattern = "mops_ratio_matrix.tsv", full.names = TRUE))
viscap <- read.delim(list.files(path = file.path(path, "VisCap/VisCap_run1"), pattern = "log2_ratio_table.xls", full.names = TRUE))

### Format Viscap
colnames(viscap) <- gsub(".all.unique.dcs.sorted.bam.", "", colnames(viscap))
viscap$Chromosome <- sapply(strsplit(viscap$X, ":"), `[`, 1)
viscap$gene_exon <- sapply(strsplit(viscap$X, ":"), `[`, 2)
viscap$Start <- sapply(strsplit(viscap$gene_exon, "-"), `[`, 1)
viscap$End <- sapply(strsplit(viscap$gene_exon, "-"), `[`, 2)

viscap <- viscap[, !(colnames(viscap) %in% c("X", "gene_exon"))]
viscap <- viscap %>%
  select(Chromosome, Start, End, everything())
viscap$Start <- as.numeric(viscap$Start)
viscap$End <- as.numeric(viscap$End)
viscap <- viscap[order(factor(viscap$Chromosome, levels = c(paste0("chr", c(1:22)), "chrX")),
                       viscap$Start), ]

buffy_coats <- colnames(viscap)[colnames(viscap) %like% "Ly"]
viscap_bc <- viscap[, colnames(viscap) %in% c("Chromosome", "Start", "End", buffy_coats)]
viscap <- viscap[, !(colnames(viscap) %in% buffy_coats)]

### Format panelcnmops
colnames(panel_buffy) <- gsub("_tm", "", colnames(panel_buffy))

### Write tables
write.table(panel_buffy, file.path(outdir, paste0(project, "_panelCNmops_buffycoat.txt")), sep = "\t", row.names = FALSE)
write.table(viscap_bc, file.path(outdir, paste0(project, "_VisCap_buffycoat.txt")), sep = "\t", row.names = FALSE)

write.table(panelcnmops, file.path(outdir, paste0(project, "_panelCNmops.txt")), sep = "\t", row.names = FALSE)
write.table(viscap, file.path(outdir, paste0(project, "_VisCap.txt")), sep = "\t", row.names = FALSE)
