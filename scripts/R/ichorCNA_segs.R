library(tidyverse)
library(plyr)
library(data.table)

### Set working variables
path <- ""
samples <- "sample_list.txt"
outdir <- ""
project <- "CHARM_LFS"
ichorCNA <- file.path(path, "ichorCNA")

### Make outdir
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "TEMP"), showWarnings = FALSE)

### Get file names
names <- list.files(path = ichorCNA, pattern = ".correctedDepth.txt", full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

for (name in names) {
  ### Get files
  rows <- read.delim(file.path(ichorCNA, paste0(name, ".correctedDepth.txt")))
  rows <- rows[, colnames(rows) %in% c("chr", "start", "end")]
  rows$chr <- factor(rows$chr, levels = c(1:22, "X"))
  
  segs <- read.delim(file.path(ichorCNA, paste0(name, ".seg")))
  segs_loc <- segs[, colnames(segs) %in% c("chr", "start", "end")]
  segs_loc$chr <- factor(segs_loc$chr, levels = c(1:22, "X"))

  ### Find overlaps
  setDT(segs_loc)
  setDT(rows)
  setkey(segs_loc)
  overlap <- foverlaps(rows, segs_loc, type="within", nomatch=0L)

  ### Merge in seg information
  seg_table <- merge(overlap, segs[, colnames(segs) %in% c("chr", "start", "end", "median")], by = c("chr", "start", "end"))
  seg_table <- as.data.frame(seg_table)
  seg_table <- seg_table[, colnames(seg_table) %in% c("chr", "i.start", "i.end", "median")]
  colnames(seg_table) <- c("chr", "start", "end", "median")
  write.table(seg_table, file.path(outdir, "TEMP", paste0(name, "_seg.txt")), sep = "\t", row.names = FALSE)
}

### Read in sample information and seg files
samples <- read.delim(samples)
files <- list.files(file.path(outdir, "TEMP"), pattern = "_seg.txt", full.names = TRUE)

rows <- names[1]
rows <- read.delim(file.path(ichorCNA, paste0(rows, ".correctedDepth.txt")), header = TRUE)
data_seg <- rows[, 1:3]

for (file in files) {
  name <- gsub(pattern = paste0(outdir, "/TEMP/"), "", x = file)
  name <- sub("(_WG).*", '\\1', name)
  file <- read.delim(file)
  colnames(file) <- c("chr", "start", "end", name)
  data_seg <- merge(data_seg, file, by = c("chr", "start", "end"), all = TRUE)
}

data_seg <- data_seg[order(factor(data_seg$chr, levels = c(1:22, "X")),
                           data_seg$start), ]

write.table(data_seg, file.path(outdir, paste0(project, "_segs.txt")), row.names = FALSE, sep = "\t")
