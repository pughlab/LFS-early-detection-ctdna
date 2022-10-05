library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM_LFS_cfmedip_medremix/LFS_QC"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/insert_size"
project <- "CHARM_LFS"

filenames <- list.files(path = path, pattern = "picardInsertSize_metrics.txt", full.names = TRUE, recursive = TRUE)
names <- list.files(path = path, pattern = "picardInsertSize_metrics.txt", full.names = FALSE, recursive = TRUE)
names <- gsub(".*/", "", names)
names <- sub("(_CM).*", '\\1', names)

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Read in data
datalist <- lapply(filenames, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                     header = TRUE, col.names = c("length", "freq", "X", "Y"))})
data <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)

### Format fragment data and save
row.names(data) <- data$length
data <- data[,-1]
colnames(data) <- names
lengths <- c(30:600)
data <- data[row.names(data) %in% lengths,]
data[is.na(data)] <- 0
data$length <- row.names(data)
data <- data %>% dplyr::select(length, everything())
write.table(data, file.path(outdir, paste0(project, "_fragment_cfmedip.txt")), row.names = FALSE, sep = "\t")

### Calculate fragment frequencies
freq <- data[, -1]
sums <- colSums2(as.matrix(freq[row.names(freq) %in% c(30:220), ]))
freq <- as.data.frame(t(t(freq)/sums*100))
freq$length <- row.names(freq)
freq <- freq %>% dplyr::select(length, everything())
write.table(freq, file.path(outdir, paste0(project, "_fragment_freq_cfmedip.txt")), row.names = FALSE, sep = "\t")

### Check frequency distributions
data_melt <- melt(freq, value = "length")
data_melt$length <- as.numeric(data_melt$length)

freq_plot <- ggplot(data_melt, aes(x = length, y = value, group = variable)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  xlab("Fragment Size") + 
  ylab("Frequency (%)") +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

ggsave(file.path(outdir, "Frequency_check_cfmedip.pdf"), freq_plot, device = "pdf", width = 12, height = 100, units = "in", limitsize = FALSE)
