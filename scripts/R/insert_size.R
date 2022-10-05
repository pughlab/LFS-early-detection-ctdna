library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/insert_size/output/swgs"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/insert_size"
mouliere <- "/Users/derekwong/OneDrive - UHN/Post-Doc/External_data/Mouliere_fragment/Mouliere_fragment.txt"
project <- "CHARM_LFS"

filenames <- list.files(path = path, pattern = "*picard.txt", full.names = TRUE)
names <- list.files(path = path, pattern = "*.txt", full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)
mouliere <- read.delim(mouliere)

datalist <- lapply(filenames, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                     header = TRUE, col.names = c("length", "freq", "X", "Y"))})
data <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)

### Format fragment data and save
row.names(data) <- data$length
data <- data[,-1]
colnames(data) <- names
lengths <- c(10:600)
data <- data[row.names(data) %in% lengths,]
data[is.na(data)] <- 0
data$length <- row.names(data)
data <- data %>% dplyr::select(length, everything())
write.table(data, file.path(outdir, paste0(project, "_fragment.txt")), row.names = FALSE, sep = "\t")

### Calculate fragment frequencies
freq <- data[, -1]
sums <- colSums2(as.matrix(freq[row.names(freq) %in% c(10:250), ]))
freq <- as.data.frame(t(t(freq)/sums*100))
freq$length <- row.names(freq)
freq <- freq %>% dplyr::select(length, everything())
write.table(freq, file.path(outdir, paste0(project, "_fragment_freq.txt")), row.names = FALSE, sep = "\t")

### Calculate fragment proportions
short <- freq[ , -1]
short <- short[rownames(short) %in% c(20:150), ]
prop <- colSums2(as.matrix(short))/100
prop <- as.data.frame(prop)

### Format fragment proportions and save
prop$sample <- colnames(short)
prop <- merge(prop, samples, by.x = "sample", by.y = "sWGS")
prop$cancer <- prop$cancer_status
prop$type <- "LFS"
prop <- prop[ , c("sample", "cancer", "type", "prop")]
colnames(prop) <- c("sample", "cancer", "type", "P.20_150.")
prop <- rbind(prop, mouliere)
write.table(prop, file.path(outdir, paste0(project, "_proportions.txt")), row.names = FALSE, sep = "\t")

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

ggsave(file.path(outdir, "Frequency_check.pdf"), freq_plot, device = "pdf", width = 12, height = 100, units = "in", limitsize = FALSE)
