library(dplyr)
library(matrixStats)
library(tidyverse)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/fragment_proportion"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

## Import data
data_prop <- read.delim(list.files(path, "prop", full.names = TRUE))
data_samples <- read.delim(samples)

### Restrict size (10-320bp)
data_prop <- data_prop[data_prop$type == "LFS", ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data_prop$sample, ]
data_prop <- data_prop[data_prop$sample %in% data_samples$sWGS, ]

### Keep only longitudinal patients
data_samples <- data_samples[data_samples$cancer_status == "negative" &
                               data_samples$previous_cancer == "no", ]
data_samples <- data_samples[data_samples$sample_parent %in% data_samples$sample_parent[duplicated(data_samples$sample_parent)], ]
positives <- data_samples$sample_parent[data_samples$cancer_status == "positive"]
seros <- c("LIB-04-0022", "LIB-04-0047", "LIB-04-0065")
data_samples <- data_samples[!(data_samples$sample_parent %in% c(positives, seros)), ]
data_samples <- data_samples %>%
  group_by(sample_parent) %>%
  arrange(timepoint)

### Merge clinical and fragment proportions
data_prop <- data_prop[data_prop$sample %in% data_samples$sWGS, ]
data <- merge(data_prop, data_samples, by.x = "sample", by.y = "sWGS")

### Set theme and colors
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13),
               axis.line = element_line(colour = "black"),
               axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.background = element_blank(),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

### Plot fragment props
plot_patients <- ggplot(data, aes(sample_parent, `P.20_150.`, group = sample_parent)) +
  geom_line(position = position_dodge2(0.5)) +
  geom_point(aes(color = cancer_status), 
             alpha = 0.75, position = position_dodge2(0.5), size = 1.5, pch = 16) +
  scale_color_manual(labels = c("Negative", "Positive"), values = c("black", "red")) +
  xlab("Patient") + 
  ylab("Proportion of fragments >150bp") +
  ggtitle("LFS Previvors") + 
  theme +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0,0))
plot_patients

ggsave(file.path(outdir, "fragment_proportions_patient.pdf"), plot_patients, device = "pdf", width = 6, height = 4, units = "in")


