library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA"
outdir <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/ichorCNA"

ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", full.names = TRUE))

### Format data
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
ichorCNA <- ichorCNA[!(ichorCNA$sample %in% exclude), ]
ichorCNA$ichor <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA$cancer_status <- factor(ichorCNA$cancer_status, levels = c("negative", "positive"),
                                 labels = c("Cancer Free", "Active Cancer"))

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               strip.background = element_rect(fill = NA),
               strip.text = element_text(size = 10),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.background = element_blank(),
               legend.position = c(0.75,0.25),
               legend.direction="vertical",
               axis.text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.title = element_text(size = 10))

### Plot data
comp_plot <- ggplot(ichorCNA, aes(TF, TF_short)) + 
  geom_point(aes(fill = cancer_status), stroke = 0, pch = 21, alpha = 0.5, size = 2) +
  ggtitle("") +
  xlab("All Fragments (TF)") +
  ylab("Short Fragments (TF)") +
  labs(fill = "Cancer Status") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 0.4, aes(label = ..rr.label..)) +
  theme + 
  scale_fill_manual(values = c("grey", "#fa4e0c")) +
  scale_y_continuous(limits = c(-0.05, 0.41)) + 
  scale_x_continuous(limits = c(-0.05, 0.41)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
comp_plot 

ggsave(file.path(outdir, paste0("ichorCNA_TFs.pdf")), comp_plot, device = "pdf", width = 3.25, height = 3.25, units = "in")

