library(tidyverse)
library(ggplot2)

# Set variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/mutation_fragment/output"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/mutation_fragment"
options(scipen=999)

sample <- "TGL49_0025_Pb_n_PE_452_TS"
gene <- "TP53"
mutation <- "X187_splice"

file <- list.files(path, pattern = sample, full.names = TRUE)
name <- list.files(path, pattern = sample, full.names = FALSE)
wildtype <- file.path(file, paste0(name, "_", gene, "_", mutation, "_wildtype"))
mutant <- file.path(file, paste0(name, "_", gene, "_", mutation, "_mutant"))

# Read in files
wildtype_reads <- read.delim(wildtype, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
wildtype_reads <- abs(wildtype_reads[complete.cases(wildtype_reads), ])
wildtype_reads <- as.data.frame(wildtype_reads)
mutant_reads <- read.delim(mutant, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
mutant_reads <- abs(mutant_reads[complete.cases(mutant_reads), ])
mutant_reads <- as.data.frame(mutant_reads)
colnames(mutant_reads) <- "mutant_reads"

min <- ifelse(min(mutant_reads) < min(wildtype_reads), min(mutant_reads), min(wildtype_reads))
max <- ifelse(max(mutant_reads) > max(wildtype_reads), max(mutant_reads), max(wildtype_reads))
p <- ks.test(wildtype_reads$wildtype_reads, mutant_reads$mutant_reads)$p.value

wt <- nrow(wildtype_reads)
mut <- nrow(mutant_reads)
median1 <- median(wildtype_reads$wildtype_reads)
median2 <- median(mutant_reads$mutant_reads)
sd1 <- sd(wildtype_reads$wildtype_reads)
sd2 <- sd(mutant_reads$mutant_reads)
t_test <- (median1 - median2)/sqrt(((sd1^2)/wt) + ((sd2^2)/mut))
t_test <- pnorm(t_test, mean = 0, sd = 1, lower.tail = FALSE)

p <- ifelse(p < t_test, p, t_test)
p <- round(p, 4)

# Plot cumulative distributions
plot <- ggplot() + 
  stat_ecdf(data = wildtype_reads, aes(wildtype_reads, color = "black"), geom = "step") +
  stat_ecdf(data = mutant_reads, aes(mutant_reads, color = "red"), geom = "step") +
  xlab("Fragment Length") + 
  ylab("Cumulative Distribution") +
  ggtitle(mutation) + 
  scale_colour_manual(name = 'Fragments', values =c("black" = "black", "red" = "red"), labels = c("Wildtype", mutation)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_x_continuous(limits=c(75, 250, expand = c(0,0))) +
  geom_vline(xintercept = median1, linetype = "dashed", color = "black", size=0.5) +
  geom_vline(xintercept = median2, linetype = "dashed", color = "red", size=0.5) +
  annotate("text", x = 215, y = c(0.40, 0.27, 0.14), 
           label = c(paste0("p = ", p), paste0("WT = ", wt), paste0(mutation, " = ", mut)),
           col = c("black", "black", "red"),
           size = 4)
plot

ggsave(file.path(outdir, paste0(sample, "_", gene, "_", mutation, ".pdf")), plot, device = "pdf", width = 3.25, height = 2.5, units = "in")

