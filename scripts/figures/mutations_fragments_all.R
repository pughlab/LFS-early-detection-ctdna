library(tidyverse)
library(ggplot2)

# Set variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/mutation_fragment/output"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/mutation_fragment"
options(scipen=999)

### Import mutation list and format into germline/somatic
mutation_list <- read.delim(file.path(outdir, "mutation_matrix.txt"))

germline_mutations <- mutation_list[ , c("sample_ID", "germline")]
germline_mutations <- germline_mutations[!(germline_mutations$germline == ""), ]

somatic_mutations <- mutation_list[ , c("sample_ID", "somatic")]
somatic_mutations <- somatic_mutations[!(somatic_mutations$somatic == ""), ]

### Loop over germline mutations and get fragment lengths
wildtype_reads <- vector()
mutant_reads <- vector()

for(i in 1:nrow(germline_mutations)) {
  data <- germline_mutations[i, ]
  sample <- data[, 1]
  gene <- "TP53"
  mutation <- data[, 2]
  
  file <- list.files(path, pattern = sample, full.names = TRUE)
  name <- list.files(path, pattern = sample, full.names = FALSE)
  wildtype <- file.path(file, paste0(name, "_", gene, "_", mutation, "_wildtype"))
  mutant <- file.path(file, paste0(name, "_", gene, "_", mutation, "_mutant"))
  if (file.exists(wildtype) == FALSE | file.exists(mutant) == FALSE) next
  
  # Read in files
  wildtype <- read.delim(wildtype, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
  wildtype <- abs(wildtype[complete.cases(wildtype), ])
  wildtype_reads <- c(wildtype_reads, wildtype)
  
  mutant <- read.delim(mutant, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
  mutant <- abs(mutant[complete.cases(mutant), ])
  mutant_reads <- c(mutant_reads, mutant)
}
  
### Loop over somatic mutations and get fragment lengths
somatic_wildtype <- vector()
somatic_mutant <- vector()

for(i in 1:nrow(somatic_mutations)) {
  data <- somatic_mutations[i, ]
  sample <- data[, 1]
  gene <- "TP53"
  mutation <- data[, 2]
  
  file <- list.files(path, pattern = sample, full.names = TRUE)
  name <- list.files(path, pattern = sample, full.names = FALSE)
  wildtype <- file.path(file, paste0(name, "_", gene, "_", mutation, "_wildtype"))
  mutant <- file.path(file, paste0(name, "_", gene, "_", mutation, "_mutant"))
  if (file.exists(wildtype) == FALSE | file.exists(mutant) == FALSE) next
  
  # Read in files
  wildtype <- read.delim(wildtype, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
  wildtype <- abs(wildtype[complete.cases(wildtype), ])
  somatic_wildtype <- c(somatic_wildtype, wildtype)
  
  mutant <- read.delim(mutant, colClasses = c(rep("NULL", 8), rep("integer", 1), rep("NULL", 3)), header = FALSE)
  mutant <- abs(mutant[complete.cases(mutant), ])
  somatic_mutant <- c(somatic_mutant, mutant)
}

### Calculate statistics for plotting (germline)
wildtype_reads <- as.data.frame(wildtype_reads)
mutant_reads <- as.data.frame(mutant_reads)

wt <- nrow(wildtype_reads)
mut <- nrow(mutant_reads)
median1 <- median(wildtype_reads$wildtype_reads)
median2 <- median(mutant_reads$mutant_reads)
sd1 <- sd(wildtype_reads$wildtype_reads)
sd2 <- sd(mutant_reads$mutant_reads)
t_test <- (median1 - median2)/sqrt(((sd1^2)/wt) + ((sd2^2)/mut))
t_test <- pnorm(t_test, mean = 0, sd = 1, lower.tail = FALSE)

p <- format(round(t_test, 4), nsmall = 4)

# Plot cumulative distributions (germline)
plot <- ggplot() + 
  stat_ecdf(data = wildtype_reads, aes(wildtype_reads, color = "black"), geom = "step") +
  stat_ecdf(data = mutant_reads, aes(mutant_reads, color = "red"), geom = "step") +
  xlab("Fragment Length") + 
  ylab("Cumulative Distribution") +
  ggtitle("TP53 (Germline)") + 
  scale_colour_manual(name = 'Fragments', values =c("black" = "black", "red" = "red"), labels = c("Wildtype", "Mutant")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15), 
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
  scale_x_continuous(limits=c(50, 275, expand = c(0,0))) +
  geom_vline(xintercept = median1, linetype = "dashed", color = "black", size=0.5) +
  geom_vline(xintercept = median2, linetype = "dashed", color = "red", size=0.5) +
  annotate("text", x = 215, y = c(0.40, 0.27, 0.14), 
           label = c(paste0("p = ", p), paste0("Wildtype = ", wt), paste0("Mutant = ", mut)),
           col = c("black", "black", "red"),
           size = 4)
plot

ggsave(file.path(outdir, "Fragment_mutation_germline.pdf"), plot, device = "pdf", width = 4, height = 2.5, units = "in")

### Calculate statistics for plotting (somatic)
somatic_wildtype <- as.data.frame(somatic_wildtype)
somatic_mutant <- as.data.frame(somatic_mutant)
somatic_mutant <- somatic_mutant[!(somatic_mutant$somatic_mutant == 118), ]
somatic_mutant <- as.data.frame(somatic_mutant)

wt <- nrow(somatic_wildtype)
mut <- nrow(somatic_mutant)
median1 <- median(somatic_wildtype$somatic_wildtype)
median2 <- median(somatic_mutant$somatic_mutant)
sd1 <- sd(somatic_wildtype$somatic_wildtype)
sd2 <- sd(somatic_mutant$somatic_mutant)
t_test <- (median1 - median2)/sqrt(((sd1^2)/wt) + ((sd2^2)/mut))
t_test <- pnorm(t_test, mean = 0, sd = 1, lower.tail = FALSE)

p <- format(round(t_test, 4), nsmall = 4)

# Plot cumulative distributions (germline)
plot <- ggplot() + 
  stat_ecdf(data = somatic_wildtype, aes(somatic_wildtype, color = "black"), geom = "step") +
  stat_ecdf(data = somatic_mutant, aes(somatic_mutant, color = "red"), geom = "step") +
  xlab("Fragment Length") + 
  ylab("Cumulative Distribution") +
  ggtitle("TP53 (Somatic)") + 
  scale_colour_manual(name = 'Fragments', values =c("black" = "black", "red" = "red"), labels = c("Wildtype", "Mutant")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15), 
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
  scale_x_continuous(limits=c(50, 275, expand = c(0,0))) +
  geom_vline(xintercept = median1, linetype = "dashed", color = "black", size=0.5) +
  geom_vline(xintercept = median2, linetype = "dashed", color = "red", size=0.5) +
  annotate("text", x = 215, y = c(0.40, 0.27, 0.14), 
           label = c(paste0("p = ", p), paste0("Wildtype = ", wt), paste0("Mutant = ", mut)),
           col = c("black", "black", "red"),
           size = 4)
plot

ggsave(file.path(outdir, "Fragment_mutation_somatic.pdf"), plot, device = "pdf", width = 4, height = 2.5, units = "in")

