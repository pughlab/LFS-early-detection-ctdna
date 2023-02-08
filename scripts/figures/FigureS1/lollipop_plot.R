library(tidyverse)
library(trackViewer)
library(data.table)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/cohort stats"

data_samples <- read.delim(file.path(path, "samples/sample_list.txt"))

### Format sample list into list of features
data_samples <- data_samples[!(duplicated(data_samples$ext_ID)), ]
data_samples <- data_samples[, c("germline_mutation", "mutation_type")]
data_samples <- data_samples[order(data_samples$germline_mutation), ]
data_samples <- data_samples[!(data_samples$germline_mutation %like% "del" |
                                 data_samples$germline_mutation %like% "dup" |
                                 data_samples$germline_mutation == "unknown" |
                                 data_samples$germline_mutation == ""), ]

data_mutations <- as.data.frame(table(data_samples$germline_mutation))
data_mutations <- data_mutations[!(data_mutations$Var1 == "intron"), ]
data_mutations$location <- gsub("[^0-9.-]", "", data_mutations$Var1)
data_mutations$location <- ifelse(data_mutations$Var1 == "c.97-11C>G", 32, data_mutations$location)
data_mutations$type <- "missense"
data_mutations$type <- ifelse(data_mutations$Var1 %like% "fs", "frameshift",
                              ifelse(data_mutations$Var1 %like% "\\*", "stop",
                                     ifelse(data_mutations$Var1 %like% "splice", "splice",
                                            ifelse(data_mutations$Var1 %like% "c.", "splice", data_mutations$type))))

data_mutations <- data_mutations[order(data_mutations$type), ]
frameshift <- nrow(data_mutations[data_mutations$type == "frameshift", ])
missense <- nrow(data_mutations[data_mutations$type == "missense", ])
splice <- nrow(data_mutations[data_mutations$type == "splice", ])
stop <- nrow(data_mutations[data_mutations$type == "stop", ])

SNP <- data_mutations$location
SNP <- as.numeric(SNP)

### Define features
# SNP features
Patients <- GRanges("TP53", IRanges(SNP, width = 1,
                                     ## the name of GRanges object will be used as label
                                     names = data_mutations$Var1),
                     ## score value will be used to for the height of lollipop
                     score = as.numeric(data_mutations$Freq),
                     ## set the color for lollipop node.
                     color = c(rep("black", frameshift), rep("#4DAF4A", missense), rep("#984EA3", splice), rep("#A65628", stop)),
                     ## set the lollipop stem color
                     border = sample(rep("black", length(SNP)))
)

Patients$label.parameter.gp <- as.list(c(rep(gpar(col = "black"), frameshift), 
                                          rep(gpar(col = "#4DAF4A"), missense), 
                                          rep(gpar(col = "#984EA3"), splice), 
                                          rep(gpar(col = "#A65628"), stop)))

# Gene features
features.gr <- GRanges("TP53", IRanges(c(1, 1, 102, 305, 326, 364),
                                       width = c(393, 91, 190, 17, 20, 29),
                                       names = c("TP53", "Transactivation Domain", "DNA-Binding Domain", 
                                                 "Nuclear Localization Signal", "Tetramerization Domain", "Basic Domain")),
                       fill = c("grey65", "#FFFFCC", "#CCEBC5", "#B3CDE3", "#FBB4AE", "#DECBE4"), ## color for exon
                       height = c(0.02, 0.05, 0.05, 0.03, 0.05, 0.05) ## height for exon
)

# Other features
xaxis <- c(1, 100, 200, 300, 393)

legend <-list(labels = c("Frameshift", "Missense", "Splice Site", "Stop"),
              col = c("black", "#4DAF4A", "#984EA3", "#A65628"), 
              fill = c("black", "#4DAF4A", "#984EA3", "#A65628"))

### Plot lolliplot
pdf(file.path(outdir, "TP53_mutations.pdf"), height = 5, width = 7)
lolliplot(Patients, 
          features.gr, 
          xaxis = xaxis,
          legend = legend)
dev.off()



