library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- ""
path2 <- ""
outdir <- ""

### Find files
ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", recursive = TRUE, full.names = TRUE))
mutation <- read.delim(file.path(path2, "oncoplot/Oncoplot_full.txt"))
tp53 <- read.delim(file.path(path, "fragment_score/CHARM_LFS_panel_score.txt"))
fragment <- read.delim(file.path(path, "fragment_score/CHARM_LFS_genome_score.txt"))
breast <- read.delim(file.path(path2, "cfMeDIP/breast_matrix.tsv"))
all <- read.delim(file.path(path2, "cfMeDIP/pancancer_score_Vrba.tsv"))
samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Calculate cumulative breast score
breast <- data.frame(One = colSums(breast),
                     cfMeDIP = colnames(breast))

### Calculate methylation limits
breast_limit <- quantile(breast$One[!(breast$cfMeDIP %in% samples$cfMeDIP)], 0.99)
all_limit <- quantile(all$score[!(all$cfMeDIP %in% samples$cfMeDIP)], 0.99)
tp53_limit <- quantile(tp53$score[tp53$sample %in% samples$TS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)
fragment_limit <- quantile(fragment$score[fragment$sample %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)

### Remove failed and unknown samples and format 
exclude_wg <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
exclude_ts <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
exclude_cm <- c("TGL49_0010_Cf_U_PE_319_CM", "TGL49_0035_Cf_U_PE_313_CM", "TGL49_0041_Cf_U_PE_311_CM", "TGL49_0209_Cf_U_PE_355_CM")

samples <- samples[!(samples$sWGS %in% exclude_wg &
                       samples$TS %in% exclude_ts &
                       samples$cfMeDIP %in% exclude_cm), ]
ichorCNA <- ichorCNA[!(ichorCNA$sample %in% exclude_wg), ]
fragment <- fragment[!(fragment$sample %in% exclude_wg), ]
mutation <- mutation[!(mutation$sample_ID %in% exclude_ts), ]
tp53 <- tp53[!(tp53$sample %in% exclude_ts), ]
breast <- breast[!(breast$cfMeDIP %in% exclude_cm), ]
all <- all[!(all$cfMeDIP %in% exclude_cm), ]

### Subset Cancer negative only
samples <- samples[samples$cancer_status == "negative", ]
ichorCNA <- ichorCNA[ichorCNA$sample %in% samples$sWGS, ]
fragment <- fragment[fragment$sample %in% samples$sWGS, ]
mutation <- mutation[mutation$sample_ID %in% samples$TS, ]
tp53 <- tp53[tp53$sample %in% samples$TS, ]
breast <- breast[breast$cfMeDIP %in% samples$cfMeDIP, ]
all <- all[all$cfMeDIP %in% samples$cfMeDIP, ]

### Merge data together
colnames(fragment) <- c("sample", "genome")
colnames(tp53) <- c("sample", "tp53")
ichorCNA$ichor <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA <- ichorCNA[, c("sample", "ichor")]
mutation <- mutation[, c("sample_ID", "TP53_somatic")]
breast <- breast[, c("cfMeDIP", "One")]
colnames(breast) <- c("cfMeDIP", "breast")
all <- all[, c("cfMeDIP", "score")]
colnames(all) <- c("cfMeDIP", "pancancer")

samples <- samples[!(samples$sWGS == "" & samples$TS == "" & samples$cfMeDIP == ""), ]
samples <- merge(samples, fragment, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, ichorCNA, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, mutation, by.x = "TS", by.y = "sample_ID", all = TRUE)
samples <- merge(samples, tp53, by.x = "TS", by.y = "sample", all = TRUE)
samples <- merge(samples, breast, by = "cfMeDIP", all = TRUE)
samples <- merge(samples, all, by = "cfMeDIP", all = TRUE)
samples <- samples[!(is.na(samples$sWGS)), ]
samples <- samples[rowSums(is.na(samples[, c(26:31)])) != 6, ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

## Set samples to graph and order factors
samples$Age <- factor(samples$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
muts <- unique(samples$TP53_somatic)
muts <- muts[!(muts %in% c("", "NS"))]

### Make detection table
samples$ichor <- ifelse(samples$ichor > 0.03, "pos", "neg")
samples$ichor[is.na(samples$ichor)] <- "NS"

samples$genome <- ifelse(samples$genome > fragment_limit, "pos", "neg")
samples$genome[is.na(samples$genome)] <- "NS"

samples$tp53 <- ifelse(samples$tp53 > tp53_limit, "pos", "neg")
samples$tp53[is.na(samples$tp53)] <- "NS"

samples$breast <- ifelse(samples$breast > breast_limit, "pos", "neg")
samples$breast[is.na(samples$breast)] <- "NS"

samples$pancancer <- ifelse(samples$pancancer > all_limit, "pos", "neg")
samples$pancancer[is.na(samples$pancancer)] <- "NS"

samples$TP53_somatic <- ifelse(is.na(samples$TP53_somatic), "NS", samples$TP53_somatic)
data_mutation <- as.data.frame(samples[, c("TP53_somatic","tp53", "ichor", "genome", "breast", "pancancer")])
row.names(data_mutation) <- samples$sWGS

data_mutation <- data_mutation[order(factor(data_mutation$TP53_somatic, levels = c(muts, "", "NS")),
                                     factor(data_mutation$tp53, levels = c("pos", "neg")),
                                     factor(data_mutation$ichor, levels = c("pos", "neg")),
                                     factor(data_mutation$genome, levels = c("pos", "neg")),
                                     factor(data_mutation$breast, levels = c("pos", "neg")),
                                     factor(data_mutation$pancancer, levels = c("pos", "neg"))), ]

colnames(data_mutation) <- c("Somatic TP53 Mutation", "TP53 Fragment Score", "Copy Number/Tumour Fraction", "Genome Fragment Score", "Breast Methylation", "Pan-cancer Methylation")
data_mutation <- t(data_mutation)

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65", pos = "#fb9a99", neg = "grey95")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*1, h*1, gp = gpar(fill = "grey95", col = "white"))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = "white"), just = "top")
}

## Set additional annotations
top_annotation <- HeatmapAnnotation("Patient Type" = samples$Age,
                                    "Cancer History" = samples$previous_cancer,
                                    border = FALSE,
                                    annotation_name_side = "left",
                                    show_annotation_name = TRUE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    col = list("Patient Type" = col_age,
                                               "Cancer History" = col_history),
                                    show_legend = FALSE,
                                    simple_anno_size = unit(0.3, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Alterations",
                                                  at = c("Missense", "Frameshift", "Stop", "Splice Site", "Deletion"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#984EA3", "#377EB8")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 1),
                                           Legend(title = "Detection",
                                                  at = c("Positive", "Negative", "Insufficient/Failed\nSequencing"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "grey95", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 1)),
                               max_height = unit(10, "cm"),
                               direction = "horizontal")

## Set orders
row_order <- row.names(data_mutation)
col_order <- colnames(data_mutation)

## Generate heatmap
pdf(file.path(outdir, "integration_negative.pdf"), width = 7, height = 2.5)
oncoPrint <- oncoPrint(data_mutation,
                       alter_fun = alter_fun, 
                       col = col,
                       column_order = col_order,
                       show_heatmap_legend = FALSE,
                       show_pct = FALSE,
                       top_annotation = top_annotation,
                       right_annotation = NULL,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 8),
                       show_column_names = FALSE,
                       row_order = row_order,
                       column_title_gp = gpar(fontsize = 8),
                       border = FALSE,
                       border_gp = gpar(col = "black"),
                       height = unit(0.3*nrow(data_mutation), "cm"),
                       width = unit(0.1*ncol(data_mutation), "cm"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()

