library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
path2 <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures"
outdir <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/integration"

### Find files
ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", recursive = TRUE, full.names = TRUE))
mutation <- read.delim(file.path(path2, "oncoplot/Oncoplot_full.txt"))
fragment <- read.delim(file.path(path2, "fragment_frequency/fragment_scores.txt"))
breast <- read.delim(file.path(path2, "cfMeDIP/methylation_score_breast.tsv"))
all <- read.delim(file.path(path2, "cfMeDIP/pancancer_score_Vrba.tsv"))
samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Calculate methylation limits
breast_limit <- quantile(breast$score[breast$cfMeDIP == "HBC"], 0.95)
all_limit <- quantile(all$score[all$cfMeDIP == "HBC"], 0.95)

### Remove failed and unknown samples and format 
exclude_wg <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
exclude_ts <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0010_Cf_U_PE_334_TS", "TGL49_0209_Cf_U_PE_378_TS")
exclude_cm <- c("TGL49_0010_Cf_U_PE_319_CM", "TGL49_0209_Cf_U_PE_355_CM")

samples <- samples[!(samples$sWGS %in% exclude_wg &
                       samples$TS %in% exclude_ts &
                       samples$cfMeDIP %in% exclude_cm), ]
ichorCNA <- ichorCNA[!(ichorCNA$sample %in% exclude_wg), ]
fragment <- fragment[!(fragment$sample %in% exclude_wg), ]
mutation <- mutation[!(mutation$sample_ID %in% exclude_ts), ]
breast <- breast[!(breast$cfMeDIP %in% exclude_cm), ]
all <- all[!(all$cfMeDIP %in% exclude_cm), ]

### Merge data together
fragment <- fragment[, c("sample", "zscores1", "limit1")]
ichorCNA$ichor <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA <- ichorCNA[, c("sample", "ichor")]
mutation <- mutation[, c("sample_ID", "TP53_somatic")]
breast <- breast[, c("cfMeDIP", "score")]
all <- all[, c("cfMeDIP", "score")]

samples <- samples[!(samples$sWGS == "" & samples$TS == "" & samples$cfMeDIP == ""), ]
samples <- merge(samples, fragment, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, ichorCNA, by.x = "sWGS", by.y = "sample", all = TRUE)
samples <- merge(samples, mutation, by.x = "TS", by.y = "sample_ID", all = TRUE)
samples <- merge(samples, breast, by = "cfMeDIP", all = TRUE)
samples <- merge(samples, all, by = "cfMeDIP", all = TRUE)
samples <- samples[!(is.na(samples$sWGS)), ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

### Keep serial samples only
samples <- samples[samples$ext_ID %in% samples$ext_ID[duplicated(samples$ext_ID)], ]

## Set samples to graph and order factors
samples$cancer_status <- factor(samples$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
samples$Age <- factor(samples$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
samples <- samples[order(samples$ext_ID,
                         samples$timepoint), ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)

### Make mutation table
samples$TP53_somatic <- ifelse(is.na(samples$TP53_somatic), "NS", samples$TP53_somatic)
data_mutation <- as.matrix(samples[, c("TP53_somatic")])
colnames(data_mutation) <- "Mutation"
row.names(data_mutation) <- samples$sWGS

### Make ichorCNA table
data_ichorCNA <- as.matrix(samples$ichor)
row.names(data_ichorCNA) <- samples$sWGS
data_ichorCNA <- cbind(data_ichorCNA, c(rep(0.03, nrow(data_ichorCNA))))

### Make fragment table
data_fragment <- as.matrix(samples[ , c("zscores1", "limit1")])
row.names(data_fragment) <- samples$sWGS

### Make breast score
data_breast <- as.matrix(cbind(samples$score.x, breast_limit))
row.names(data_breast) <- samples$sWGS

### Make pan cancer score
data_all <- as.matrix(cbind(samples$score.y, all_limit))
row.names(data_all) <- samples$sWGS

## Set Age
data_age <- as.matrix(samples$Age)
row.names(data_age) <- samples$sWGS

## Set Cancer status
data_cancer <- as.matrix(samples$cancer_status)
row.names(data_cancer) <- samples$sWGS

## Set Previous cancers
data_history <- as.matrix(samples$previous_cancer)
row.names(data_history) <- samples$sWGS

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65")
col_cancer <- c(Positive = "#fb9a99", Negative = "#a6cee3")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*1, h*1, gp = gpar(fill = "grey95", col = "black"))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x - 0.5*w + 1:n/n*w, y, 1/n*w, h, 
                  gp = gpar(fill = col[names(which(v))], col = "black"), just = "right")
}

## Set additional annotations
left_annotation <- rowAnnotation("Patient Type" = data_age,
                                 "Cancer History" = data_history,
                                 "Cancer Status" = data_cancer,
                                 border = FALSE,
                                 annotation_name_side = "top",
                                 show_annotation_name = TRUE,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 7),
                                 col = list("Patient Type" = col_age,
                                            "Cancer History" = col_history,
                                            "Cancer Status" = col_cancer),
                                 show_legend = FALSE,
                                 simple_anno_size = unit(0.3, "cm"))

right_annotation <- rowAnnotation("Copy Number\n(Tumor Fraction)" = anno_lines(data_ichorCNA,
                                                             add_points = TRUE,
                                                             pch = c(16, NA),
                                                             gp = gpar(col = c("black", "red"),
                                                                       lty = c("solid", "dashed")),
                                                             pt_gp = gpar(col = c("black", NA)),
                                                             size = unit(1, "mm"),
                                                             axis_param = list(side = "top",
                                                                               labels_rot = 90,
                                                                               gp = gpar(fontsize = 6)),
                                                             border = FALSE),
                                  "Fragment Size\n(Z-score)" = anno_lines(data_fragment,
                                                               add_points = TRUE,
                                                               pch = c(16, NA),
                                                               gp = gpar(col = c("black", "red"),
                                                                         lty = c("solid", "dashed")),
                                                               pt_gp = gpar(col = c("black", NA)),
                                                               size = unit(1, "mm"),
                                                               axis_param = list(side = "top",
                                                                                 labels_rot = 90,
                                                                                 gp = gpar(fontsize = 6)),
                                                               border = FALSE),
                                  "Breast\nMethylation\nScore" = anno_lines(data_breast,
                                                                          add_points = TRUE,
                                                                          pch = c(16, NA),
                                                                          gp = gpar(col = c("black", "red"),
                                                                                    lty = c("solid", "dashed")),
                                                                          pt_gp = gpar(col = c("black", NA)),
                                                                          size = unit(1, "mm"),
                                                                          axis_param = list(side = "top",
                                                                                            labels_rot = 90,
                                                                                            gp = gpar(fontsize = 6)),
                                                                          border = FALSE),
                                  "Pan-Cancer\nMethylation\nScore" = anno_lines(data_all,
                                                                           add_points = TRUE,
                                                                           pch = c(16, NA),
                                                                           gp = gpar(col = c("black", "red"),
                                                                                     lty = c("solid", "dashed")),
                                                                           pt_gp = gpar(col = c("black", NA)),
                                                                           size = unit(1, "mm"),
                                                                           axis_param = list(side = "top",
                                                                                             labels_rot = 90,
                                                                                             gp = gpar(fontsize = 6)),
                                                                           border = FALSE),
                                  gap = unit(2, "mm"),
                                  border = FALSE,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "top",
                                  annotation_name_rot = 90,
                                  annotation_name_gp = gpar(fontsize = 7),
                                  show_legend = FALSE,
                                  simple_anno_size = unit(0.3, "cm"),
                                  width = unit(4, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Patient Type", 
                                                  at = c("Adult", "Pediatric"),
                                                  legend_gp = gpar(fill = c("#6A3D9A", "#CAB2D6")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer History",
                                                  at = c("Yes"),
                                                  legend_gp = gpar(fill = c("#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Status",
                                                  at = c("Positive", "Negative"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "#a6cee3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Mutations",
                                                  at = c("Missense", "Frameshift", "Stop", "Deletion", "Not Detected", "Not Sequenced"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#377EB8", "white", "grey65")),
                                                  nrow = 2,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")), 
                                           Legend(title = "Analyses",
                                                  type = "lines",
                                                  at = c("Detection Threshold"),
                                                  legend_gp = gpar(col = c("red"), 
                                                                   lty = c("dashed")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white")),
                               direction = "horizontal")

## Set orders
row_order <- row.names(data_mutation)
col_order <- colnames(data_mutation)
row_split <- samples$ext_ID
order <- unique(row_split)
row_split <- factor(row_split, levels = order)

## Generate heatmap
pdf(file.path(outdir, "integration.pdf"), width = 3, height = 10)
oncoPrint <- oncoPrint(data_mutation,
                       alter_fun = alter_fun, 
                       col = col,
                       show_heatmap_legend = FALSE,
                       show_pct = FALSE,
                       row_order = row_order,
                       top_annotation = NULL,
                       left_annotation = left_annotation,
                       right_annotation = right_annotation,
                       show_row_names = FALSE,
                       show_column_names = TRUE,
                       row_split = row_split,
                       row_gap = unit(0.35, "mm"),
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize = 6),
                       column_names_side = "top",
                       column_names_rot = 90,
                       column_names_gp = gpar(fontsize = 8),
                       border = FALSE,
                       border_gp = gpar(col = "black"),
                       width = unit(0.3, "cm"))
draw(oncoPrint, show_annotation_legend = FALSE)
dev.off()

pdf(file.path(outdir, "legend.pdf"), width = 8, height = 3)
draw(annotation_legend)
dev.off()
