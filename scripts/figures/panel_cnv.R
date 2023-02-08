library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(data.table)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/panel_cnv"
hbc_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/panel_cnv"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/panel_cnv"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
samples_bc <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list_TS_BC.txt"

### Read in files
samples <- read.delim(samples)
samples_bc <- read.delim(samples_bc)
tp53 <- read.delim(list.files(path, "TP53", full.names = TRUE))
cnmops <- read.delim(list.files(path, "CHARM_LFS_panelCNmops.txt", full.names = TRUE))
cnmops_bc <- read.delim(list.files(path, "CHARM_LFS_panelCNmops_buffycoat.txt", full.names = TRUE))
cnmops_hbc <- read.delim(list.files(hbc_path, "panelCNmops", full.names = TRUE))
viscap <- read.delim(list.files(path, "CHARM_LFS_VisCap.txt", full.names = TRUE))
viscap_bc <- read.delim(list.files(path, "CHARM_LFS_VisCap_buffycoat.txt", full.names = TRUE))
viscap_hbc <- read.delim(list.files(hbc_path, "VisCap", full.names = TRUE))

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0035_Cf_U_PE_370_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS", "TGL49_0025_Cf_U_PE_334_TS")
samples <- samples[!(samples$TS) %in% exclude, ]
samples_bc <- samples_bc[!(samples_bc$TS == "TGL49_0213_Ly_R_PE_311_TS"), ]

### Format TP53 file
tp53 <- tp53[order(tp53$Start, decreasing = TRUE), ]
row.names(tp53) <- paste0(tp53$Chromosome, "_", tp53$Start)

### Format cnmops data and merge
cnmops <- cnmops[cnmops$Chromosome == "chr17" & cnmops$Gene == "TP53", ]
row.names(cnmops) <- paste0(cnmops$Chromosome, "_", cnmops$Start)
cnmops <- cnmops[row.names(tp53), colnames(cnmops) %in% samples$TS]
samples_cnmops <- samples[samples$TS %in% colnames(cnmops), ]

cnmops_bc <- cnmops_bc[cnmops_bc$Chromosome == "chr17" & cnmops_bc$Gene == "TP53", ]
row.names(cnmops_bc) <- paste0(cnmops_bc$Chromosome, "_", cnmops_bc$Start)
cnmops_bc <- cnmops_bc[row.names(tp53), colnames(cnmops_bc) %in% samples_bc$TS]
samples_cnmops_bc <- samples_bc[samples_bc$TS %in% colnames(cnmops_bc), ]

cnmops_hbc <- cnmops_hbc[cnmops_hbc$Chromosome == "chr17" & cnmops_hbc$Gene == "TP53", ]
row.names(cnmops_hbc) <- paste0(cnmops_hbc$Chromosome, "_", cnmops_hbc$Start)
cnmops_hbc <- cnmops_hbc[row.names(tp53), !(colnames(cnmops_hbc) %in% c("Chromosome", "Gene", "Exon", "Start", "End"))]

### Format viscap data
viscap <- viscap[viscap$Chromosome == "chr17" & viscap$Start < 8000000, ]
row.names(viscap) <- paste0(viscap$Chromosome, "_", viscap$Start - 1)
viscap <- viscap[row.names(tp53), colnames(viscap) %in% samples$TS]
samples_viscap <- samples[samples$TS %in% colnames(viscap), ]

viscap_bc <- viscap_bc[viscap_bc$Chromosome == "chr17" & viscap_bc$Start < 8000000, ]
row.names(viscap_bc) <- paste0(viscap_bc$Chromosome, "_", viscap_bc$Start - 1)
viscap_bc <- viscap_bc[row.names(tp53), colnames(viscap_bc) %in% samples_bc$TS]
samples_viscap_bc <- samples_bc[samples_bc$TS %in% colnames(viscap_bc), ]

viscap_hbc <- viscap_hbc[viscap_hbc$Chromosome == "chr17" & viscap_hbc$Start < 8000000, ]
row.names(viscap_hbc) <- paste0(viscap_hbc$Chromosome, "_", viscap_hbc$Start - 1)
viscap_hbc <- viscap_hbc[row.names(tp53), !(colnames(viscap_hbc) %in% c("Chromosome", "Start", "End"))]

### Check samples are same and in same order
samples <- samples[samples$TS %in% samples_cnmops$TS &
                     samples$TS %in% samples_viscap$TS, ]
samples$mutation_type <- ifelse(samples$germline_mutation %like% " del", "Deletion",
                                ifelse(samples$germline_mutation %like% "dup", "Duplication", "Mutation"))
samples$cancer_status <- factor(samples$cancer_status, levels = c("negative", "positive"),
                                labels = c("Negative", "Positive"))
samples <- samples[order(factor(samples$mutation_type, levels = c("Deletion", "Duplication", "Mutation")),
                         factor(samples$ext_ID, levels = paste0("LFS", c(1:92))),
                         samples$timepoint), ]

cnmops <- cnmops[ , samples$TS]
cnmops <- as.matrix(cnmops)
cnmops_bc <- as.matrix(cnmops_bc)
viscap <- viscap[ , samples$TS]
viscap <- as.matrix(viscap)
viscap_bc <- as.matrix(viscap_bc)

cnmops_hbc <- cnmops_hbc[, colnames(viscap_hbc)]
cnmops_hbc <- as.matrix(cnmops_hbc)
viscap_hbc <- as.matrix(viscap_hbc)

### Normalize cnmops
median_cn <- cnmops[, colnames(cnmops) %in% samples$TS[samples$cancer_status == "Negative" & samples$previous_cancer == "no" & samples$mutation_type == "Mutation"]]
sd_cn <- rowSds(median_cn)
median_cn <- rowMeans2(median_cn)

cnmops <- (cnmops - median_cn)/sd_cn
cnmops <- t(cnmops)

median_cn <- rowMeans2(cnmops_hbc)
sd_cn <- rowSds(cnmops_hbc)

cnmops_hbc <- (cnmops_hbc - median_cn)/sd_cn
cnmops_hbc <- t(cnmops_hbc)

median_cn <- rowMeans2(cnmops_bc)
sd_cn <- rowSds(cnmops_bc)

cnmops_bc <- (cnmops_bc - median_cn)/sd_cn
cnmops_bc <- t(cnmops_bc)

### Normalize viscap to healthy controls
median_vis <- viscap[, colnames(viscap) %in% samples$TS[samples$cancer_status == "Negative" & samples$previous_cancer == "no" & samples$mutation_type == "Mutation"]]
sd_vis <- rowSds(median_vis)
median_vis <- rowMeans2(median_vis)

viscap <- (viscap - median_vis)/sd_vis
viscap <- t(viscap)

median_vis <- rowMeans2(viscap_hbc)
sd_vis <- rowSds(viscap_hbc)

viscap_hbc <- (viscap_hbc - median_vis)/sd_vis
viscap_hbc <- t(viscap_hbc)

median_vis <- rowMeans2(viscap_bc)
sd_vis <- rowSds(viscap_bc)

viscap_bc <- (viscap_bc - median_cn)/sd_cn
viscap_bc <- t(viscap_bc)

### Score cnmops for somatic deletions
cnmops_norm <- rowSums2(cnmops_hbc)
cnmops_mean <- mean(cnmops_norm)
cnmops_sd <- sd(cnmops_norm)

cnmops_score <- rowSums2(cnmops)
cnmops_score <- (cnmops_score - cnmops_mean)/cnmops_sd

### Score viscap for deletions
viscap_norm <- rowSums2(viscap_hbc)
viscap_mean <- mean(viscap_norm)
viscap_sd <- sd(viscap_norm)

viscap_score <- rowSums2(viscap)
viscap_score <- (viscap_score - viscap_mean)/viscap_sd

### Make score table for whole gene
cnmops_score <- as.matrix(cnmops_score)
row.names(cnmops_score) <- row.names(cnmops)

viscap_score <- as.matrix(viscap_score)
row.names(viscap_score) <- row.names(viscap)

scores <- merge(cnmops_score, viscap_score, by = "row.names")
colnames(scores) <- c("sample", "cnmops", "viscap")

#scores$cnmops <- ifelse(scores$cnmops > 1.96 | scores$cnmops < -1.96, scores$cnmops, NA)
#scores$viscap <- ifelse(scores$viscap > 1.96 | scores$viscap < -1.96, scores$viscap, NA)

### Merge cnmops data
cnmops <- rbind(cnmops, cnmops_bc, cnmops_hbc)
viscap <- rbind(viscap, viscap_bc, viscap_hbc)

empty <- data.frame(matrix(nrow = nrow(samples_cnmops_bc), ncol = ncol(samples)))
colnames(empty) <- colnames(samples)
empty$TS <- samples_cnmops_bc$TS
empty$ext_ID <- "Lymphocytes"
empty$cancer_status <- "Lymphocytes"
empty$mutation_type <- "Lymphocytes"
samples <- rbind(samples, empty)

empty <- data.frame(matrix(nrow = nrow(cnmops_hbc), ncol = ncol(samples)))
colnames(empty) <- colnames(samples)
empty$TS <- rownames(cnmops_hbc)
empty$ext_ID <- "Healthy\nControls"
empty$cancer_status <- "Healthy\nControls"
empty$mutation_type <- "Healthy\nControls"
samples <- rbind(samples, empty)

samples <- samples[order(factor(samples$mutation_type, levels = c("Deletion", "Duplication", "Mutation", "Lymphocytes", "Healthy\nControls")),
                         samples$ext_ID,
                         samples$timepoint), ]
cnmops <- cnmops[samples$TS, ]
viscap <- viscap[samples$TS, ]

## Set Cancer status
data_cancer <- as.matrix(samples$cancer_status)
row.names(data_cancer) <- samples$TS

## Mutation Type
data_mutation <- as.matrix(samples$mutation_type)
row.names(data_mutation) <- samples$TS

## Set TP53 annotation
data_tp53 <- as.matrix(tp53$type)
row.names(data_tp53) <- row.names(tp53)

## Set colours
col_fun <- colorRamp2(c(-3, -2, 0, 2, 3), 
                      c("#1f78b4", "white", "white", "white", "#e31a1c"))
col_cancer <- c(Positive = "#fb9a99", Negative = "#a6cee3", Lymphocytes = "#FDBF6F", "Healthy\nControls" = "#B2DF8A")
col_mutation <- c(Duplication = "#fb9a99", Deletion = "#a6cee3", Mutation = "grey65", Lymphocytes = "#FDBF6F", "Healthy\nControls" = "#B2DF8A")
col_tp53 <- c(exon = "grey65", intron = "grey95")

## Set additional annotations
left_annotation <- rowAnnotation("Germline Variant" = data_mutation,
                                 "Cancer Status" = data_cancer,
                                 border = FALSE,
                                 annotation_name_side = "top",
                                 show_annotation_name = FALSE,
                                 col = list("Germline Variant" = col_mutation,
                                            "Cancer Status" = col_cancer),
                                 show_legend = FALSE,
                                 simple_anno_size = unit(0.3, "cm"))

bottom_annotation <- HeatmapAnnotation("TP53" = data_tp53,
                                       border = FALSE,
                                       annotation_name_side = "left",
                                       show_annotation_name = TRUE,
                                       annotation_name_gp = gpar(fontface = "italic",
                                                                 fontsize = 8),
                                       col = list("TP53" = col_tp53),
                                       show_legend = FALSE,
                                       simple_anno_size = unit(0.3, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Germline Variant", 
                                                  at = c("Deletion", "Duplication", "Small Sequence Mutation"),
                                                  legend_gp = gpar(fill = c("#a6cee3", "#fb9a99", "grey65")),
                                                  #nrow = 1,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Status", 
                                                  at = c("Negative", "Positive", "Lymphocytes", "HBC"),
                                                  legend_gp = gpar(fill = c("#a6cee3", "#fb9a99", "#FDBF6F", "#B2DF8A")),
                                                  nrow = 2,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "TP53", 
                                                  at = c("Exon", "Intron"),
                                                  legend_gp = gpar(fill = c("grey65", "grey95")),
                                                  #nrow = 1,
                                                  title_gp = gpar(fontface = "italic",
                                                                  fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Z-score",
                                                  at = c(-3, 0, 3),
                                                  col_fun = col_fun,
                                                  direction = "horizontal",
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(2, "mm"))),
                               direction = "horizontal")

## Set orders
row_order <- samples$TS
col_order_cn <- colnames(cnmops)
col_order_vis <- colnames(viscap)

## Set splits
row_split <- samples$ext_ID
order <- unique(row_split)
row_split <- factor(row_split, levels = order)

col_split <- tp53$number
order <- unique(col_split)
col_split <- factor(col_split, levels = order,
                    labels = c("1", " ", "2", "3", "4", "  ", "5", "6", "7", "8", "9", "   ", "10", "    ", "11"))

### Set white line function
white.line <- function(j, i, x, y, w, h, fill) { grid.lines(x = c(x - w/2, x + w / 2), y = c(y + h / 2, y + h / 2), gp = gpar(col = 'white', lwd = 0.5)) }

## Generate heatmap
heatmap_cn <- Heatmap(cnmops,
                      cell_fun = white.line,
                      col = col_fun,
                      show_heatmap_legend = FALSE,
                      row_order = row_order,
                      column_order = col_order_cn,
                      left_annotation = left_annotation,
                      bottom_annotation = bottom_annotation,
                      row_split = row_split,
                      row_gap = unit(0.5, "mm"),
                      column_split = col_split,
                      row_title_rot = 0,
                      row_title_side = "left",
                      row_title_gp = gpar(fontsize = 4),
                      show_row_names = FALSE,
                      column_title_side = "bottom",
                      column_title_gp = gpar(fontsize = 8),
                      show_column_names = FALSE,
                      border = FALSE)
#heatmap_cn

heatmap_vis <- Heatmap(viscap,
                       cell_fun = white.line,
                       col = col_fun,
                       show_heatmap_legend = FALSE,
                       row_order = row_order,
                       column_order = col_order_cn,
                       left_annotation = left_annotation,
                       bottom_annotation = bottom_annotation,
                       row_split = row_split,
                       row_gap = unit(0.5, "mm"),
                       column_split = col_split,
                       row_title_rot = 0,
                       show_row_names = FALSE,
                       column_title_side = "bottom",
                       column_title_gp = gpar(fontsize = 8),
                       row_title_gp = gpar(fontsize = 0),
                       show_column_names = FALSE,
                       border = FALSE)
#heatmap_vis

### Combine
heatmap <- heatmap_cn + heatmap_vis

pdf(file.path(outdir, "combined_2.pdf"), width = 5, height = 7)
draw(heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE,
     column_title = "PanelCNmops                                     VisCap", 
     column_title_gp = gpar(fontsize = 10), heatmap_legend_side = "bottom")
dev.off()

write.table(scores, file.path(outdir, "panel_cnv_whole_gene.txt"), sep = "\t", row.names = FALSE)

### Calculate scores for focal CNVs
cnmops <- ifelse(cnmops > 1.96, 1,
                 ifelse(cnmops < -1.96, -1, 0))
cnmops_norm <- cnmops[row.names(cnmops) %in% samples$TS[samples$cancer_status == "Negative" & samples$previous_cancer == "no" & samples$mutation_type == "Mutation"], ]
cnmops_sum <- rowSums2(cnmops_norm)
cnmops_upper <- quantile(cnmops_sum, 0.95)
cnmops_lower <- quantile(cnmops_sum, 0.05)
cnmops_score <- rowSums2(cnmops)

viscap <- ifelse(viscap > 1.96, 1,
                 ifelse(viscap < -1.96, -1, 0))
viscap_norm <- viscap[row.names(viscap) %in% samples$TS[samples$cancer_status == "Negative" & samples$previous_cancer == "no" & samples$mutation_type == "Mutation"], ]
viscap_sum <- rowSums2(viscap_norm)
viscap_upper <- quantile(viscap_sum, 0.95)
viscap_lower <- quantile(viscap_sum, 0.05)
viscap_score <- rowSums2(viscap)

### Make score table
cnmops_score <- as.matrix(cnmops_score)
row.names(cnmops_score) <- row.names(cnmops)
viscap_score <- as.matrix(viscap_score)
row.names(viscap_score) <- row.names(viscap)

scores <- merge(cnmops_score, viscap_score, by = "row.names")
colnames(scores) <- c("sample", "cnmops", "viscap")

scores$cnmops <- ifelse(scores$cnmops > cnmops_upper | scores$cnmops < cnmops_lower, scores$cnmops, NA)
scores$viscap <- ifelse(scores$viscap > viscap_upper | scores$viscap < viscap_lower, scores$viscap, NA)

write.table(scores, file.path(outdir, "panel_cnv_partial_gene.txt"), sep = "\t", row.names = FALSE)
