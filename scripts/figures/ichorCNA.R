library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA"
path2 <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/ichorCNA"

### Find files
mutations <- read.delim(file.path(path2, "Oncoplot_full.txt"))
ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", full.names = TRUE))
segs <- read.delim(list.files(path, "segs.txt", full.names = TRUE))
segs <- segs[complete.cases(segs), ]
samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
ichorCNA <- ichorCNA[!(ichorCNA$sample %in% exclude), ]

### Format mutations table
mutations <- mutations[mutations$sample_ID %in% ichorCNA$TS, ]
ichorCNA <- merge(ichorCNA[,1:6], samples, by.x = "sample", by.y = "sWGS")
ichorCNA <- merge(ichorCNA, mutations, by.x = "TS", by.y = "sample_ID", all = TRUE)
ichorCNA$TP53_somatic[is.na(ichorCNA$TP53_somatic)] <- "NS"

### Get chromosome information
chr <- segs[, 1:3]
row.names(chr) <- paste0(chr$chr, "_", chr$start)

## Set samples to graph and order factors
ichorCNA$sample <- as.factor(ichorCNA$sample)
ichorCNA$cancer_status <- factor(ichorCNA$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
ichorCNA$Age <- factor(ichorCNA$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
ichorCNA$previous_cancer <- factor(ichorCNA$previous_cancer, levels = c("yes", "no", ""),
                                       labels = c("Yes", "No", "Unknown"))
ichorCNA$stage <- factor(ichorCNA$stage, levels = c("high", "low", "", "unknown"),
                         labels = c("Stage III/IV", "Stage 0/I/II", "none", "Unknown"))
ichorCNA$ichorCNA <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA <- ichorCNA[order(ichorCNA$cancer_status,
                           ichorCNA$ichorCNA,
                           factor(ichorCNA$TP53_somatic, levels = c("NS", "", unique(ichorCNA$TP53_somatic)[!(unique(ichorCNA$TP53_somatic) %in% c("NS", ""))])),
                           ichorCNA$vaf), ]
row.names(ichorCNA) <- ichorCNA$sample
order <- row.names(ichorCNA)

### Format and order segs
row.names(segs) <- paste0(segs$chr, "_", segs$start)
segs <- segs[, colnames(segs) %in% ichorCNA$sample]
segs <- segs[, ichorCNA$sample]
segs <- as.matrix(segs)
segs <- t(segs)

### Make Mutation table
data_mutation <- merge(mutations[,c("sample_ID", "TP53_somatic")], ichorCNA[, c("TS", "sample")], by.x = "sample_ID", by.y = "TS", all = TRUE)
data_mutation <- data_mutation[!(is.na(data_mutation$sample)), c("TP53_somatic", "sample")]

data_mutation$TP53_somatic <- ifelse(is.na(data_mutation$TP53_somatic), "NS", data_mutation$TP53_somatic)
data_mutation <- data_mutation[order(factor(data_mutation$sample, levels = order)), ]
data_mutation <- as.matrix(data_mutation[, c("TP53_somatic")])
colnames(data_mutation) <- "Mutation"
row.names(data_mutation) <- row.names(ichorCNA)

### Make ichorCNA table and apphend mutation information
data_ichorCNA <- as.matrix(ichorCNA$ichorCNA)
row.names(data_ichorCNA) <- ichorCNA$sample
data_ichorCNA <- cbind(data_ichorCNA, c(rep(0.03, nrow(data_ichorCNA))))

mutations <- merge(mutations[,c("sample_ID", "vaf")], ichorCNA[, c("TS", "sample")], by.x = "sample_ID", by.y = "TS", all = TRUE)
mutations <- mutations[!(is.na(mutations$sample)), c("vaf", "sample")]
data_ichorCNA <- merge(data_ichorCNA, mutations, by.x = "row.names", by.y = "sample")

row.names(data_ichorCNA) <- data_ichorCNA$Row.names
data_ichorCNA <- data_ichorCNA[ichorCNA$sample, c("V1", "V2", "vaf")]

## Set Age
data_age <- as.matrix(ichorCNA$Age)
row.names(data_age) <- ichorCNA$sample

## Set Cancer status
data_cancer <- as.matrix(ichorCNA$cancer_status)
row.names(data_cancer) <- ichorCNA$sample

## Set Previous cancers
data_history <- as.matrix(ichorCNA$previous_cancer)
row.names(data_history) <- ichorCNA$sample

## Set Tumor stage
data_stage <- as.matrix(ichorCNA$stage)
row.names(data_stage) <- ichorCNA$sample

## Set colours
col_fun <- colorRamp2(c(-0.5, 0, 0.5), 
                      c("#1f78b4", "white", "#e31a1c"))
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65")
col_cancer <- c("Positive" = "#fb9a99", "Negative" = "#a6cee3")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_stage <- c("Stage III/IV" = "#fb9a99", "Stage 0/I/II" = "#ebdc78", "none" = "grey95", Unknown = "grey65")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*1, h*1, gp = gpar(fill = "grey95", col = "white", lwd = 0.5))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x - 0.5*w + 1:n/n*w, y, 1/n*w, h, 
                  gp = gpar(fill = col[names(which(v))], col = "white", lwd = 0.5), just = "right")
}

## Set additional annotations
right_annotation <- rowAnnotation("Patient Type" = data_age,
                                  "Cancer History" = data_history,
                                  "Tumor Stage" = data_stage,
                                  "Cancer Status" = data_cancer,
                                  "Tumour Fraction" = anno_lines(data_ichorCNA,
                                                                 add_points = TRUE,
                                                                 pch = c(16, NA, 4),
                                                                 gp = gpar(col = c(NA, "red", NA),
                                                                           lty = c("solid", "dashed", "solid")),
                                                                 pt_gp = gpar(col = c("black", NA, "red")),
                                                                 size = unit(1, "mm"),
                                                                 axis_param = list(side = "top",
                                                                                   labels_rot = 0,
                                                                                   gp = gpar(fontsize = 8)),
                                                                 border = FALSE),
                                  gap = unit(1, "mm"),
                                  border = FALSE,
                                  annotation_name_side = "top",
                                  show_annotation_name = TRUE,
                                  annotation_name_gp = gpar(fontsize = 8),
                                  col = list("Patient Type" = col_age,
                                             "Cancer History" = col_history,
                                             "Cancer Status" = col_cancer,
                                             "Tumor Stage" = col_stage),
                                  show_legend = FALSE,
                                  simple_anno_size = unit(0.3, "cm"),
                                  width = unit(3.5, "cm"))

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
                                           Legend(title = "Tumor Stage",
                                                  at = c("Stage 0/I/II", "Stage III/IV", "Unknown/None"),
                                                  legend_gp = gpar(fill = c("#ebdc78", "#fb9a99", "grey65")),
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
                                                  ncol = 1,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Copy Number",
                                                  at = c(-0.5, 0, 0.5),
                                                  col_fun = col_fun,
                                                  direction = "horizontal",
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(2, "mm")),
                                           Legend(title = "Scores",
                                                  type = "points",
                                                  at = c("Tumor Fraction", "Mutation VAF"),
                                                  pch = c(16, 4),
                                                  legend_gp = gpar(col = c("black", "red")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white"),
                                           Legend(title = "",
                                                  type = "lines",
                                                  at = c("IchorCNA LOD"),
                                                  legend_gp = gpar(col = c("red"),
                                                                   lty = c("dashed")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white")
                               ))

## Set chromosome order
chr_levels <- c(1:22, "X")
chr$chr <- factor(chr$chr, levels = chr_levels)
chr$chr <- factor(chr$chr, levels = chr_levels,
                       labels = c(1:17, " ", 19, "  ", 21, "   ", "X"))


## Set orders
row_order <- row.names(segs)
col_order <- colnames(segs)
row_split <- ichorCNA$cancer_status
row_split <- factor(row_split, levels = c("Negative", "Positive"))
col_split <- chr$chr

## Generate heatmap
ichor <- Heatmap(segs,
                 col = col_fun,
                 show_heatmap_legend = FALSE,
                 row_order = row_order,
                 column_order = col_order,
                 row_split = row_split,
                 column_split = col_split,
                 show_row_names = FALSE,
                 column_title_gp = gpar(fontsize = 8),
                 row_title = NULL,
                 show_column_names = FALSE,
                 border = FALSE)
ichor

muts <- oncoPrint(data_mutation,
                  alter_fun = alter_fun, 
                  col = col,
                  show_heatmap_legend = FALSE,
                  row_order = row_order,
                  top_annotation = NULL,
                  right_annotation = right_annotation,
                  row_split = row_split,
                  show_row_names = FALSE,
                  row_title_gp = gpar(fontsize = 0),
                  show_column_names = TRUE,
                  column_names_side = "top",
                  column_names_rot = 90,
                  column_names_gp = gpar(fontsize = 8),
                  show_pct = FALSE,
                  border = FALSE,
                  width = unit(0.3, "cm"))
muts
Figure <- ichor + muts
Figure

pdf(file.path(outdir, "ichorCNA.pdf"), width = 10, height = 5)
draw(Figure, heatmap_legend_list = annotation_legend, heatmap_legend_side = "left", show_annotation_legend = FALSE) 
dev.off()
