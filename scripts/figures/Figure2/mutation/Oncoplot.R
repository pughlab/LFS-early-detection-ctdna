library(ComplexHeatmap)
library(dplyr)
library(stringr)

path <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

## Import data
data_onco <- read.delim(file.path(path, "Oncoplot_full.txt"))
data_samples <- read.delim(samples)

## Remove failed samples and format sample sheet
failed <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0010_Cf_U_PE_334_TS", "TGL49_0209_Cf_U_PE_378_TS")
data_onco <- data_onco[!(data_onco$sample_ID %in% failed), ]
data_onco <- data_onco[data_onco$TP53_germline != "failed", ]
data_samples <- data_samples[data_samples$TS %in% data_onco$sample_ID, ]

## Set clinical information and sample order
data_samples$cancer_status <- factor(data_samples$cancer_status, levels = c("positive", "negative", "unknown"),
                                  labels = c("Cancer Positive", "Cancer Negative", "Unknown"))
data_samples <- data_samples[order(factor(data_samples$cancer_status, levels = c("Cancer Positive", "Cancer Negative", "Unknown"))), ]

## Set Age
data_samples$Age <- factor(data_samples$Age, levels = c("adult", "pediatric"),
                        labels = c("Adult", "Pediatric"))
data_age <- as.matrix(data_samples$Age)
row.names(data_age) <- data_samples$TS

## Set cancer status
data_cancer <- as.matrix(data_samples$cancer_status)
row.names(data_cancer) <- data_samples$TS

## Set cancer type
data_samples$cancer_type <- factor(data_samples$cancer_type, levels = c("appendiceal_adenocarcinoma", "bladder", "breast", "endometrial", "lung",
                                                                        "lymphoma", "osteosarcoma", "prostate", "sarcoma", ""),
                                labels = c("Adenocarcinoma", "Bladder", "Breast", "Endometrial", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma", "NA"))
data_type <- as.matrix(data_samples$cancer_type)
row.names(data_type) <- data_samples$TS

## Set cancer previous
data_samples$previous_cancer <- factor(data_samples$previous_cancer, levels = c("yes", "no", ""),
                                      labels = c("Yes", "No", "NA"))
data_previous <- as.matrix(data_samples$previous_cancer)
row.names(data_previous) <- data_samples$TS

## Set Stage
data_samples$stage <- factor(data_samples$stage, levels = c("low", "high", "unknown", ""),
                             labels = c("Stage 0/I/II", "Stage III/IV", "NA", "NA"))
data_stage <- as.matrix(data_samples$stage)
row.names(data_stage) <- data_samples$TS

## Format heatmap
data_onco[is.na(data_onco)] <- ""
row.names(data_onco) <- data_onco$sample_ID
data_onco <- data_onco[, -1]
data_onco <- as.data.frame(t(data_onco))
data_onco <- data_onco[!(row.names(data_onco) == "vaf"), ]
row.names(data_onco) <- c("TP53 (Germline)", "TP53 (Somatic)", "PALB2", "MSH2", "MSH6", "APC", "BRCA1", "BRCA2", "MLH1", "EPCAM", "PMS2")
data_onco <- as.matrix(data_onco)
data_onco <- data_onco[ , data_samples$TS]

## Make oncoprint barplot (minus GL)
data_bar <- data_onco[-1, ]
missense <- matrix(str_count(data_bar, "missense"), nrow = nrow(data_bar), ncol = ncol(data_bar))
missense <- colSums(missense)
missense2 <- matrix(str_count(data_bar, "missense2"), nrow = nrow(data_bar), ncol = ncol(data_bar))
missense2 <- colSums(missense2)
frameshift <- matrix(str_count(data_bar, "frameshift"), nrow = nrow(data_bar), ncol = ncol(data_bar))
frameshift <- colSums(frameshift)
stop <- matrix(str_count(data_bar, "stop"), nrow = nrow(data_bar), ncol = ncol(data_bar))
stop <- colSums(stop)
splice_site <- matrix(str_count(data_bar, "splice_site"), nrow = nrow(data_bar), ncol = ncol(data_bar))
splice_site <- colSums(splice_site)
translation <- matrix(str_count(data_bar, "translation"), nrow = nrow(data_bar), ncol = ncol(data_bar))
translation <- colSums(translation)
data_bar <- rbind(missense, missense2, frameshift, stop, splice_site, translation)
colnames(data_bar) <- colnames(data_onco)
data_bar <- t(data_bar)

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         splice_site = "#984EA3", translation = "#FF7F00",
         deletion = "#377EB8", duplication = "#E41A1C")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", 
              "Ewing Sarcoma" = "#ffee65", Glioma = "#beb9db", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Osteosarcoma = "#ebdc78", 
              Prostate = "#b3d4ff", Sarcoma = "grey65", Unknown = "grey95", "NA" = "grey95")
col_previous <- c(Yes = "#fb9a99", No = "grey95", "NA" = "grey95")
col_stage <- c("Stage 0/I/II" = "#ebdc78", "Stage III/IV" = "#fb9a99", "NA" = "grey95")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*0.75, h*0.9, gp = gpar(fill = "grey95", col = NA))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = "white"), just = "top")
}

## Set additional annotations
top_annotation <- HeatmapAnnotation("# of Somatic\nMutations" = anno_barplot(data_bar, 
                                                                       gp = gpar(fill = col,
                                                                                 col = NA),
                                                                       border = FALSE,
                                                                       height = unit(2, "cm"),
                                                                       axis_param = list(gp = gpar(fontsize = 8))),
                                    "Patient Type" = data_age,
                                    "Cancer History" = data_previous,
                                    "Cancer Stage" = data_stage,
                                    "Cancer Type" = data_type,
                                    col = list("Somatic Mutations" = col, "Patient Type" = col_age, "Cancer Type" = col_type, 
                                               "Cancer History" = col_previous, "Cancer Stage" = col_stage),
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp= gpar(fontsize = 10),
                                    simple_anno_size = unit(0.35, "cm"),
                                    border = FALSE,
                                    show_legend = c(FALSE, FALSE, FALSE, FALSE, FALSE))

right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(c("missense", "frameshift", "stop", "splice_site",
                                                                        "deletion", "gain", "inframe_deletion", "inframe_insertion", "translation"),
                                                                      border = TRUE, width = unit(1.5, "cm"), 
                                                                      axis_param = list(side = "top", labels_rot = 90, gp = gpar(fontsize = 8))),
                                 annotation_name_gp = gpar(fontsize = 10))
  

## Set labels
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
                                           Legend(title = "Cancer Stage",
                                                  at = c("Stage 0/I/II", "Stage III/IV"),
                                                  legend_gp = gpar(fill = c("#ebdc78", "#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Type",
                                                  at = c(c("Adenocarcinoma", "Bladder", "Breast", "Endometrial", 
                                                           "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma")),
                                                  legend_gp = gpar(fill = c(c("#fd7f6f","#b2e061","#bd7ebe","#ffb55a", 
                                                                              "#fdcce5","#8bd3c7","#ebdc78","#b3d4ff","grey65"))),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Alterations",
                                                  at = c("Missense", "Frameshift", "Stop", "Splice Site", "Deletion"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#984EA3", "#377EB8")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"))), 
                               max_height = unit(4, "in"),
                               row_gap = unit(0.25, "cm")
                               )

## Set splits
column_split <- data_cancer
column_split <- ifelse(column_split == "Unknown", "Cancer Negative", column_split)
column_split <- factor(column_split, levels = c("Cancer Positive", "Cancer Negative"))
row_order <- c("TP53 (Germline)", "TP53 (Somatic)", "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC")

## Generate oncoprint
pdf(file.path(path,"Oncoprint_full.pdf"), width = 10, height = 3.6)
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       show_heatmap_legend = FALSE,
                       row_order = row_order,
                       row_names_side = "left",
                       pct_side = "right",
                       top_annotation = top_annotation,
                       right_annotation = right_annotation,
                       row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                       column_split = column_split,
                       border = FALSE,
                       border_gp = gpar(col = "black"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE)
dev.off()
