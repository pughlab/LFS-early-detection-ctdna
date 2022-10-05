library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(readr)
library(data.table)

### Set paths
path <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/cfMeDIP"
path2 <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures"
path3 <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"

### Find files
data_matrix <- read_tsv(file.path(path, "pancancer_matrix_Vrba.tsv"))
data_scores <- read_tsv(file.path(path, "pancancer_score_Vrba.tsv"))
data_samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

ichorCNA <- read.delim(list.files(path3, "summary_reviewed.txt", recursive = TRUE, full.names = TRUE))
mutation <- read.delim(file.path(path2, "oncoplot/Oncoplot_full.txt"))
fragment <- read.delim(file.path(path2, "fragment_frequency/fragment_scores.txt"))

### Format other analyses
ichorCNA$ichor <- ifelse(ichorCNA$TF > ichorCNA$TF_short, ichorCNA$TF, ichorCNA$TF_short)
ichorCNA$ichor <- ifelse(ichorCNA$ichor >= 0.03, "pos", "neg")
ichorCNA <- ichorCNA[, c("sample", "ichor")]

fragment$frag <- ifelse(fragment$zscores1 > fragment$limit1, "pos", "neg")
fragment <- fragment[, c("sample", "frag")]

### Reformat matrices
hbc_matrix <- data_matrix[, colnames(data_matrix) %like% "HBC"]
colnames(hbc_matrix) <- paste0("HBC_", c(1:ncol(hbc_matrix)))
colnames(hbc_matrix) <- paste0("HBC_", c(1:ncol(hbc_matrix)))

hbc_scores <- data_scores[data_scores$cancer_type == "HBC", ]
hbc_scores$cfMeDIP <- colnames(hbc_matrix)

### Remove failed samples and combine scores/samples
exclude <- c("TGL49_0010_Cf_U_PE_319_CM", "TGL49_0209_Cf_U_PE_355_CM")
data_samples <- data_samples[!(data_samples$cfMeDIP %in% exclude), ]
data_matrix <- data_matrix[, colnames(data_matrix) %in% data_samples$cfMeDIP]
data_scores <- data_scores[data_scores$cfMeDIP %in% data_samples$cfMeDIP, ]
data_scores <- merge(data_scores[, c(1,3)], data_samples, by = "cfMeDIP")
data_scores <- bind_rows(data_scores, hbc_scores)
data_scores[is.na(data_scores)] <- "HBC"

data_matrix <- bind_cols(data_matrix, hbc_matrix)

### Bind other data to cfMeDIP data
data_scores <- merge(data_scores, ichorCNA, by.x = "sWGS", by.y = "sample", all = TRUE)
data_scores <- merge(data_scores, fragment, by.x = "sWGS", by.y = "sample", all = TRUE)
data_scores <- merge(data_scores, mutation[, c("sample_ID", "TP53_somatic")], by.x = "TS", by.y = "sample_ID", all = TRUE)
data_scores <- data_scores[!(is.na(data_scores$cfMeDIP)), ]

data_scores$ichor <- ifelse(data_scores$cancer_status == "HBC", "HBC", data_scores$ichor)
data_scores$frag <- ifelse(data_scores$cancer_status == "HBC", "HBC", data_scores$frag)
data_scores$TP53_somatic <- ifelse(data_scores$cancer_status == "HBC", "HBC", data_scores$TP53_somatic)
data_scores[, c("ichor", "frag", "TP53_somatic")][is.na(data_scores[, c("ichor", "frag", "TP53_somatic")])] <- "NS"

### Factor and order tables
data_scores$cancer_type <- factor(data_scores$cancer_type, levels = c("appendiceal_adenocarcinoma", "adrenal", "astrocytoma", "bladder", "breast", "endometrial",
                                                                      "ewing_sarcoma", "glioma", "lung", "lymphoma", "osteosarcoma", "prostate", "sarcoma", "", "HBC"),
                                  labels = c("Adenocarcinoma", "Adrenal", "Glioma", "Bladder", "Breast", "Endometrial",
                                             "Sarcoma", "Glioma", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma", "Unknown", "Healthy"))
data_scores$cancer_status <- factor(data_scores$cancer_status, levels = c("HBC", "negative", "positive"),
                                    labels = c("Healthy", "Cancer Negative", "Cancer Positive"))
data_scores$Age <- factor(data_scores$Age, levels = c("HBC", "adult", "pediatric"),
                          labels = c("Healthy", "Adult", "Pediatric"))
data_scores$previous_cancer <- factor(data_scores$previous_cancer, levels = c("yes", "no", "", "HBC"),
                                      labels = c("Yes", "No", "Unknown", "Healthy"))
data_scores$stage <- factor(data_scores$stage, levels = c("HBC", "high", "low", "", "unknown"),
                            labels = c("Healthy", "Stage III/IV", "Stage 0/I/II", "none", "Unknown"))
data_scores$ichor <- factor(data_scores$ichor, levels = c("HBC", "NS", "neg", "pos"),
                            labels = c("Healthy", "Not Sequenced", "Negative", "Positive"))
data_scores$frag <- factor(data_scores$frag, levels = c("HBC", "NS", "neg", "pos"),
                            labels = c("Healthy", "Not Sequenced", "Negative", "Positive"))
data_scores$TP53_somatic <- factor(data_scores$TP53_somatic, levels = c("HBC", "NS", "", unique(data_scores$TP53_somatic[!(data_scores$TP53_somatic %in% c("HBC", "NS", ""))])),
                                   labels = c("Healthy", "NS", "", unique(data_scores$TP53_somatic[!(data_scores$TP53_somatic %in% c("HBC", "NS", ""))])))

data_scores <- data_scores[order(data_scores$cancer_status,
                           data_scores$score), ]
data_matrix <- data_matrix[, data_scores$cfMeDIP]
data_matrix <- as.matrix(data_matrix)

### Make mutation table
data_mutation <- as.matrix(data_scores[, c("TP53_somatic")])
data_mutation <- t(data_mutation)
row.names(data_mutation) <- "Mutation"
colnames(data_mutation) <- data_scores$cfMeDIP

### Calculate threshold
threshold <- quantile(hbc_scores$score, 0.95)
data_anno <- as.matrix(cbind(data_scores$score,
                             threshold))
row.names(data_anno) <- data_scores$cfMeDIP

## Set colours
col_fun <- colorRamp2(c(0, 1), 
                      c("white", "#e31a1c"))
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65", Healthy = "#B2DF8A")
col_cancer <- c(Healthy = "#B2DF8A", Positive = "#fb9a99", Negative = "#a6cee3")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", 
              Glioma = "#beb9db", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Osteosarcoma = "#ebdc78", 
              Prostate = "#b3d4ff", Sarcoma = "grey65", Unknown = "grey95", "NA" = "grey95", Healthy = "#B2DF8A")
col_stage <- c(Healthy = "#B2DF8A", "Stage III/IV" = "#fb9a99", "Stage 0/I/II" = "#ebdc78", "none" = "grey95", Unknown = "grey65")
col_other <- c(Healthy = "#B2DF8A", Positive = "#fb9a99", "Not Sequenced" = "grey65", Negative = "grey95")

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
                  gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
}

## Set additional annotations
top_annotation <- HeatmapAnnotation("Pan-Cancer\nScore" = anno_lines(data_anno,
                                                                        add_points = TRUE,
                                                                        pch = c(16, NA),
                                                                        gp = gpar(col = c(NA, "red"),
                                                                                  lty = c("solid", "dashed")),
                                                                        pt_gp = gpar(col = c("black", NA)),
                                                                        size = unit(1, "mm"),
                                                                        axis_param = list(side = "left",
                                                                                          labels_rot = 0,
                                                                                          gp = gpar(fontsize = 8)),
                                                                        border = FALSE,
                                                                        height = unit(2, "cm")),
                                    "Cancer Type" = data_scores$cancer_type,
                                    "Cancer Stage" = data_scores$stage,
                                    "Cancer History" = data_scores$previous_cancer,
                                    "Patient Type" = data_scores$Age,
                                    "Fragment Size" = data_scores$frag,
                                    "Copy Number" = data_scores$ichor,
                                    border = FALSE,
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 10),
                                    annotation_name_rot = 0,
                                    col = list("Patient Type" = col_age,
                                               "Cancer History" = col_history,
                                               "Cancer Type" = col_type,
                                               "Cancer Stage" = col_stage,
                                               "Fragment Size" = col_other,
                                               "Copy Number" = col_other),
                                    show_legend = FALSE,
                                    simple_anno_size = unit(0.3, "cm"))

## Set orders
col_order <- colnames(data_matrix)
col_split <- data_scores$cancer_status
order <- unique(col_split)
col_split <- factor(col_split, levels = order)

## Generate heatmap
Heatmap <- Heatmap(data_matrix,
                   col = col_fun,
                   show_heatmap_legend = FALSE,
                   column_order = col_order,
                   column_split = col_split,
                   row_title_rot = 0,
                   show_row_names = FALSE,
                   column_title_gp = gpar(fontsize = 8),
                   row_title_gp = gpar(fontsize = 8),
                   show_column_names = FALSE,
                   border = FALSE)

oncoPrint <- oncoPrint(data_mutation,
                       alter_fun = alter_fun, 
                       col = col,
                       show_heatmap_legend = FALSE,
                       column_order = col_order,
                       column_split = col_split,
                       show_pct = FALSE,
                       top_annotation = top_annotation,
                       right_annotation = NULL,
                       show_column_names = FALSE,
                       row_names_gp = gpar(fontsize = 10),
                       row_names_side = "left",
                       column_title_gp = gpar(fontsize = 10),
                       border = FALSE,
                       border_gp = gpar(col = "black"),
                       height = unit(0.3, "cm"))
Figure <- oncoPrint %v% Heatmap

pdf(file.path(path, "cfMeDIP_pancancer.pdf"), width = 10, height = 3.5)
draw(Figure, show_annotation_legend = FALSE) 
dev.off()
