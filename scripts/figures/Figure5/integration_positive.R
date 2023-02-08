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
mutation <- read.delim(file.path(path2, "Oncoplot_full.txt"))
fragment <- read.delim(file.path(path2, "fragment_scores.txt"))
breast <- read.delim(file.path(path2, "methylation_score_breast.tsv"))
all <- read.delim(file.path(path2, "pancancer_score_Vrba.tsv"))
samples <- read.delim("sample_list.txt")

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

### Subset Cancer positive only
samples <- samples[samples$cancer_status == "positive", ]
ichorCNA <- ichorCNA[ichorCNA$sample %in% samples$sWGS, ]
fragment <- fragment[fragment$sample %in% samples$sWGS, ]
mutation <- mutation[mutation$sample_ID %in% samples$TS, ]
breast <- breast[breast$cfMeDIP %in% samples$cfMeDIP, ]
all <- all[all$cfMeDIP %in% samples$cfMeDIP, ]

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
samples <- samples[rowSums(is.na(samples[, c(26:31)])) != 6, ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

## Set samples to graph and order factors
samples$cancer_type <- factor(samples$cancer_type, levels = c("appendiceal_adenocarcinoma", "adrenal", "astrocytoma", "bladder", "breast", "endometrial",
                                                              "ewing_sarcoma", "glioma", "lung", "lymphoma", "leukemia", "osteosarcoma", "prostate", "sarcoma", ""),
                              labels = c("Adenocarcinoma", "Adrenal", "Glioma", "Bladder", "Breast", "Endometrial",
                                         "Ewing Sarcoma", "Glioma", "Lung", "Lymphoma", "Leukemia", "Osteosarcoma", "Prostate", "Sarcoma", "Unknown"))
samples$cancer_status <- factor(samples$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
samples$Age <- factor(samples$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
samples$stage <- factor(samples$stage, levels = c("low", "high", "unknown", ""),
                        labels = c("Stage 0/I/II", "Stage III/IV", "Unknown", "Unknown"))
samples$ext_ID <- ifelse(samples$ext_ID %in% samples$ext_ID[duplicated(samples$ext_ID)], samples$ext_ID, "none")
samples$b1 <- "none"
samples$b2 <- "none"
samples <- samples[order(samples$stage,
                         samples$cancer_type), ]

### Make detection table
samples$ichor <- ifelse(samples$ichor > 0.03, "pos", "neg")
samples$ichor[is.na(samples$ichor)] <- "NS"

samples$limit1[is.na(samples$limit1)] <- unique(samples$limit1)[[2]]
samples$zscores1 <- ifelse(samples$zscores1 > samples$limit1, "pos", "neg")
samples$zscores1[is.na(samples$zscores1)] <- "NS"

samples$score.x <- ifelse(samples$score.x > breast_limit, "pos", "neg")
samples$score.x[is.na(samples$score.x)] <- "NS"

samples$score.y <- ifelse(samples$score.y > all_limit, "pos", "neg")
samples$score.y[is.na(samples$score.y)] <- "NS"

samples$TP53_somatic <- ifelse(is.na(samples$TP53_somatic), "NS", samples$TP53_somatic)
data_mutation <- as.matrix(samples[, c("TP53_somatic", "ichor", "zscores1", "score.x", "score.y")])
colnames(data_mutation) <- c("Mutation", "Copy Number", "Fragment Size", "Breast Methylation", "Pan-cancer Methylation")
row.names(data_mutation) <- samples$sWGS
data_mutation <- t(data_mutation)

### Set patient annotation
data_patient <- cbind(samples$b1, samples$ext_ID, samples$b2)

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65", pos = "#fb9a99", neg = "grey95")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", "Ewing Sarcoma" = "#ffee65", 
              Glioma = "#beb9db", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Leukemia = "grey65", Osteosarcoma = "#ebdc78", Prostate = "#b3d4ff", Sarcoma = "#8be04e", Unknown = "grey95")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_patient <- c(LFS3 = "grey35", LFS2 =  "grey65", LFS72 = "grey35", LFS58 = "grey65", LFS74 = "grey35", LFS79 = "grey65", LFS57 = "grey35", none = "grey95")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*1, h*1, gp = gpar(fill = "grey95", col = "black"))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = "black"), just = "top")
}

## Set additional annotations
top_annotation <- HeatmapAnnotation("Patient Type" = samples$Age,
                                 "Cancer History" = samples$previous_cancer,
                                 "Cancer Type" = samples$cancer_type,
                                 "Patient" = anno_simple(data_patient,
                                                         height = unit(0.3, "cm"),
                                                         col = col_patient),
                                 border = FALSE,
                                 annotation_name_side = "left",
                                 show_annotation_name = TRUE,
                                 annotation_name_rot = 0,
                                 annotation_name_gp = gpar(fontsize = 8),
                                 col = list("Patient Type" = col_age,
                                            "Cancer History" = col_history,
                                            "Cancer Type" = col_type,
                                            "Patient" = col_patient),
                                 show_legend = FALSE,
                                 simple_anno_size = unit(0.3, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Cancer Type",
                                                  at = c(c("Adenocarcinoma", "Adrenal", "Bladder", "Breast", "Endometrial", "Ewing Sarcoma", 
                                                           "Glioma", "Leukemia", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma")),
                                                  legend_gp = gpar(fill = c(c("#fd7f6f","#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65", 
                                                                              "#beb9db", "grey65", "#fdcce5","#8bd3c7","#ebdc78","#b3d4ff","#8be04e"))),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 7),
                                           Legend(title = "Alterations",
                                                  at = c("Missense", "Frameshift", "Stop", "Splice Site", "Deletion"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#984EA3", "#377EB8")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 5),
                                           Legend(title = "Detection",
                                                  at = c("Positive", "Negative", "Not Sequenced/\nFailed"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "grey95", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 3)),
                               max_height = unit(10, "cm"),
                               direction = "horizontal")

## Set orders
row_order <- row.names(data_mutation)
col_order <- colnames(data_mutation)
col_split <- samples$stage

## Generate heatmap
pdf(file.path(outdir, "integration_positive.pdf"), width = 6, height = 5)
oncoPrint <- oncoPrint(data_mutation,
                       alter_fun = alter_fun, 
                       col = col,
                       column_split = col_split,
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
                       column_order = col_order,
                       border = FALSE,
                       border_gp = gpar(col = "black"),
                       height = unit(0.3*nrow(data_mutation), "cm"),
                       width = unit(0.2*ncol(data_mutation), "cm"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()

