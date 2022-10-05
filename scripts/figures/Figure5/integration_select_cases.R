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

## Set samples to graph and order factors
samples$cancer_status <- factor(samples$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
samples$Age <- factor(samples$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
samples <- samples[order(samples$sample_parent,
                         samples$timepoint), ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$TP53_somatic <- ifelse(is.na(samples$TP53_somatic), "NS", samples$TP53_somatic)

### Samples to plot
cases <- c("LFS2", "LFS3", "LFS5", "LFS45", "LFS59", "LFS68", "LFS78")

### Loop over cases
for (case in cases) {
  
  ### Subset samples
  samples_case <- samples[samples$ext_ID == case, ]
  
  ### Make mutation table
  data_mutation <- as.matrix(samples_case[, c("TP53_somatic")])
  colnames(data_mutation) <- "TP53 Alteration"
  row.names(data_mutation) <- samples_case$sWGS

  ### Make ichorCNA table
  data_ichorCNA <- as.matrix(samples_case$ichor)
  row.names(data_ichorCNA) <- samples_case$sWGS
  data_ichorCNA <- cbind(data_ichorCNA, c(rep(0.03, nrow(data_ichorCNA))))

  ### Make fragment table
  data_fragment <- as.matrix(samples_case[ , c("zscores1", "limit1")])
  row.names(data_fragment) <- samples_case$sWGS
  
  ### Make breast score
  data_breast <- as.matrix(cbind(samples_case$score.x, breast_limit))
  row.names(data_breast) <- samples_case$sWGS
  
  ### Make all score
  data_all <- as.matrix(cbind(samples_case$score.y, all_limit))
  row.names(data_all) <- samples_case$sWGS

  ## Set Age
  data_age <- as.matrix(samples_case$Age)
  row.names(data_age) <- samples_case$sWGS

  ## Set Cancer status
  data_cancer <- as.matrix(samples_case$cancer_status)
  row.names(data_cancer) <- samples_case$sWGS

  ## Set Previous cancers
  data_history <- as.matrix(samples_case$previous_cancer)
  row.names(data_history) <- samples_case$sWGS

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
  bottom_annotation <- HeatmapAnnotation("Cancer Status" = data_cancer,
                                         #"Patient Type" = data_age,
                                         gp = gpar(col = "black"),
                                         annotation_name_side = "left",
                                         show_annotation_name = TRUE,
                                         annotation_name_rot = 0,
                                         annotation_name_gp = gpar(fontsize = 10),
                                         col = list("Patient Type" = col_age,
                                                    "Cancer Status" = col_cancer),
                                         show_legend = FALSE,
                                         simple_anno_size = unit(0.4, "cm"))

  top_annotation <- HeatmapAnnotation("Breast Methylation\nScore" = anno_lines(data_breast,
                                                                               add_points = TRUE,
                                                                               pch = c(16, NA),
                                                                               gp = gpar(col = c("black", "red"),
                                                                                         lty = c("solid", "dashed")),
                                                                               pt_gp = gpar(col = c("black", NA)),
                                                                               size = unit(2, "mm"),
                                                                               axis_param = list(side = "left",
                                                                                                 labels_rot = 0,
                                                                                                 gp = gpar(fontsize = 8)),
                                                                               
                                                                               border = FALSE),
                                      "Pan-Cancer Methylation\nScore" = anno_lines(data_all,
                                                                               add_points = TRUE,
                                                                               pch = c(16, NA),
                                                                               gp = gpar(col = c("black", "red"),
                                                                                         lty = c("solid", "dashed")),
                                                                               pt_gp = gpar(col = c("black", NA)),
                                                                               size = unit(2, "mm"),
                                                                               axis_param = list(side = "left",
                                                                                                 labels_rot = 0,
                                                                                                 gp = gpar(fontsize = 8)),
                                                                               border = FALSE),
                                      "Fragment Size\n(Z-score)" = anno_lines(data_fragment,
                                                                              add_points = TRUE,
                                                                              pch = c(16, NA),
                                                                              gp = gpar(col = c("black", "red"),
                                                                                        lty = c("solid", "dashed")),
                                                                              pt_gp = gpar(col = c("black", NA)),
                                                                              size = unit(2, "mm"),
                                                                              axis_param = list(side = "left",
                                                                                                labels_rot = 0,
                                                                                                gp = gpar(fontsize = 8)),
                                                                              border = FALSE),
                                      "Copy Number\n(Tumor Fraction)" = anno_lines(data_ichorCNA,
                                                                               add_points = TRUE,
                                                                               pch = c(16, NA),
                                                                               gp = gpar(col = c("black", "red"),
                                                                                         lty = c("solid", "dashed")),
                                                                               pt_gp = gpar(col = c("black", NA)),
                                                                               size = unit(2, "mm"),
                                                                               axis_param = list(side = "left",
                                                                                                 labels_rot = 0,
                                                                                                 gp = gpar(fontsize = 8)),
                                                                               border = FALSE),
                                  
                                  gap = unit(2, "mm"),
                                  border = TRUE,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "left",
                                  annotation_name_rot = 0,
                                  annotation_name_gp = gpar(fontsize = 10),
                                  show_legend = FALSE,
                                  simple_anno_size = unit(0.4, "cm"),
                                  height = unit(6, "cm"))

  ## Set orders
  data_mutation <- t(data_mutation)
  row_order <- row.names(data_mutation)
  col_order <- colnames(data_mutation)

  ## Generate heatmap
  pdf(file.path(outdir, paste0(case, "_integration.pdf")), width = ncol(data_mutation)/2 + 2, height = 3.5)
  oncoPrint <- oncoPrint(data_mutation,
                         alter_fun = alter_fun, 
                         col = col,
                         show_heatmap_legend = FALSE,
                         show_pct = FALSE,
                         row_order = row_order,
                         column_order = col_order,
                         top_annotation = top_annotation,
                         bottom_annotation = bottom_annotation,
                         left_annotation = NULL,
                         right_annotation = NULL,
                         show_row_names = TRUE,
                         show_column_names = TRUE,
                         row_names_gp = gpar(fontsize = 10),
                         row_names_side = "left",
                         column_title = case,
                         column_title_gp = gpar(fontsize = 12),
                         column_names_gp = gpar(fontsize = 0),
                         border = FALSE,
                         border_gp = gpar(col = "black"),
                         height = unit(0.4, "cm"),
                         width = unit(ncol(data_mutation)*1.25, "cm"))
  draw(oncoPrint, show_annotation_legend = FALSE)
  dev.off()
}