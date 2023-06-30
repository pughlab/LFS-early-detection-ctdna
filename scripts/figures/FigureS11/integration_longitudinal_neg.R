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
integrated <- readRDS(file.path(outdir, "integrated_scores.rds"))
samples <- read.delim("sample_list.txt")

### Calculate cumulative breast score
breast <- data.frame(One = colSums(breast),
                     cfMeDIP = colnames(breast))

### Amalgamate integrated scores
integrated <- unlist(integrated, recursive = FALSE)
integrated <- bind_rows(integrated)

integrated <- integrated %>%
  group_by(sample) %>%
  dplyr::summarise(score=mean(positive),
                   sd=sd(positive),
                   N=n())
integrated$CI <- integrated$sd/sqrt(integrated$N)*qt(0.975, integrated$N -1)

### Calculate methylation limits
breast_limit <- quantile(breast$One[!(breast$cfMeDIP %in% samples$cfMeDIP)], 0.99)
all_limit <- quantile(all$score[!(all$cfMeDIP %in% samples$cfMeDIP)], 0.99)
tp53_limit <- quantile(tp53$score[tp53$sample %in% samples$TS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)
fragment_limit <- quantile(fragment$score[fragment$sample %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)
integrated_limit <- quantile(integrated$score[integrated$sample %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)

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
samples <- merge(samples, integrated[,1:2], by.x = "sWGS", by.y = "sample")
samples <- samples[!(is.na(samples$sWGS)), ]
#samples <- samples[rowSums(is.na(samples[, c(26:31)])) != 6, ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)
samples$sWGS <- ifelse(samples$sWGS == "", samples$cfMeDIP, samples$sWGS)

### Keep serial samples and phenoconverters only
samples <- samples[samples$ext_ID %in% samples$ext_ID[duplicated(samples$ext_ID)], ]
samples <- samples[!(samples$ext_ID %in% samples$ext_ID[samples$cancer_status == "positive"]), ]

## Set samples to graph and order factors
samples$cancer_status <- factor(samples$cancer_status,
                                     levels = c("negative", "positive"),
                                     labels = c("Negative", "Positive"))
samples$Age <- factor(samples$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
samples$developed <- ifelse(samples$ext_ID %in% c("LFS28", "LFS45", "LFS55", "LFS68", "LFS71"), "yes", "no")
samples <- samples[order(samples$ext_ID,
                         samples$timepoint), ]
samples$sWGS <- ifelse(samples$sWGS == "", samples$TS, samples$sWGS)

### Make mutation table
samples$TP53_somatic <- ifelse(is.na(samples$TP53_somatic), "NS", samples$TP53_somatic)
samples$tp53 <- ifelse(is.na(samples$tp53), "NS", samples$tp53)
data_mutation <- as.matrix(samples[, c("TP53_somatic")])
colnames(data_mutation) <- c("TP53 Somatic Mutation")
row.names(data_mutation) <- samples$sWGS
data_mutation <- t(data_mutation)

### Make TP53 fragment score
data_tp53 <- data.frame(score = samples$tp53,
                        limit = tp53_limit)
data_tp53$score <- ifelse(data_tp53$score == "NS", NA, data_tp53$score)
data_tp53$score <- as.numeric(data_tp53$score)
data_tp53 <- as.matrix(data_tp53)
row.names(data_tp53) <- samples$sWGS

### Make ichorCNA table
data_ichorCNA <- as.matrix(samples$ichor)
row.names(data_ichorCNA) <- samples$sWGS
data_ichorCNA <- cbind(data_ichorCNA, c(rep(0.03, nrow(data_ichorCNA))))

### Make fragment table
data_fragment <- data.frame(score = samples$genome,
                            limit = fragment_limit)
data_fragment <- as.matrix(data_fragment)
row.names(data_fragment) <- samples$sWGS

### Make breast score
data_breast <- as.matrix(cbind(samples$breast, breast_limit))
row.names(data_breast) <- samples$sWGS

### Make pan cancer score
data_all <- as.matrix(cbind(samples$pancancer, all_limit))
row.names(data_all) <- samples$sWGS

### Make integrated score
data_int <- as.matrix(cbind(samples$score, integrated_limit))
row.names(data_int) <- samples$sWGS

## Set Age
data_age <- as.matrix(samples$Age)
row.names(data_age) <- samples$sWGS

# Set cancer development
data_dev <- as.matrix(samples$developed)
row.names(data_dev) <- samples$sWGS

## Set Previous cancers
data_history <- as.matrix(samples$previous_cancer)
row.names(data_history) <- samples$sWGS

## Set Cancer Type
data_type <- as.matrix(samples$cancer_type)
row.names(data_type) <- samples$sWGS

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         deletion = "#377EB8", splice_site = "#984EA3", NS = "grey65", pos = "#fb9a99")
col_dev <- c(yes = "#fb9a99", no = "grey95")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", 
              Leukemia = "#add3ac", Glioma = "grey65", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Osteosarcoma = "#ebdc78", 
              Prostate = "#b3d4ff", Sarcoma = "#beb9db", Unknown = "grey95", "NA" = "grey95")

## Set variables
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = "grey95", col = NA))
  # alterations
  n = sum(v)  # how many alterations for current gene in current sample
  h = h
  w = w
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = "white"), just = "top")
}

## Set additional annotations
top_annotation <- HeatmapAnnotation("Developed Cancer" = data_dev,
                                    "Integrated" = anno_lines(data_int,
                                                              add_points = TRUE,
                                                              pch = c(16, NA),
                                                              gp = gpar(col = c("black", "red"),
                                                                        lty = c("solid", "dashed")),
                                                              pt_gp = gpar(col = c("black", NA)),
                                                              size = unit(1, "mm"),
                                                              axis_param = list(side = "left",
                                                                                labels_rot = 90,
                                                                                gp = gpar(fontsize = 5)),
                                                              border = FALSE),
                                    "Breast\nMethylation" = anno_lines(data_breast,
                                                                       add_points = TRUE,
                                                                       pch = c(16, NA),
                                                                       gp = gpar(col = c("black", "red"),
                                                                                 lty = c("solid", "dashed")),
                                                                       pt_gp = gpar(col = c("black", NA)),
                                                                       size = unit(1, "mm"),
                                                                       axis_param = list(side = "left",
                                                                                         labels_rot = 90,
                                                                                         gp = gpar(fontsize = 5)),
                                                                       border = FALSE),
                                    "Pan-Cancer\nMethylation" = anno_lines(data_all,
                                                                           add_points = TRUE,
                                                                           pch = c(16, NA),
                                                                           gp = gpar(col = c("black", "red"),
                                                                                     lty = c("solid", "dashed")),
                                                                           pt_gp = gpar(col = c("black", NA)),
                                                                           size = unit(1, "mm"),
                                                                           axis_param = list(side = "left",
                                                                                             labels_rot = 90,
                                                                                             gp = gpar(fontsize = 5)),
                                                                           border = FALSE),
                                    "Copy Number\n(Tumor Fraction)" = anno_lines(data_ichorCNA,
                                                                                 add_points = TRUE,
                                                                                 pch = c(16, NA),
                                                                                 gp = gpar(col = c("black", "red"),
                                                                                           lty = c("solid", "dashed")),
                                                                                 pt_gp = gpar(col = c("black", NA)),
                                                                                 size = unit(1, "mm"),
                                                                                 axis_param = list(side = "left",
                                                                                                   labels_rot = 90,
                                                                                                   gp = gpar(fontsize = 5)),
                                                                                 border = FALSE),
                                    "Genome Fragment" = anno_lines(data_fragment,
                                                                   add_points = TRUE,
                                                                   pch = c(16, NA),
                                                                   gp = gpar(col = c("black", "red"),
                                                                             lty = c("solid", "dashed")),
                                                                   pt_gp = gpar(col = c("black", NA)),
                                                                   size = unit(1, "mm"),
                                                                   axis_param = list(side = "left",
                                                                                     labels_rot = 90,
                                                                                     gp = gpar(fontsize = 5)),
                                                                   border = FALSE),
                                    "TP53 Fragment" = anno_lines(data_tp53,
                                                                 add_points = TRUE,
                                                                 pch = c(16, NA),
                                                                 gp = gpar(col = c("black", "red"),
                                                                           lty = c("solid", "dashed")),
                                                                 pt_gp = gpar(col = c("black", NA)),
                                                                 size = unit(1, "mm"),
                                                                 axis_param = list(side = "left",
                                                                                   labels_rot = 90,
                                                                                   gp = gpar(fontsize = 5)),
                                                                 border = FALSE),
                                    col = list("Developed Cancer" = col_dev),
                                    gap = unit(2, "mm"),
                                    border = FALSE,
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 7),
                                    show_legend = FALSE,
                                    simple_anno_size = unit(0.3, "cm"),
                                    height = unit(5, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Cancer Status",
                                                  at = c("Positive", "Negative"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "#a6cee3")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Type",
                                                  at = c("Adenocarcinoma", "Adrenal", "Bladder", "Breast", "Endometrial", "Leukemia", 
                                                         "Glioma", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma"),
                                                  legend_gp = gpar(fill = c("#fd7f6f", "#7eb0d5", "#b2e061","#bd7ebe","#ffb55a", "#add3ac",
                                                                            "grey65", "#fdcce5","#8bd3c7","#ebdc78","#b3d4ff","#beb9db")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm"),
                                                  ncol = 1),
                                           Legend(title = "Mutations",
                                                  at = c("Missense", "Frameshift", "Stop", "Deletion", "Not Detected", "Not Sequenced"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#377EB8", "white", "grey65")),
                                                  nrow = 6,
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm")), 
                                           Legend(title = "Analyses",
                                                  type = "lines",
                                                  at = c("Detection\nThreshold"),
                                                  legend_gp = gpar(col = c("red"), 
                                                                   lty = c("dashed")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white"),
                                           Legend(title = "Developed Cancer",
                                                  at = c("Yes", "No"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "grey95")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm"))),
                               direction = "vertical",
                               max_height = unit(3, "in"))

## Set orders
row_order <- row.names(data_mutation)
col_order <- colnames(data_mutation)
col_split <- samples$ext_ID
order <- c("LFS1", "LFS4", "LFS6", "LFS10", "LFS13", "LFS14", "LFS16", "LFS19", "LFS42", "LFS52", "LFS75", "LFS28", "LFS45", "LFS55", "LFS68", "LFS71")
col_split <- factor(col_split, levels = order)

## Generate heatmap
pdf(file.path(outdir, "integration_longitudinal_negs.pdf"), width = 8.5, height = 3)
oncoPrint <- oncoPrint(data_mutation,
                       alter_fun = alter_fun, 
                       col = col,
                       show_heatmap_legend = FALSE,
                       show_pct = FALSE,
                       row_order = row_order,
                       column_order = col_order,
                       top_annotation = top_annotation,
                       right_annotation = NULL,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       show_column_names = FALSE,
                       column_split = col_split,
                       row_gap = unit(0.2, "mm"),
                       row_title_rot = 0,
                       row_names_gp = gpar(fontsize = 7),
                       column_title_side = "top",
                       column_title_rot = 90,
                       column_title_gp = gpar(fontsize = 6),
                       border = FALSE,
                       border_gp = gpar(col = "black"),
                       width = unit(0.15*ncol(data_mutation), "cm"),
                       height = unit(0.3, "cm"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "right")
dev.off()
