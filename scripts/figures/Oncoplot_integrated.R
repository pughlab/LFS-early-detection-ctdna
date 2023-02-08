library(ComplexHeatmap)
library(dplyr)
library(stringr)

path <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot_integrated"
path_onco <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/oncoplot/Oncoplot_full.txt"
path_ichor <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/ichorCNA/CHARM_LFS_ichorCNA_summary_reviewed.txt"
path_samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
path_panel <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/fragment_score/CHARM_LFS_panel_score.txt"
path_swgs <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/fragment_score/CHARM_LFS_genome_score.txt"

## Import data
onco <- read.delim(path_onco)
ichor <- read.delim(path_ichor)
samples <- read.delim(path_samples)
panel <- read.delim(path_panel)
swgs <- read.delim(path_swgs)

## Remove failed samples and format sample sheet
failed_TS <- c("TGL49_0025_Cf_U_PE_334_TS", "TGL49_0035_Cf_U_PE_370_TS", "TGL49_0041_Cf_U_PE_327_TS", "TGL49_0209_Cf_U_PE_378_TS", "TGL49_0010_Cf_U_PE_334_TS")
failed_WG <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")

onco <- onco[!(onco$sample_ID %in% failed_TS), ]
panel <- panel[panel$sample %in% onco$sample_ID, ]

ichor <- ichor[!(ichor$sample %in% failed_WG), ]
swgs <- swgs[swgs$sample %in% ichor$sample, ]

samples <- samples[samples$TS %in% onco$sample_ID |
                    samples$sWGS %in% ichor$sample, ]

### Merge all data together
onco[, c(5:13)][is.na(onco[, c(5:13)])] <- ""
colnames(panel) <- c("sample", "tp53")
colnames(swgs) <- c("sample", "genome")

data <- merge(ichor[, c("sample", "TF", "TF_short")], samples, by.x = "sample", by.y = "sWGS")
data <- merge(data, onco, by.x = "TS", by.y = "sample_ID", all = TRUE)
data <- merge(data, panel, by.x = "TS", by.y = "sample", all = TRUE)
data <- merge(data, swgs, by = "sample", all = TRUE)

### Calculate scores threshold
panel_limit <- quantile(data$tp53[data$cancer_status == "negative" & data$previous_cancer == "no"], 0.9, na.rm = TRUE)
swgs_limit <- quantile(data$genome[data$cancer_status == "negative" & data$previous_cancer == "no"], 0.9, na.rm = TRUE)

data$tp53 <- ifelse(data$tp53 > panel_limit, "pos", "")
data$tp53[is.na(data$tp53)] <- "NS"
data$genome <- ifelse(data$genome > swgs_limit, "pos", "")

### Format data
data[is.na(data)] <- "NS"
data$vaf <- ifelse(data$vaf == "NS", NA, data$vaf)
data$vaf <- as.numeric(data$vaf)

### Set ichorCNA
data$ichor <- ifelse(data$TF > data$TF_short, data$TF, data$TF_short)

## Set clinical information and sample order
data$cancer_status <- factor(data$cancer_status, levels = c("positive", "negative", "unknown"),
                             labels = c("Cancer Positive", "Cancer Negative", "Unknown"))

data <- data[order(factor(data$cancer_status, levels = c("Cancer Negative", "Cancer Positive")),
                   data$ichor,
                   factor(data$TP53_somatic, levels = c("NS", "", unique(data$TP53_somatic)[!(unique(data$TP53_somatic) %in% c("NS", ""))])),
                   factor(data$tp53, levels = c("NS", "neg", "pos")),
                   factor(data$genome, levels = c("neg", "pos")),
                   data$vaf, decreasing = TRUE), ]

## Set Age
data$Age <- factor(data$Age, levels = c("adult", "pediatric"),
                        labels = c("Adult", "Pediatric"))
data_age <- as.matrix(data$Age)
row.names(data_age) <- data$TS

## Set cancer type
data$cancer_type <- factor(data$cancer_type, levels = c("adrenal", "appendiceal_adenocarcinoma", "astrocytoma", "bladder", "breast", "endometrial", 
                                                        "glioma", "leukemia", "lung", "lymphoma", "osteosarcoma", "prostate", "sarcoma", ""),
                                labels = c("Adrenal", "Adenocarcinoma", "Glioma", "Bladder", "Breast", "Endometrial", "Glioma",
                                           "Leukemia", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma", "NA"))
data_type <- as.matrix(data$cancer_type)
row.names(data_type) <- data$sample

## Set cancer previous
data$previous_cancer <- factor(data$previous_cancer, levels = c("yes", "no", ""),
                                      labels = c("Yes", "No", "NA"))
data_previous <- as.matrix(data$previous_cancer)
row.names(data_previous) <- data$sample

## Set Stage
data$stage <- factor(data$stage, levels = c("low", "high", "unknown", ""),
                             labels = c("Stage 0/I/II", "Stage III/IV", "Unknown", "NA"))
data_stage <- as.matrix(data$stage)
row.names(data_stage) <- data$sample

## Format ichorCNA Tumour Fraction
data_ichor <- data.frame(TF = data$ichor,
                         limit = 0.03,
                         vaf = data$vaf,
                         zero = 0)
row.names(data_ichor) <- data$sample
data_ichor <- as.matrix(data_ichor)

## Format oncoprint heatmap
data_onco <- data[, c("genome", "tp53", "TP53_germline", "TP53_somatic", "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC")]
data_onco[is.na(data_onco)] <- ""
row.names(data_onco) <- data$sample
data_onco <- as.data.frame(t(data_onco))
row.names(data_onco) <- c("Genome Fragmentation Score", "TP53 Fragmentation Score", "TP53 (Germline)", "TP53 (Somatic)", "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC")
data_onco <- as.matrix(data_onco)
data_onco <- data_onco[ , data$sample]

## Make oncoprint barplot (minus GL)
data_bar <- data_onco[-c(1:3), ]
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
data_bar <- rbind(missense, missense2, frameshift, stop, splice_site)
colnames(data_bar) <- colnames(data_onco)
data_bar <- t(data_bar)

### Calculate Alteration Percents
data_percent <- data_onco
data_percent[data_percent == "NS"] <- NA
data_percent <- ifelse(data_percent == "", 0, 1)
data_percent <- rowSums(data_percent, na.rm = TRUE)/c(nrow(onco), nrow(ichor), rep(nrow(onco), 11))*100
data_percent <- round(data_percent, 1)
data_percent <- paste0(data_percent, "%")

## Set colours
col <- c(missense = "#4DAF4A", missense2 = "#4DAF4A", frameshift = "black", stop = "#A65628",
         splice_site = "#984EA3", deletion = "#377EB8", duplication = "#E41A1C", 
         pos = "#fb9a99", NS = "grey65")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", 
              Leukemia = "#add3ac", Glioma = "grey65", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Osteosarcoma = "#ebdc78", 
              Prostate = "#b3d4ff", Sarcoma = "#beb9db", Unknown = "grey95", "NA" = "grey95")
col_previous <- c(Yes = "#fb9a99", No = "grey95", "NA" = "grey95")
col_stage <- c("Stage 0/I/II" = "#ebdc78", "Stage III/IV" = "#fb9a99", "Unknown" = "grey65", "NA" = "grey95")

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
top_annotation <- HeatmapAnnotation("Tumour Fraction" = anno_lines(data_ichor,
                                                                   add_points = TRUE,
                                                                   pch = c(16, NA, 4, NA),
                                                                   gp = gpar(col = c(NA, "red", NA, NA),
                                                                             lty = c("solid", "dashed", "solid", "solid")),
                                                                   pt_gp = gpar(col = c("black", NA, "red", NA)),
                                                                   size = unit(1.25, "mm"),
                                                                   axis_param = list(side = "left",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 8)),
                                                                   border = FALSE,
                                                                   height = unit(2, "cm")),
                                    "# of Somatic\nMutations" = anno_barplot(data_bar, 
                                                                             gp = gpar(fill = col,
                                                                                       col = NA),
                                                                             border = FALSE,
                                                                             height = unit(1.5, "cm"),
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
                                    show_legend = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                                    gap = unit(c(2,1,0,0,0), "mm"))

right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(c("missense", "frameshift", "stop", "splice_site",
                                                                        "deletion", "gain", "inframe_deletion", "inframe_insertion", "translation"),
                                                                      border = FALSE, 
                                                                      width = unit(1.5, "cm"), 
                                                                      axis_param = list(side = "top", labels_rot = 0, gp = gpar(fontsize = 7))),
                                 "Percent" = anno_text(data_percent, 
                                                       gp = gpar(fontsize = 7)))
  

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
                                                  at = c("Adenocarcinoma", "Adrenal", "Bladder", "Breast", "Endometrial", "Leukemia", 
                                                         "Glioma", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma"),
                                                  legend_gp = gpar(fill = c("#fd7f6f", "#7eb0d5", "#b2e061","#bd7ebe","#ffb55a", "#add3ac",
                                                                            "grey65", "#fdcce5","#8bd3c7","#ebdc78","#b3d4ff","#beb9db")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 4),
                                           Legend(title = "IchorCNA",
                                                  graphics = list(
                                                    function(x, y, w, h) grid.points(x, y, w*0.5, h*0.5, gp = gpar(col = "black"), pch = 16),
                                                    function(x, y, w, h) grid.points(x, y, w*0.5, h*0.5, gp = gpar(col = "red"), pch = 4),
                                                    function(x, y, w, h) grid.rect(x, y, w, h*0, gp = gpar(col = "red"))),
                                                  labels = c("Tumor Fraction", "Somatic TP53 VAF", "IchorCNA LOD"),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  background = "white"),
                                           Legend(title = "Fragment Score",
                                                  at = c("Positive"),
                                                  legend_gp = gpar(fill = c("#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")), 
                                           Legend(title = "Alterations",
                                                  at = c("Missense", "Frameshift", "Stop", "Splice Site", "Deletion", "Insufficient/Failed\nSequencing"),
                                                  legend_gp = gpar(fill = c("#4DAF4A", "black", "#A65628", "#984EA3", "#377EB8", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 3)), 
                               max_height = unit(6.5, "in"),
                               row_gap = unit(0.25, "cm"),
                               direction = "horizontal")

## Set splits
column_split <- data$cancer_status
column_split <- factor(column_split, levels = c("Cancer Positive", "Cancer Negative"))
row_order <- row.names(data_onco)
col_order <- colnames(data_onco)

## Generate oncoprint
pdf(file.path(path,"Oncoprint_integrated.pdf"), width = 10, height = 5.5)
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       show_heatmap_legend = FALSE,
                       row_order = row_order,
                       row_names_side = "left",
                       column_order = col_order,
                       show_pct = FALSE,
                       top_annotation = top_annotation,
                       right_annotation = right_annotation,
                       row_names_gp = gpar(fontsize = 10, fontface = c(rep("plain", 2), rep("italic", 12))),
                       column_split = column_split,
                       border = FALSE,
                       height = unit(0.35*nrow(data_onco), "cm"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()
