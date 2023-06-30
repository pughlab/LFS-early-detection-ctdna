library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- ""
outdir <- ""

data_samples <- read.delim(file.path(path, "sample_list.txt"))
data_samples <- data_samples[!(data_samples$sWGS %in% c("TGL49_0041_Cf_U_PE_317_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0209_Cf_U_PE_373_WG")), ] ### these patients excluded due to CHIP/sample swaps
data_samples <- data_samples[!(data_samples$cfMeDIP %in% c("TGL49_0316_Pl_T_PE_301_CM", "TGL49_0060_Cf_n_PE_301_CM")), ] # These samples excluded (unsure origins)

data_samples <- data_samples[, c("ext_ID", "timepoint", "Age", "cancer_status", "cfMeDIP", "sWGS", "TS")]

### Order samples
order <- paste0("LFS", 1:92)
data_samples <- data_samples[order(factor(data_samples$ext_ID, levels = order),
                                   data_samples$timepoint), ]

### Format data
data_samples$cfMeDIP <- ifelse(data_samples$cfMeDIP == "", "B", "A")
data_samples$sWGS <- ifelse(data_samples$sWGS == "", "B", "A")
data_samples$TS <- ifelse(data_samples$TS == "", "B", "A")

### Split samples dataframe
datalist <- list(a = data_samples[1:24, ],
                 b = data_samples[25:48, ],
                 c = data_samples[49:72, ],
                 d = data_samples[73:98, ],
                 e = data_samples[99:123, ],
                 f = data_samples[124:147, ],
                 g = data_samples[148:171, ],
                 h = data_samples[172: nrow(data_samples), ])

names <- c("a","b","c","d","e","f","g","h")

### Make oncoplot for each dataframe
for (i in c(1:8)) {
  data <- datalist[[i]]
  name <- names[[i]]
  
  ### Make sequencing table
  data_seq <- as.matrix(data[, c("cfMeDIP", "sWGS", "TS")])
  row.names(data_seq) <- paste0(data$ext_ID, "_", data$timepoint)
  
  ### Set cancer status
  data_cancer <- as.matrix(data$cancer_status)
  row.names(data_cancer) <- paste0(data$ext_ID, "_", data$timepoint)
  
  ### Set patient age
  data_age <- as.matrix(data$Age)
  row.names(data_age) <- paste0(data$ext_ID, "_", data$timepoint)
  
  ### Set data splits
  row_split <- data$ext_ID
  row_split <- factor(row_split, levels = order)
  
  ### Set orders
  column_order <- colnames(data_seq)
  row_order <- row.names(data_seq)
  
  ## Set colours
  col <- c(A = "#B2DF8A", B = "grey65")
  col_status <- c(negative = "#A6CEE3", positive = "#FB9A99")
  col_age <- c(adult = "#6A3D9A", pediatric = "#CAB2D6")
  
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
  left_annotation <- rowAnnotation("Patient Age" = data_age,
                                   "Cancer Status" = data_cancer,
                                   col = list("Patient Age" = col_age,
                                              "Cancer Status" = col_status),
                                   show_annotation_name = TRUE,
                                   annotation_name_side = "bottom",
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = 7),
                                   show_legend = FALSE,
                                   simple_anno_size = unit(0.25, "cm"))
  
  ## Generate heatmap
  oncoPrint <- oncoPrint(data_seq,
                         alter_fun = alter_fun, 
                         col = col,
                         show_heatmap_legend = FALSE,
                         show_pct = FALSE,
                         row_order = row_order,
                         column_order = column_order,
                         left_annotation = left_annotation,
                         top_annotation = NULL,
                         right_annotation = NULL,
                         show_row_names = FALSE,
                         row_names_side = "left",
                         show_column_names = TRUE,
                         row_split = row_split,
                         row_gap = unit(0.5, "mm"),
                         row_title_rot = 0,
                         row_title_gp = gpar(fontsize = 7),
                         column_names_gp = gpar(fontsize = 7),
                         border = FALSE,
                         border_gp = gpar(col = "black"),
                         width = unit(0.25*ncol(data_seq), "cm"),
                         height = unit(0.25*nrow(data_seq), "cm"))
  assign(name, oncoPrint)
}

### Make the Legend=
annotation_legend = packLegend(list = list(Legend(title = "Patient Age",
                                                  at = c("Adult", "Pediatric"),
                                                  legend_gp = gpar(fill = c("#6A3D9A", "#CAB2D6")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Status",
                                                  at = c("Positive", "Negative"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "#a6cee3")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Sequencing",
                                                  at = c("Sequenced", "Not Sequebnced"),
                                                  legend_gp = gpar(fill = c("#B2DF8A", "grey65")),
                                                  title_gp = gpar(fontsize = 7),
                                                  labels_gp = gpar(fontsize = 6),
                                                  grid_height = unit(1, "mm"),
                                                  ncol = 1)),
                               direction = "horizontal",
                               max_height = unit(3, "in"))

pdf(file.path(outdir, "sample_seq_legend.pdf"), height = 1, width = 3)
draw(annotation_legend)
dev.off()

### Save plots
pdf(file.path(outdir, "sample_seq.pdf"), height = 3.5, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 8)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(a, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(b, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(c, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(d, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
draw(e, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 6))
draw(f, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7))
draw(g, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 8))
draw(h, newpage = FALSE)
upViewport()
dev.off()

