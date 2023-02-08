library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(lemon)
library(ggh4x)
library(ggplotify)

### Set paths
path <- ""
path2 <- ""
outdir <- ""

### Find files
ichorCNA <- read.delim(list.files(path, "summary_reviewed.txt", recursive = TRUE, full.names = TRUE))
mutation <- read.delim(file.path(path2, "Oncoplot_full.txt"))
tp53 <- read.delim(file.path(path, "CHARM_LFS_panel_score.txt"))
fragment <- read.delim(file.path(path, "CHARM_LFS_genome_score.txt"))
breast <- read.delim(file.path(path2, "breast_matrix.tsv"))
all <- read.delim(file.path(path2, "pancancer_score_Vrba.tsv"))
samples <- read.delim("sample_list.txt")

### Calculate cumulative breast score
breast <- data.frame(One = colSums(breast),
                     cfMeDIP = colnames(breast))

### Calculate methylation limits
breast_limit <- quantile(breast$One[!(breast$cfMeDIP %in% samples$cfMeDIP)], 0.99)
all_limit <- quantile(all$score[!(all$cfMeDIP %in% samples$cfMeDIP)], 0.99)
tp53_limit <- quantile(tp53$score[tp53$sample %in% samples$TS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)
fragment_limit <- quantile(fragment$score[fragment$sample %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]], 0.9)

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
samples <- samples[!(is.na(samples$sWGS)), ]
samples <- samples[rowSums(is.na(samples[, c(26:31)])) != 6, ]
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
samples$TP53_somatic <- ifelse(samples$TP53_somatic == "", "Not Detected",
                               ifelse(samples$TP53_somatic == "NS", "Not Sequenced", "Somatic Mutation"))
samples <- samples[order(samples$ext_ID,
                         samples$timepoint),]

### Samples to plot
cases <- c("LFS3", "LFS5", "LFS15", "LFS78")
samples <- samples[samples$ext_ID %in% cases, ]

### Convert dates to days
samples$date_blood <- as.Date(samples$date_blood)
data_date <- samples$date_blood[samples$timepoint == "0"]
counts <- samples %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())
data_date <- rep(data_date, counts$N)
data_date <- as.Date(data_date)

samples$days <- difftime(samples$date_blood, data_date, units = "days")
samples$days <- as.numeric(gsub(" days", "", samples$days))/(365/12)

### Format for plotting
samples <- samples[, c("ext_ID", "days", "genome", "ichor", "tp53", "TP53_somatic", "breast", "pancancer", "cancer_status")]
colnames(samples) <- c("ext_ID", "days", "Genome-wide\nFragmentation Score", "IchorCNA\n(Tumor Fraction)","TP53 Gene-level\nFragmentation Score", 
                       "TP53 Mutation", "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score", "cancer_status")
samples$Clinical <- NA
data_plot <- reshape2::melt(samples, id = c("ext_ID", "days", "cancer_status", "TP53 Mutation"))
data_plot <- data_plot[complete.cases(data_plot),]
data_plot$variable <- factor(data_plot$variable, levels = c("Clinical", "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score", "IchorCNA\n(Tumor Fraction)", 
                                                            "Genome-wide\nFragmentation Score", "TP53 Gene-level\nFragmentation Score"))

### Format datapoint colors
data_plot$color <- "Negative"
data_plot$color <- ifelse(data_plot$variable == "Breast Cancer\nMethylation Score" & data_plot$value > breast_limit, "Positive", data_plot$color)
data_plot$color <- ifelse(data_plot$variable == "Pan-Cancer\nMethylation Score" & data_plot$value > all_limit, "Positive", data_plot$color)
data_plot$color <- ifelse(data_plot$variable == "IchorCNA\n(Tumor Fraction)" & data_plot$value > 0.03, "Positive", data_plot$color)
data_plot$color <- ifelse(data_plot$variable == "Genome-wide\nFragmentation Score" & data_plot$value > fragment_limit, "Positive", data_plot$color)
data_plot$color <- ifelse(data_plot$variable == "TP53 Gene-level\nFragmentation Score" & data_plot$value > tp53_limit, "Positive", data_plot$color)

### Make clinical table
data_clinical <- data.frame(ext_ID = c("LFS15", "LFS3", "LFS3", "LFS5", "LFS5", "LFS5", "LFS78"),
                            days = c(39.3, 37.9, 85.97, 0, 18.19, 57.5, 11.8),
                            variable = c(rep("Clinical", 7)),
                            value = c(rep(0, 7)),
                            label = c("Leukemia", "Astrocytoma", "Osteosarcoma", "Glioma", "Melanoma", "Osteosarcoma", "Metastatic\nBreast Cancer"))
data_clinical$variable <- factor(data_clinical$variable, levels = c("Clinical", "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score", "IchorCNA\n(Tumor Fraction)", 
                                                                    "Genome-wide\nFragmentation Score", "TP53 Gene-level\nFragmentation Score"))

### Make imaging table
data_imaging <- data.frame(ext_ID = c("LFS3", "LFS3", "LFS3", "LFS5", "LFS15", "LFS15", "LFS78", "LFS78"),
                           variable = rep("Clinical", 8),
                           value = c(36.1, 75.9, 85.9, 0, 27.3, 39.3, 3, 11.8))
data_imaging$variable <- factor(data_imaging$variable, levels = c( "Clinical", "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score", "IchorCNA\n(Tumor Fraction)", 
                                                                  "Genome-wide\nFragmentation Score", "TP53 Gene-level\nFragmentation Score"))

### Make limits table
limits <- data.frame(variable = c("Clinical", "Genome-wide\nFragmentation Score", "IchorCNA\n(Tumor Fraction)", "TP53 Gene-level\nFragmentation Score", 
                                  "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score"),
                     limit = c(0, fragment_limit, 0.03, tp53_limit, breast_limit, all_limit),
                     type = c("solid", "dashed", "dashed", "dashed", "dashed", "dashed"))
limits$variable <- factor(limits$variable, levels = c("Clinical", "Breast Cancer\nMethylation Score", "Pan-Cancer\nMethylation Score", "IchorCNA\n(Tumor Fraction)", 
                                                      "Genome-wide\nFragmentation Score", "TP53 Gene-level\nFragmentation Score"))

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               panel.spacing = unit(1, "lines"),
               plot.background = element_blank(),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.position = "none",
               legend.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 10),
               strip.text.x = element_text(size = 12),
               strip.text.y.left = element_text(angle = 0),
               strip.placement = "outside",
               axis.text = element_text(size = 10),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

### Make loop to plot each patient individually
for (case in cases) {
  data_case <- data_plot[data_plot$ext_ID == case, ]
  clinical_case <- data_clinical[data_clinical$ext_ID == case, ]
  imaging_case <- data_imaging[data_imaging$ext_ID == case, ]
  
  plot <- ggplot(data_case, aes(days, value)) +
    geom_line(alpha = 1, color = "black", size = 0.5) +
    geom_point(aes(color = color, shape = `TP53 Mutation`), size = 2) +
    geom_point(data = clinical_case, shape = 19, size = 5, color = "red", fill = "red") +
    geom_label(data = clinical_case, aes(label = label), label.size = NA, fill = NA, size = 3.5) +
    geom_vline(data = clinical_case, aes(xintercept = days, linetype = "dashed"), color = "red", size = 0.5, alpha = 0.5) +
    geom_hline(data = limits, aes(yintercept = limit, linetype = type, size = type, alpha = type)) +
    geom_vline(data = imaging_case, aes(xintercept = value, linetype = "solid"), color = "black", size = 1, alpha = 0.25) +
    scale_color_manual(name = "Molecular Signal", values = c(Negative = "black", Positive = "red")) +
    scale_shape_manual(values = c(`Not Sequenced` = 1, `Not Detected` = 19, `Somatic Mutation` = 4)) +
    scale_linetype_manual(values = c(dashed = "dashed", solid = "solid"), guide = "none") +
    scale_size_manual(values = c(dashed = 0.5, solid = 5), guide = 'none') +
    scale_alpha_manual(values = c(dashed = 0.5, solid = 0.25), guide = 'none') +
    guides(linetype = guide_legend(title = "Line Type")) +
    xlab("Months from Baseline") +
    ylab("") +
    labs(color = "Cancer Status") +
    facet_grid2(variable~ext_ID, scales = "free", space = "free_x", switch = "y", axes = "y") +
    facetted_pos_scales(y = list(scale_y_continuous(limits = c(-1, 1)),
                                 scale_y_continuous(limits = c(0, 15)),
                                 scale_y_continuous(limits = c(15, 300)),
                                 scale_y_continuous(limits = c(-0.1, 0.55)),
                                 scale_y_continuous(limits = c(-0.95, -0.15)),
                                 scale_y_continuous(limits = c(-0.75, -0.25)))) +
    theme +
    coord_cartesian(expand = TRUE) +
    scale_x_continuous(position = "top")
  plot
  assign(paste0(case, "_plot"), plot)
}

### Merge plots
plot <- grid.arrange(arrangeGrob(LFS3_plot + theme(plot.margin=unit(c(0.2,-0.4,0.4,-0.4), "cm")), 
                                 LFS15_plot + theme(strip.text.y = element_text(size = 0),
                                                    axis.text.y = element_blank(),
                                                    plot.margin=unit(c(0.2,-0.4,0.4,-0.4), "cm")),
                                 LFS5_plot + theme(strip.text.y = element_text(size = 0),
                                                   axis.text.y = element_blank(),
                                                   plot.margin=unit(c(0.2,-0.4,0.4,-0.4), "cm")),
                                 LFS78_plot + theme(strip.text.y = element_text(size = 0),
                                                    axis.text.y = element_blank(),
                                                    plot.margin=unit(c(0.2,0,0.4,-0.4), "cm"),
                                                    legend.position = "right"), 
                                 widths = c(1.2, 0.55, 0.75, 0.9),
                                 nrow = 1))
plot
ggsave(file.path(outdir, "integration_select_cases.pdf"), plot, width = 12, height = 8, units = "in")




