library(tidyverse)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(ComplexHeatmap)

### Set variables
path <- "/Users/derekwong/Documents/GitHub/TGL49_CHARM_LFS/data"
outdir <- "/Users/derekwong/Documents/GitHub/TGL49_CHARM_LFS/scripts/figures/Figure1"

data_samples <- read.delim(file.path(path, "sample_list_github.txt"))
data_samples <- data_samples[!(data_samples$notes == "exclude"), ]

### Plot Ages
label_names <- list("male" = "Male",
                    "female" = "Female")
labeller <- function(variable,value){
  return(label_names[value])
}

plot_age <- ggplot(data_samples[data_samples$Sex %in% c("male", "female"), ], aes(x = years, fill = Sex)) + 
  geom_histogram(position = "identity") + 
  ggtitle("") +
  xlab("Age (years) at Time of Blood Draw") +
  ylab("# of Samples") +
  scale_fill_manual(values = c("#FB9A99","#A6CEE3")) +
  facet_grid(.~Sex, labeller = labeller) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 10),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))
plot_age
ggsave(file.path(outdir, "age_histogram.pdf"), width = 2.5, height = 3)

### Make table of families
data_family <- data_samples %>%
  group_by(family) %>%
  dplyr::summarise(samples=n(),
                   patients=length(unique(ext_ID)))
data_family <- data_family[!(data_family$family == "") &
                             data_family$patients > 1, ]
data_family <- data_family[order(data_family$patients), ]
data_family$family <- factor(data_family$family, levels = data_family$family, 
                             labels = c(1:nrow(data_family)))

### Plot Families
plot_fam <- ggplot(data_family, aes(family, patients)) + 
  geom_bar(stat = "identity", width = 0.75, fill = "grey65") + 
  geom_text(aes(family, patients + 0.25, label = samples), size = 3) +
  ggtitle("") +
  xlab("Family") +
  ylab("# of Patients") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5.5))
plot_fam
ggsave(file.path(outdir, "LFS_family.pdf"), width = 4, height = 3)

### Find Number of patients and demographics
samples <- nrow(data_samples)

excluded_patients <- nrow(data_samples[data_samples$notes == "exclude", ])

data_samples <- data_samples[!(data_samples$notes == "exclude"), ]

unique_patients <- length(unique(data_samples$ext_ID))

ped_patients <- data_samples[data_samples$Age == "pediatric", ]
adult_patients <- data_samples[data_samples$Age == "adult", ]

ped_patients <- ped_patients[!(ped_patients$ext_ID %in% adult_patients$ext_ID), ]
ped_patients <- length(unique(ped_patients$ext_ID))

adult_patients <- length(unique(adult_patients$ext_ID))

UHN_patients <- data_samples[grep("LIB-04", data_samples$sample_parent), ]
UHN_patients <- length(unique(UHN_patients$ext_ID))

SK_patients <- data_samples[-grep("LIB-04", data_samples$sample_parent), ]
SK_patients <- length(unique(SK_patients$ext_ID))

males <- data_samples[data_samples$Sex == "male", ]
male_samples <- nrow(males)
males <- length(unique(males$ext_ID))

females <- data_samples[data_samples$Sex == "female", ]
female_samples <- nrow(females)
females <- length(unique(females$ext_ID))

age_min <- min(data_samples$years, na.rm = TRUE)
age_max <- max(data_samples$years, na.rm = TRUE)

### Find number of patients with serial samples
data_samples <- data_samples[!(data_samples$notes == "exclude"), ]
data_serial <- data_samples[data_samples$ext_ID %in% data_samples$ext_ID[duplicated(data_samples$ext_ID)],]
serial_patients <- length(unique(data_serial$ext_ID))

serial_ped <- data_serial[data_serial$Age == "pediatric", ]
serial_adult <- data_serial[data_serial$Age == "adult", ]

serial_ped <- serial_ped[!(serial_ped$ext_ID %in% serial_adult$ext_ID), ]
serial_ped <- length(unique(serial_ped$ext_ID))

serial_adult <- length(unique(serial_adult$ext_ID))

serial_sums <- data_serial[, 1:7] %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())

serial_upper <- max(serial_sums$N)
serial_median <- median(serial_sums$N)

### Find number of seroconverter patients
data_sero <- distinct(data_serial, ext_ID, cancer_status, .keep_all= TRUE)
data_sero <- data_sero[data_sero$ext_ID %in% data_sero$ext_ID[duplicated(data_sero$ext_ID)],]
data_sero <- data_sero[order(data_sero$ext_ID,
                             data_sero$cancer_status), ]

phenoconverters <- length(unique(data_sero$ext_ID))
data_sero$timepoint <- as.numeric(data_sero$timepoint)

sero_neg_pos <- data_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
sero_neg_pos <- sero_neg_pos[sero_neg_pos$difference > 0, ]
sero_neg_pos <- length(unique(sero_neg_pos$ext_ID))

sero_pos_neg <- data_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
sero_pos_neg <- sero_pos_neg[sero_pos_neg$difference < 0, ]
sero_pos_neg <- length(unique(sero_pos_neg$ext_ID)) - 1

### Find number of positive patients and cancer types
data_pos <- data_samples[data_samples$cancer_status == "positive", ]
positive_samples <- nrow(data_pos)
positive_samples_serial <- nrow(data_pos[data_pos$ext_ID %in% data_pos$ext_ID[duplicated(data_pos$ext_ID)],])
positive_samples_single <- positive_samples - positive_samples_serial
positive_patients <- length(unique(data_pos$ext_ID))

data_pos <- distinct(data_pos, ext_ID, cancer_status, cancer_type, .keep_all= TRUE)
data_pos <- as.data.frame(table(data_pos$cancer_type))
data_pos <- data_pos[order(data_pos$Freq, decreasing = TRUE), ]

### Find number of survivors and previvors
data_surv <- data_samples[!(data_samples$ext_ID %in% data_sero$ext_ID), ]
data_surv <- data_surv[!(data_surv$cancer_status == "positive"), ]

data_neg <- data_samples[data_samples$cancer_status == "negative", ]
negative_samples <- nrow(data_neg)
negative_samples_serial <- nrow(data_neg[data_neg$ext_ID %in% data_neg$ext_ID[duplicated(data_neg$ext_ID)],])
negative_samples_single <- negative_samples - negative_samples_serial
data_neg <- length(unique(data_surv$ext_ID))

data_previvor <- data_surv[data_surv$previous_cancer == "no", ]
data_previvor <- length(unique(data_previvor$ext_ID))

data_survivor <- data_surv[data_surv$previous_cancer == "yes", ]
data_survivor <- length(unique(data_survivor$ext_ID))

data_unknown <- data_neg - data_previvor - data_survivor

### Find number of patients excluded from sequencing
data_ns <- data_samples[data_samples$notes == "NS", ]
ns_total <- nrow(data_ns)

ns_ped <- nrow(data_ns[data_ns$Age == "pediatric", ])
ns_adult <- nrow(data_ns[data_ns$Age == "adult", ])

total_seq <- samples - ns_total
total_patients <- data_samples[!(data_samples$notes == "NS"), ]
total_patients <- length(unique(total_patients$ext_ID))

### Find number of targeted sequencing patients
data_seq <- data_samples[!(data_samples$notes == "NS"), ]

data_ts <- data_seq[!(data_seq$TS == ""), ]
ts_samples <- nrow(data_ts)
ts_pediatric <- nrow(data_ts[data_ts$Age == "pediatric", ])
ts_adult <- nrow(data_ts[data_ts$Age == "adult", ])

ts_patients <- length(unique(data_ts$ext_ID))

ts_excluded <- total_seq - ts_samples
data_excluded <- data_seq[!(row.names(data_seq) %in% row.names(data_ts)), ]
ts_excluded_ped <- nrow(data_excluded[data_excluded$Age == "pediatric", ])
ts_excluded_adult <- nrow(data_excluded[data_excluded$Age == "adult", ])

ts_pos <- nrow(data_ts[data_ts$cancer_status == "positive", ])
ts_neg <- nrow(data_ts[data_ts$cancer_status == "negative", ])

ts_serial <- data_ts[data_ts$ext_ID %in% data_ts$ext_ID[duplicated(data_ts$ext_ID)],]

ts_serial_samples <- nrow(ts_serial)
ts_serial_sums <- ts_serial[, 1:7] %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())

ts_median <- median(ts_serial_sums$N)
ts_upper <- max(ts_serial_sums$N)

ts_serial <- length(unique(ts_serial$ext_ID))

ts_sero <- distinct(data_ts, ext_ID, cancer_status, .keep_all= TRUE)
ts_sero <- ts_sero[ts_sero$ext_ID %in% ts_sero$ext_ID[duplicated(ts_sero$ext_ID)],]
ts_sero <- ts_sero[order(ts_sero$ext_ID,
                         ts_sero$cancer_status), ]

ts_neg_pos <- ts_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
ts_neg_pos <- ts_neg_pos[ts_neg_pos$difference > 0, ]
ts_neg_pos <- length(unique(ts_neg_pos$ext_ID))

ts_pos_neg <- ts_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
ts_pos_neg <- ts_pos_neg[ts_pos_neg$difference < 0, ]
ts_pos_neg <- length(unique(ts_pos_neg$ext_ID))

ts_sero <- length(unique(ts_sero$ext_ID))

### Find number of sWGS samples
data_wg <- data_seq[!(data_seq$sWGS == ""), ]
wg_samples <- nrow(data_wg)
wg_pediatric <- nrow(data_ts[data_ts$Age == "pediatric", ])
wg_adult <- nrow(data_ts[data_ts$Age == "adult", ])

wg_patients <- length(unique(data_wg$ext_ID))

wg_excluded <- total_seq - wg_samples
data_excluded <- data_seq[!(row.names(data_seq) %in% row.names(data_wg)), ]
wg_excluded_ped <- nrow(data_excluded[data_excluded$Age == "pediatric", ])
wg_excluded_adult <- nrow(data_excluded[data_excluded$Age == "adult", ])

wg_pos <- nrow(data_wg[data_wg$cancer_status == "positive", ])
wg_neg <- nrow(data_wg[data_wg$cancer_status == "negative", ])

wg_serial <- data_wg[data_wg$ext_ID %in% data_wg$ext_ID[duplicated(data_wg$ext_ID)],]

wg_serial_samples <- nrow(wg_serial)
wg_serial_sums <- wg_serial[, 1:7] %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())

wg_median <- median(wg_serial_sums$N)
wg_upper <- max(wg_serial_sums$N)

wg_serial <- length(unique(wg_serial$ext_ID))

wg_sero <- distinct(data_wg, ext_ID, cancer_status, .keep_all= TRUE)
wg_sero <- wg_sero[wg_sero$ext_ID %in% wg_sero$ext_ID[duplicated(wg_sero$ext_ID)],]
wg_sero <- wg_sero[order(wg_sero$ext_ID,
                         wg_sero$cancer_status), ]

wg_neg_pos <- wg_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
wg_neg_pos <- wg_neg_pos[wg_neg_pos$difference > 0, ]
wg_neg_pos <- length(unique(wg_neg_pos$ext_ID))

wg_pos_neg <- wg_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
wg_pos_neg <- wg_pos_neg[wg_pos_neg$difference < 0, ]
wg_pos_neg <- length(unique(wg_pos_neg$ext_ID))

wg_sero <- length(unique(wg_sero$ext_ID))

### Find number of cfMeDIP samples
data_cm <- data_seq[!(data_seq$cfMeDIP == ""), ]
cm_samples <- nrow(data_wg)
cm_pediatric <- nrow(data_ts[data_ts$Age == "pediatric", ])
cm_adult <- nrow(data_ts[data_ts$Age == "adult", ])

cm_patients <- length(unique(data_wg$ext_ID))

cm_excluded <- total_seq - cm_samples
data_excluded <- data_seq[!(row.names(data_seq) %in% row.names(data_wg)), ]
cm_excluded_ped <- nrow(data_excluded[data_excluded$Age == "pediatric", ])
cm_excluded_adult <- nrow(data_excluded[data_excluded$Age == "adult", ])

cm_pos <- nrow(data_wg[data_wg$cancer_status == "positive", ])
cm_neg <- nrow(data_wg[data_wg$cancer_status == "negative", ])

cm_serial <- data_wg[data_wg$ext_ID %in% data_wg$ext_ID[duplicated(data_wg$ext_ID)],]

cm_serial_samples <- nrow(cm_serial)
cm_serial_sums <- cm_serial[, 1:7] %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())

cm_median <- median(cm_serial_sums$N)
cm_upper <- max(cm_serial_sums$N)

cm_serial <- length(unique(cm_serial$ext_ID))

cm_sero <- distinct(data_wg, ext_ID, cancer_status, .keep_all= TRUE)
cm_sero <- cm_sero[cm_sero$ext_ID %in% cm_sero$ext_ID[duplicated(cm_sero$ext_ID)],]
cm_sero <- cm_sero[order(cm_sero$ext_ID,
                         cm_sero$cancer_status), ]

cm_neg_pos <- cm_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
cm_neg_pos <- cm_neg_pos[cm_neg_pos$difference > 0, ]
cm_neg_pos <- length(unique(cm_neg_pos$ext_ID))

cm_pos_neg <- cm_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
cm_pos_neg <- cm_pos_neg[cm_pos_neg$difference < 0, ]
cm_pos_neg <- length(unique(cm_pos_neg$ext_ID))

cm_sero <- length(unique(cm_sero$ext_ID))

### Sum up statistics into summaries
# general stats
descriptors <- c("samples", "patients", "pediatric", "adult", "UHN", "SickKids", "male", "male samples", "female", "female samples", "min age", "max age", 
                 "positive samples", "positive patients", "not sequenced", "ns pediatric", "ns adult")
stats <- c(samples, unique_patients, ped_patients, adult_patients, UHN_patients, SK_patients, males, male_samples, females, female_samples, age_min, age_max, 
           positive_samples, positive_patients, ns_total, ns_ped, ns_adult)
stats_cohort <- as.data.frame(cbind(descriptors, stats))

#s eroconverters and serial samples
descriptors <- c("serial patients", "pediatric", "adult", "median", "min", "max", "phenoconverters", "negative to positive", "positive to negative",
                 "negative patients", "previvors", "survivors", "unknown")
stats <- c(serial_patients, serial_ped, serial_adult, serial_median, 2, serial_upper, phenoconverters, sero_neg_pos, sero_pos_neg,
           data_neg, data_previvor, data_survivor, data_unknown)
stats_serial <- as.data.frame(cbind(descriptors, stats))

# Targeted Sequencing
descriptors <- c("samples", "pediatric", "adult", "excluded", "excluded pediatric", "excluded adult", "patients", "positive", "negative", 
                 "serial samples", "serial patients", "median", "min", "max","phenoconverters", "negative to positive", "positive to negative")
stats <- c(ts_samples, ts_pediatric, ts_adult, ts_excluded, ts_excluded_ped, ts_excluded_adult, ts_patients, ts_pos, ts_neg, 
           ts_serial_samples, ts_serial, ts_median, 2, ts_upper, ts_sero, ts_neg_pos, ts_pos_neg)
stats_ts <- as.data.frame(cbind(descriptors, stats))

# sWGS
descriptors <- c("samples", "pediatric", "adult", "excluded", "excluded pediatric", "excluded adult", "patients", "positive", "negative", 
                 "serial samples", "serial patients", "median", "min", "max","phenoconverters", "negative to positive", "positive to negative")
stats <- c(wg_samples, wg_pediatric, wg_adult, wg_excluded, wg_excluded_ped, wg_excluded_adult, wg_patients, wg_pos, wg_neg, 
           wg_serial_samples, wg_serial, wg_median, 2, wg_upper, wg_sero, wg_neg_pos, wg_pos_neg)
stats_wg <- as.data.frame(cbind(descriptors, stats))

# cfMeDIP
descriptors <- c("samples", "pediatric", "adult", "excluded", "excluded pediatric", "excluded adult", "patients", "positive", "negative", 
                 "serial samples", "serial patients", "median", "min", "max","phenoconverters", "negative to positive", "positive to negative")
stats <- c(cm_samples, cm_pediatric, cm_adult, cm_excluded, cm_excluded_ped, cm_excluded_adult, cm_patients, cm_pos, cm_neg, 
           cm_serial_samples, cm_serial, cm_median, 2, cm_upper, cm_sero, cm_neg_pos, cm_pos_neg)
stats_cm <- as.data.frame(cbind(descriptors, stats))

### Make a Sankey Diagram (Patients)
nodes = data.frame("name" = 
                     c(paste0("SickKids Hospital (", SK_patients, ")"),                      #0
                       paste0("Princess Margaret (", UHN_patients, ")"),                     #1
                       paste0("Adult (", adult_patients, ")"),                               #2
                       paste0("Pediatric (", ped_patients, ")"),                             #3
                       paste0("Patients (", unique_patients, ")"),                           #4
                       paste0("Phenoconverters (", phenoconverters, ")"),                    #5
                       paste0("Forward (", sero_neg_pos, ")"),                               #6
                       paste0("Reverse (", sero_pos_neg, ")"),                               #7
                       paste0("Active Cancer (", positive_patients, ")"),                  #8
                       paste0("Cancer Free (", data_neg, ")"),                           #9
                       paste0("LFS Past Cancer (", data_survivor, ")"),                            #10
                       paste0("LFS Healthy (", data_previvor, ")")                             #11
                       ))
nodes$group <- as.factor(c("a","b","c","d","e","f","g","h","i","j","k","l"))
my_color <- 'd3.scaleOrdinal() .domain(["a","b","c","d","e","f","g","h","i","j"]) 
.range(["#A6CEE3","#1F78B4","#FF7F00","#FDBF6F","#B2DF8A","#CAB2D6","#FB9A99","#FDBF6F","#FB9A99","#A6CEE3","#FDBF6F","#B2DF8A"])'

links = as.data.frame(matrix(c(
  0, 3, ped_patients,                  # Each row represents a link. (source, target, value)
  0, 2, adult_patients - UHN_patients, 
  1, 2, UHN_patients,                  
  2, 4, excluded_patients,                             
  2, 4, adult_patients,
  3, 4, ped_patients,
  4, 8, positive_patients,
  4, 9, data_neg,
  8, 5, phenoconverters,
  5, 6, sero_neg_pos,
  5, 7, sero_pos_neg,
  9, 10, data_survivor,
  9, 11, data_previvor
  ),
  byrow = TRUE, ncol = 3))

names(links) = c("source", "target", "value")

sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value",  NodeID = "name",
                        fontSize= 12, nodeWidth = 30, NodeGroup = "group",
                        colourScale = my_color, sinksRight = FALSE,
                        width = 1100, height = 200, iterations = 1)
sankey
saveWidget(sankey, file = file.path(outdir, "sankey.html"))
webshot(file.path(outdir, "sankey.html"), file.path(outdir, "cohort_sankey.pdf"), vwidth = 1100, vheight = 150)

### Make Sankey Diagram (Samples)
nodes = data.frame("name" = 
                     c(paste0("Samples (", samples, ")"),                       #0
                       paste0("Cancer Negative (", negative_samples, ")"),          #1
                       paste0("Cancer Positive (", positive_samples, ")"),        #2
                       paste0("Serial (", positive_samples_serial, ")"),        #3
                       paste0("Single (", positive_samples_single, ")"),        #4
                       paste0("Serial (", negative_samples_serial, ")"),        #5
                       paste0("Single (", negative_samples_single, ")")         #6
                     ))
nodes$group <- as.factor(c("a","b","c","d","e","f","g"))
my_color <- 'd3.scaleOrdinal() .domain(["a","b","c","d","e","f","g"]) 
.range(["#FDBF6F","#A6CEE3","#FB9A99","#B2DF8A","#CAB2D6","#B2DF8A","#CAB2D6"])'

links = as.data.frame(matrix(c(
  0, 1, negative_samples,                  # Each row represents a link. (source, target, value)
  0, 2, positive_samples, 
  1, 5, negative_samples_serial,                  
  1, 6, negative_samples_single,                             
  2, 3, positive_samples_serial,
  2, 4, positive_samples_single
),
byrow = TRUE, ncol = 3))

names(links) = c("source", "target", "value")

sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value",  NodeID = "name",
                        fontSize= 12, nodeWidth = 30, NodeGroup = "group",
                        colourScale = my_color, sinksRight = FALSE,
                        width = 700, height = 200, iterations = 1)
sankey
saveWidget(sankey, file = file.path(outdir, "samples_sankey.html"))
webshot(file.path(outdir, "samples_sankey.html"), file.path(outdir, "samples_sankey.pdf"), vwidth = 700, vheight = 175)


### Make a table for upset plot of sequencing
data_seq <- data_samples[, c("TS", "sWGS", "cfMeDIP")]
data_seq <- ifelse(data_seq == "", 0, 1)
colnames(data_seq) <- c("Targeted", "sWGS", "cfMeDIP")

### Make upset plot
m1 <- make_comb_mat(data_seq)
cs <- comb_size(m1)
upset <- make_comb_mat(data_seq)

UpSet <- UpSet(upset,
      top_annotation = HeatmapAnnotation(
        "Samples" = anno_barplot(cs,
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE, 
                                 gp = gpar(fill = "black"),
                                 height = unit(4, "cm")), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      
      comb_order = order(comb_size(m1), decreasing = TRUE))

pdf(file.path(outdir, "sequencing_upset.pdf"), width = 3, height = 3)
ht <- draw(UpSet)
od = column_order(ht)
decorate_annotation("Samples", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("bottom"), 
            gp = gpar(fontsize = 10, col = "#404040"), rot = 0)
})
dev.off()
