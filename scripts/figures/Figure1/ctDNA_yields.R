library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)
library(data.table)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/ctDNA_yields"

data_yield <- file.path(path, "ctDNA_yields/ctDNA_yields.txt")
data_healthy <- file.path(path, "ctDNA_yields/hbc_yields.txt")
data_samples <- file.path(path, "samples/sample_list.txt")

### Read in data
data_yield <- read.delim(data_yield)
data_healthy <- read.delim(data_healthy)
data_samples <- read.delim(data_samples)
source("/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Sum samples with multiple replicates
data_yield$Volume[is.na(data_yield$Volume)] <- 0
data_yield$Volume <- as.numeric(data_yield$Volume)
data_yield <- data_yield %>% group_by(LIBERATEID, Time.Point) %>%
  dplyr::summarise(LIB_ID = unique(LIBERATEID),
                   timepoint = unique(Time.Point),
                   volume = sum(Volume),
                   yield = sum(`Total.DNA.Amount..ng.`, na.rm = TRUE))

data_healthy <- data_healthy %>% group_by(Patient.ID, Time.Point) %>%
  dplyr::summarise(LIB_ID = unique(Patient.ID),
                   timepoint = unique(Time.Point),
                   volume = sum(Volume, na.rm = TRUE),
                   yield = sum(`Total.DNA.Amount..ng.`, na.rm = TRUE))

### Format Tables
data_yield <- data_yield[, 3:6]
data_yield$timepoint <- gsub("T", "", data_yield$timepoint)
data_yield$concentration <- data_yield$yield/data_yield$volume
data_yield$concentration[is.infinite(data_yield$concentration)] <- NA
colnames(data_yield) <- c("patient", "timepoint", "volume", "yield", "concentration")
data_yield$type <- "LFS"

data_healthy <- data_healthy[, 3:6]
data_healthy$timepoint <- gsub("T", "", data_healthy$timepoint)
data_healthy$concentration <- data_healthy$yield/data_healthy$volume
colnames(data_healthy) <- c("patient", "timepoint", "volume", "yield", "concentration")
data_healthy$type <- "HBC"

### Merge clinical information with yields
samples_ped <- data_samples[data_samples$LIB_ID %like% "LFS", ]
samples_adult <- data_samples[!(data_samples$LIB_ID %in% samples_ped$LIB_ID), ]

data_ped <- data_yield[, c("patient", "volume", "yield", "concentration", "type")]
data_ped <- data_ped[data_ped$patient %like% "LFS", ]
data_ped <- merge(data_ped, samples_ped, by.x = "patient", by.y = "Patient_ID", all = TRUE)
data_ped <- data_ped[, c("patient", "volume", "yield", "concentration", "type.x", "Age", "years", "Sex", "cancer_status", "previous_cancer")]
data_ped$type.x <- "LFS"

data_adult <- data_yield[!(data_yield$patient %in% data_ped$patient), ]
data_adult <- merge(data_adult, samples_adult, by.x = c("patient", "timepoint"), by.y = c("LIB_ID", "timepoint"), all = TRUE)
data_adult <- data_adult[, c("patient", "volume", "yield", "concentration", "type.x", "Age", "years", "Sex", "cancer_status", "previous_cancer")]
data_adult$type.x <- "LFS"

data_healthy$type.x <- "HBC"
data_healthy$Age <- "HBC"
data_healthy$Sex <- "HBC"
data_healthy$cancer_status <- "HBC"
data_healthy$previous_cancer <- "HBC"
data_healthy$number_previous <- "HBC"
data_healthy <- data_healthy[, c("patient", "volume", "yield", "concentration", "type.x", "Age", "Sex", "cancer_status", "previous_cancer")]

data <- bind_rows(data_ped, data_adult, data_healthy)
data$years[is.na(data$years)] <- 0
data <- data[complete.cases(data), ]
data <- data[!(data$years == 0 &
                 data$type == "LFS"), ]
names(data)[names(data) == "type.x"] <- "type"
data$type <- factor(data$type, levels = c("HBC", "LFS"))
data$Age <- factor(data$Age, levels = c("HBC", "pediatric", "adult"),
                      labels = c("HBC", "Pediatric", "Adult"))
data$Sex <- factor(data$Sex, levels = c("HBC", "male", "female"),
                   labels = c("HBC", "XY", "XX"))
data$cancer_status <- factor(data$cancer_status, levels = c("HBC", "negative", "positive"),
                             labels = c("HBC", "Negative", "Positive"))
data$previous_cancer <- factor(data$previous_cancer, levels = c("HBC", "no", "yes"),
                               labels = c("HBC", "No", "Yes"))

### Remove extreme outliers
Q <- quantile(data$yield, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(data$yield)
up <-  Q[2] + 5*iqr 
data <- data[data$yield < up, ]

Q <- quantile(data$concentration, probs=c(.05, .95), na.rm = FALSE)
iqr <- IQR(data$yield)
up <-  Q[2] + 3*iqr 
data <- data[data$concentration < up, ]

### Calculate p-values for different variables
#Total yield
stats_yield <- data[data$type == "LFS", ] %>%
  group_by(Age)%>% 
  dplyr::summarise(Median=median(yield, na.rm = TRUE),
                   Mean=mean(yield, na.rm = TRUE),
                   SD=sd(yield, na.rm = TRUE),
                   N=n())
a <- t.test(data$yield[data$Age == "Pediatric"],data$yield[data$Age == "Adult"])$p.value
t_test <- c(1, a)
stats_yield$pvalue <- p.adjust(t_test)
stats_yield$annot <- ifelse(stats_yield$pvalue < 0.05 & stats_yield$pvalue > 0.01, "*",
                           ifelse(stats_yield$pvalue < 0.01 & stats_yield$pvalue > 0.001, "**",
                                  ifelse(stats_yield$pvalue < 0.001, "***", "")))

#Type
stats_type <- data %>%
  group_by(type) %>% 
  dplyr::summarise(Median=median(concentration, na.rm = TRUE),
                   Mean=mean(concentration, na.rm = TRUE),
                   SD=sd(concentration, na.rm = TRUE),
                   N=n())
a <- t.test(data$concentration[data$type == "HBC"],data$concentration[data$type == "LFS"])$p.value
t_test <- c(1, a)
stats_type$pvalue <- t_test
stats_type$annot <- ifelse(stats_type$pvalue < 0.05 & stats_type$pvalue > 0.01, "*",
                          ifelse(stats_type$pvalue < 0.01 & stats_type$pvalue > 0.001, "**",
                                 ifelse(stats_type$pvalue < 0.001, "***", "")))

#Age
stats_age <- data[data$type == "LFS", ] %>%
  group_by(Age) %>% 
  dplyr::summarise(Median=median(concentration, na.rm = TRUE),
            Mean=mean(concentration, na.rm = TRUE),
            SD=sd(concentration, na.rm = TRUE),
            N=n())
a <- t.test(data$concentration[data$Age == "Pediatric"],data$concentration[data$Age == "Adult"])$p.value
t_test <- c(1, a)
stats_age$pvalue <- t_test
stats_age$annot <- ifelse(stats_age$pvalue < 0.05 & stats_age$pvalue > 0.01, "*",
                           ifelse(stats_age$pvalue < 0.01 & stats_age$pvalue > 0.001, "**",
                                  ifelse(stats_age$pvalue < 0.001, "***", "")))

# Cancer Status
stats_status <- data[data$type == "LFS", ] %>%
  group_by(cancer_status) %>% 
  dplyr::summarise(Median=median(concentration, na.rm = TRUE),
                   Mean=mean(concentration, na.rm = TRUE),
                   SD=sd(concentration, na.rm = TRUE),
                   N=n())
a <- t.test(data$concentration[data$cancer_status == "Negative"], data$concentration[data$cancer_status == "Positive"])$p.value
t_test <- c(1, a)
stats_status$pvalue <- t_test
stats_status$annot <- ifelse(stats_status$pvalue < 0.05 & stats_status$pvalue > 0.01, "*",
                          ifelse(stats_status$pvalue < 0.01 & stats_status$pvalue > 0.001, "**",
                                 ifelse(stats_status$pvalue < 0.001, "***", "")))

# Sex
stats_sex <- data[data$type == "LFS", ] %>%
  group_by(Sex) %>% 
  dplyr::summarise(Median=median(concentration, na.rm = TRUE),
                   Mean=mean(concentration, na.rm = TRUE),
                   SD=sd(concentration, na.rm = TRUE),
                   N=n())
a <- t.test(data$concentration[data$Sex == "XX"],data$concentration[data$Sex == "XY"])$p.value
t_test <- c(1, a)
stats_sex$pvalue <- t_test
stats_sex$annot <- ifelse(stats_sex$pvalue < 0.05 & stats_sex$pvalue > 0.01, "*",
                          ifelse(stats_sex$pvalue < 0.01 & stats_sex$pvalue > 0.001, "**",
                                 ifelse(stats_sex$pvalue < 0.001, "***", "")))
 
### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
               axis.line = element_line(colour = "black"),
               axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 12))

col <- c("grey65", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99")

### Graph yields
plot_type <- ggplot(data, aes(type, concentration, fill = type)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_type, aes(x = type, y = -2, label = N), size = 4) +
  geom_text(data = stats_type, aes(x = type, y = 35, label = annot), size = 5) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Yield (ng/mL of plasma)") +
  ggtitle("Type") + 
  theme
plot_type

data <- data[data$type == "LFS", ]
plot_age <- ggplot(data, aes(Age, concentration, fill = Age)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_age, aes(x = Age, y = -2, label = N), size = 4) +
  geom_text(data = stats_age, aes(x = Age, y = 35, label = annot), size = 5) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Yield (ng/mL of plasma)") +
  ggtitle("Age") + 
  theme
plot_age

plot_sex <- ggplot(data, aes(Sex, concentration, fill = Sex) ) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_sex, aes(x = Sex, y = -2, label = N), size = 4) +
  geom_text(data = stats_sex, aes(x = Sex, y = 35, label = annot), size = 5) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Yield (ng/mL of plasma)") +
  ggtitle("Sex") + 
  theme
plot_sex

plot_status <- ggplot(data, aes(cancer_status, concentration, fill = cancer_status) ) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_status, aes(x = cancer_status, y = -2, label = N), size = 4) +
  geom_text(data = stats_status, aes(x = cancer_status, y = 35, label = annot), size = 5) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Yield (ng/mL of plasma)") +
  ggtitle("Cancer\nStatus") + 
  theme
plot_status

plot_yield <- ggplot(data[data$type == "LFS", ], aes(Age, yield, fill = Age) ) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_yield, aes(x = Age, y = -20, label = N), size = 4) +
  geom_text(data = stats_yield, aes(x = Age, y = 430, label = annot), size = 5) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Total Yield (ng)") +
  ggtitle("Yield") + 
  theme +
  scale_y_continuous(expand = c(0,0), limits = c(-50, 450))
plot_yield

plot_years <- ggplot(data, aes(years, concentration)) + 
  geom_point(stroke = 0, pch = 21, alpha = 0.5, size = 2, fill = "black") +
  ggtitle("Patient Age") +
  xlab("Age (years)") +
  ylab("Yield (ng/mL of plasma)") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_regline_equation(label.y = 33, aes(label = ..rr.label..)) +
  theme + 
  scale_fill_manual(values = c("grey", "#fa4e0c")) +
  scale_y_continuous(limits = c(0, 35)) + 
  scale_x_continuous(limits = c(0, 75)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
plot_years

Figure <- ggarrange(plot_type, 
                    plot_age + 
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank()), 
                    plot_sex + 
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank()),
                    plot_yield,
                    plot_years,
                    nrow = 1,
                    align = "h",
                    widths = c(1.25, 1, 1, 1.25, 3))
Figure
ggsave(file.path(outdir, "ctDNA_yields.pdf"), Figure, device = "pdf", width = 10, height = 4, units = "in")
