library(tidyverse)
library(dplyr)
library(vroom)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/pipeline_output"
OICR_path <- "/Users/derekwong/Desktop/H4H/external_data/TGL49_CHARM/LFS/LFS_TS/mafs"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/pipeline_analysis"
project <- "CHARM_LFS"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Get Files
germline <- list.files(path = file.path(path, "HaplotypeCaller/CPSR"), pattern = "*.tsv", full.names = TRUE)
germline_mafs <- list.files(path = file.path(path, "HaplotypeCaller/cohort/VCF2MAF"), pattern = "*.maf", full.names = TRUE)
md5s <- germline_mafs[grep("md5", germline_mafs)]
germline_mafs <- germline_mafs[!(germline_mafs %in% md5s)]

mutect2 <- list.files(path = file.path(path, "MuTect2"), pattern = "*.maf", recursive = TRUE, full.names = TRUE)
md5s <- mutect2[grep("md5", mutect2)]
mutect2 <- mutect2[!(mutect2 %in% md5s)]

OICR <- list.files(path = file.path(OICR_path), pattern = "*.maf", full.names = TRUE)

### Import data
germline_file <- read_tsv(germline)

datalist <- lapply(germline_mafs, function(x){vroom(file = x, skip = 1)})
germline_mafs_file <- do.call("rbind", datalist)

datalist <- lapply(mutect2, function(x){vroom(file = x, skip = 1)})
mutect2_file <- do.call("rbind", datalist)

datalist <- lapply(OICR, function(x){vroom(file = x, skip = 1)})
OICR_file <- do.call("rbind", datalist)

### Parse germline (Happlotype Caller)
germline <- germline_file[, colnames(germline_file) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                           "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                           "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG", "flanking_bps")]

cohort_snps <- germline[germline$Hugo_Symbol %in% c("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "BRCA1", "BRCA2", "PALB2", "APC"), ]
cohort_snps$vaf <- cohort_snps$t_alt_count/cohort_snps$t_depth*100
cohort_snps <- cohort_snps %>% relocate(vaf, .after = t_alt_count)
cohort_snps <- cohort_snps[cohort_snps$vaf > 40, ]
cohort_snps <- cohort_snps[cohort_snps$vaf < 60, ]

germline <- germline[germline$Hugo_Symbol == "TP53", ]
germline <- germline[order(germline$Tumor_Sample_Barcode), ]
germline$vaf <- germline$t_alt_count/germline$t_depth*100
germline <- germline %>% relocate(vaf, .after = t_alt_count)
germline <- germline[germline$t_alt_count > 0, ]
write.table(germline, file.path(outdir, paste0(project, "_haplotypecaller.txt")), sep = "\t", row.names = FALSE)

### Parse germline (Raw mafs)
germline_mafs <- germline_mafs_file[, colnames(germline_mafs_file) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                                          "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                                          "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG", "flanking_bps")]
cohort_snps2 <- germline_mafs[germline_mafs$Hugo_Symbol %in% c("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "BRCA1", "BRCA2", "PALB2", "APC"), ]
cohort_snps2$vaf <- cohort_snps2$t_alt_count/cohort_snps2$t_depth*100
cohort_snps2 <- cohort_snps2 %>% relocate(vaf, .after = t_alt_count)
cohort_snps2 <- cohort_snps2[cohort_snps2$vaf > 40, ]
cohort_snps2 <- cohort_snps2[cohort_snps2$vaf < 60, ]

germline_mafs <- germline_mafs[germline_mafs$Hugo_Symbol == "TP53", ]
germline_mafs <- germline_mafs[!grepl("benign", germline_mafs$PolyPhen), ]
germline_mafs <- germline_mafs[!grepl("benign", germline_mafs$CLIN_SIG), ]
germline_mafs <- germline_mafs[order(germline_mafs$Tumor_Sample_Barcode), ]
germline_mafs$vaf <- germline_mafs$t_alt_count/germline_mafs$t_depth*100
germline_mafs <- germline_mafs %>% relocate(vaf, .after = t_alt_count)
germline_mafs <- germline_mafs[germline_mafs$vaf > 20, ]
germline_mafs <- germline_mafs[germline_mafs$vaf < 80, ]
write.table(germline_mafs, file.path(outdir, paste0(project, "_haplotypecaller_raw.txt")), sep = "\t", row.names = FALSE)

##### Make list of germline snps in cohort
cohort_snps <- rbind(cohort_snps, cohort_snps2)
cohort_snps <- distinct(cohort_snps, Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, .keep_all = TRUE)

### Parse mutations (Mutect2)
mutect2 <- mutect2_file[mutect2_file$Hugo_Symbol %in% c("TP53", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "BRCA1", "BRCA2", "PALB2", "APC"), ]
mutect2 <- mutect2[!(mutect2$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "Intron")), ]
mutect2 <- mutect2[!(mutect2$IMPACT == "LOW"), ]
mutect2 <- mutect2[!grepl("benign", mutect2$PolyPhen), ]
mutect2 <- mutect2[, colnames(mutect2) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                              "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                              "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG", "flanking_bps")]
mutect2$vaf <- mutect2$t_alt_count/mutect2$t_depth*100
mutect2 <- mutect2 %>% relocate(vaf, .after = t_alt_count)
mutect2 <- anti_join(mutect2, cohort_snps, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "HGVSc", "HGVSp_Short"))

recurrent <- mutect2[mutect2$Tumor_Sample_Barcode %like% "_all", ]
recurrent <- recurrent %>%
  group_by(Hugo_Symbol, Variant_Classification, HGVSc, HGVSp_Short, SIFT, PolyPhen, CLIN_SIG) %>%
  dplyr::summarise(N=n())
recurrent <- recurrent[!(recurrent$Hugo_Symbol == "TP53" |
                           recurrent$N < 2), ]
recurrent <- recurrent[!(recurrent$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation")), ]
mutect2 <- anti_join(mutect2, recurrent, by = c("Hugo_Symbol", "HGVSc", "HGVSp_Short"))

mutect2_tp53 <- mutect2[mutect2$Hugo_Symbol == "TP53", ]
mutect2_tp53$patient <- str_remove(mutect2_tp53$Tumor_Sample_Barcode, "_all")
mutect2_tp53$patient <- str_remove(mutect2_tp53$patient, "_dcs")
mutect2_tp53$patient <- str_remove(mutect2_tp53$patient, "_sscs")
mutect2_tp53 <- anti_join(mutect2_tp53, germline, by = c("Hugo_Symbol", "HGVSc", "HGVSp_Short", "patient" = "Tumor_Sample_Barcode"))

mutect2 <- mutect2[!(mutect2$Hugo_Symbol == "TP53"), ]
mutect2$patient <- str_remove(mutect2$Tumor_Sample_Barcode, "_all")
mutect2$patient <- str_remove(mutect2$patient, "_dcs")
mutect2$patient <- str_remove(mutect2$patient, "_sscs")
mutect2 <- mutect2 %>% 
  add_count(Hugo_Symbol, patient, HGVSc) %>% 
  dplyr::filter(n > 1)

mutect2 <- bind_rows(mutect2, mutect2_tp53)
write.table(mutect2, file.path(outdir, paste0(project, "_mutect2.txt")), sep = "\t", row.names = FALSE)

### Parse mutations (OICR)
OICR <- OICR_file[OICR_file$Hugo_Symbol %in% c("TP53", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "BRCA1", "BRCA2", "PALB2", "APC"), ]
OICR <- OICR[!(OICR$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "Intron")), ]
OICR <- OICR[!(OICR$IMPACT == "LOW"), ]
OICR <- OICR[!grepl("benign", OICR$PolyPhen), ]
OICR <- OICR[, colnames(OICR) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                              "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                              "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG", "flanking_bps")]
OICR$vaf <- OICR$t_alt_count/OICR$t_depth*100
OICR <- OICR %>% relocate(vaf, .after = t_alt_count)
OICR <- anti_join(OICR, cohort_snps, by = c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "HGVSc", "HGVSp_Short"))
OICR <- OICR[!(OICR$Tumor_Sample_Barcode %like% "_Ly_"), ]

recurrent <- OICR %>%
  group_by(Hugo_Symbol, Variant_Classification, HGVSc, HGVSp_Short, SIFT, PolyPhen, CLIN_SIG) %>%
  dplyr::summarise(N=n())
recurrent <- recurrent[!(recurrent$Hugo_Symbol == "TP53" |
                           recurrent$N < 2), ]
recurrent <- recurrent[!(recurrent$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation")), ]
OICR <- anti_join(OICR, recurrent, by = c("Hugo_Symbol", "HGVSc", "HGVSp_Short"))

write.table(OICR, file.path(outdir, paste0(project, "_OICR.txt")), sep = "\t", row.names = FALSE)

##### Make a combined mutations list
germline <- rbind(germline, germline_mafs)
germline <- germline[!duplicated(germline[ , 1:9]),]

mutect2 <- mutect2[!duplicated(mutect2[ , 1:9]),]

OICR <- OICR[!duplicated(OICR[ , 1:9]),]
mutations <- bind_rows(germline, mutect2, OICR)
mutations <- mutations[!duplicated(mutations[, 1:8]), ]
mutations <- mutations[ , !(colnames(mutations) %in% c("Tumor_Sample_Barcode", "t_depth", "t_ref", "t_alt", "patient", "n"))]

write.table(mutations, file.path(outdir, paste0(project, "_mutations_list.txt")), sep = "\t", row.names = FALSE)
