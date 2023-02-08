# This script works on H4H, find the best hyperparameters
library(arrow)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(caret)
library(randomForest)
library(limma)

setwd('')

data <- read_parquet('neg_HBC.parquet')

row_info <- data[,c('bin_chr', 'bin_start', 'bin_end')]
df_median <- apply(data[,55:68], 1, median)
df_median <- cbind(row_info, df_median)
colnames(df_median) <- c('bin_chr','bin_start','bin_end','median')

# remove EnCODE blacklist
blacklist <- read.csv('~/Project/cfMeDIP/CHARM/LFS/hg38-blacklist.v2.bed.gz',sep='\t', header=F)
colnames(blacklist) <- c('bin_chr', 'bin_start', 'bin_end','issue')
blacklist_gr <- with(blacklist, GRanges(seqnames = Rle(bin_chr) , ranges = IRanges(as.numeric(bin_start),as.numeric(bin_end))))
bin_gr <- with(row_info, GRanges(seqnames = Rle(bin_chr) , ranges = IRanges(as.numeric(bin_start),as.numeric(bin_end))))

Ov <- findOverlaps(query = blacklist_gr, subject = bin_gr)
df_black_removed <- df_median[-Ov@to,]

# compute cpg density
df_bin_cpg <- read.delim('~/Project/cfMeDIP/CHARM/LFS/bin_mapped_cpg.tsv')

cpg_density <- c()
for (i in 1:nrow(df_bin_cpg)) {
  cpg_vec <- strsplit(df_bin_cpg$cpg_list[[i]],',')[[1]]
  cpg_density <- c(cpg_density, length(cpg_vec))
}
cpg_density <- as.data.frame(cpg_density)
colnames(cpg_density) <- "density"

thres_beta <- 0.2
thres_cg <- 1
nFreq <- 50
nDMR <- 150

# thres_beta_vec <- c(0.1)
# thres_cg_vec <- c(1)
# nDMR_vec <- c(300)
# nFreq_vec <- c(30)

pbl_del <- df_black_removed[df_black_removed$median < thres_beta,]
pbl_del['median'] <- NULL

dense_cpg <- df_bin_cpg[cpg_density$density >= thres_cg, c('bin_chr', 'bin_start', 'bin_end')]

# combine blacklist removed bins to cpg dense bins
result_bin <- semi_join(pbl_del, dense_cpg, by=c('bin_chr', 'bin_start', 'bin_end'))

Data_mapped <- semi_join(data, result_bin, by=c('bin_chr', 'bin_start', 'bin_end'))

row_info <- Data_mapped[,c('bin_chr','bin_start','bin_end')]
setwd('')
write_parquet(Data_mapped, 'neg_HBC_filtered_0721.parquet')

# cross-validation
DMR.classes <- c(rep("One", 45), rep("Others", 12))

Data_mapped[,c('bin_chr','bin_start','bin_end')] <- NULL

data_cancer <- Data_mapped[, c(1:51)]
data_control <- Data_mapped[, c(52:65)]

Des <- model.matrix(~0 + DMR.classes)
colnames(Des) <- levels(factor(DMR.classes))
Features <- list()

for (k in 1:100) {
  Data_sub <- cbind(data_cancer[,sample(c(1:ncol(data_cancer)), size=45)],
                data_control[,sample(c(1:ncol(data_control)), size=12)])
  LimmaFit <- lmFit(Data_sub, Des)%>%
              contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
              eBayes(., trend = FALSE)%>%
              topTable(., number = nrow(Data_sub))

  LimmaFit <- LimmaFit%>%.[order(.$t,decreasing = T),]

  Features[[k]] <- LimmaFit[1:nDMR,] # This line we used in saving the results
}


Features_vec <- c()
for (k in 1:100) {
  Features_vec <- c(Features_vec, rownames(Features[[k]]))
}
re <- as.data.frame(table(Features_vec))
re$Features_vec <- as.numeric(as.character(re$Features_vec))
tmp <- re[re$Freq >= nFreq,]

Feature_list <- row_info[tmp$Features_vec,]
setwd('')
write.table(Feature_list, 'neg_HBC_DMRs.tsv', quote=F, row.names = F, sep='\t')
