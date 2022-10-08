cat("\014")
rm(list = ls())

library(arrow)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(caret)
library(randomForest)
library(limma)

# Compute three different types of DMRs
#BC vs HBC
#BC vs LFS_negative
#BC vs LFS_positive

# ------------------- load data -------------------
setwd('/cluster/projects/pughlab/projects/CHARM/LFS/Ping_medremix')
data <- read_parquet('LFS_0530.parquet')
HBC <- read_parquet('HBC.parquet')
data <- full_join(data, HBC, by=c('bin_chr', 'bin_start', 'bin_end'))

# pbl depleted
row_info <- data[,c('bin_chr', 'bin_start', 'bin_end')]
df_median <- apply(data[,124:137], 1, median)
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

thres_beta <- 0.1
thres_cg <- 1

pbl_del <- df_black_removed[df_black_removed$median < thres_beta,]
pbl_del['median'] <- NULL

dense_cpg <- df_bin_cpg[cpg_density$density >= thres_cg, c('bin_chr', 'bin_start', 'bin_end')]

# combine blacklist removed bins to cpg dense bins
result_bin <- semi_join(pbl_del, dense_cpg, by=c('bin_chr', 'bin_start', 'bin_end'))
Data_mapped <- semi_join(data, result_bin, by=c('bin_chr', 'bin_start', 'bin_end'))# This is the data after preprocessing

row_info <- Data_mapped[,c('bin_chr','bin_start','bin_end')]#This dataframe contains feature names
Data_mapped[,c('bin_chr','bin_start','bin_end')] <- NULL

write_parquet(Data_mapped, 'Data_mapped_0811.parquet')
write_parquet(row_info, 'row_info_0811.parquet')


###########################
###########################
###########################The following section will be run on personal computer
###########################
###########################
nFreq <- 30
nDMR <- 300

setwd('~/Project/cfMeDIP/CHARM/LFS/cfmedip_supervised')
Data_mapped <- read_parquet('Data_mapped_0811.parquet')
row_info <- read_parquet('row_info_0811.parquet')
metadata <- readRDS('~/Project/cfMeDIP/CHARM/LFS/metadata/metadata_0719.rds')

metadata <- metadata[1:123,]

LFS <- Data_mapped[,c(1:123)]
# LFS_negative - 51 negative samples
LFS_neg <- LFS[,metadata$new_status == 'Previvor']
# LFS_positive - non-BC 17 samples
LFS_pos <- LFS[,metadata$new_status == 'Positive' & metadata$cancer_type != 'breast']
# LFS BC 11 samples
LFS_BC <- LFS[,metadata$new_status == 'Positive' & metadata$cancer_type == 'breast']
# HBC 14 samples
HBC <- Data_mapped[,c(124:137)]

# cross-validation
# vs LFS_negative
DMR.classes <- c(rep("One", 11), rep("Others", 11))
Des <- model.matrix(~0 + DMR.classes)
colnames(Des) <- levels(factor(DMR.classes))
Features <- list()

for (k in 1:100) {
  Data_sub <- cbind(LFS_BC, LFS_neg[,sample(c(1:ncol(LFS_neg)), size=11)])
  LimmaFit <- lmFit(Data_sub, Des)%>%
              contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
              eBayes(., trend = FALSE)%>%
              topTable(., number = nrow(Data_sub))

  LimmaFit <- LimmaFit%>%.[order(.$t,decreasing = T),]
  Features[[k]] <- LimmaFit[1:nDMR,]
}

Features_vec <- c()
for (k in 1:100) {
  Features_vec <- c(Features_vec, rownames(Features[[k]]))
}
re <- as.data.frame(table(Features_vec))

re$Features_vec <- as.numeric(as.character(re$Features_vec)) # This is only useful for current script since Features_vec contains index of features
tmp <- re[re$Freq >= nFreq,]

# CV accuracy 0.8714286
BC_neg_DMR <- row_info[tmp$Features_vec,]

# vs LFS_positive
DMR.classes <- c(rep("One", 11), rep("Others", 17))
Des <- model.matrix(~0 + DMR.classes)
colnames(Des) <- levels(factor(DMR.classes))

Data_sub <- cbind(LFS_BC, LFS_pos)
LimmaFit <- lmFit(Data_sub, Des)%>%
  contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
  eBayes(., trend = FALSE)%>%
  topTable(., number = nrow(Data_sub))

LimmaFit <- LimmaFit%>%.[order(.$t,decreasing = T),]
Features_pos <- LimmaFit[1:150,]

# CV accuracy 0.80
BC_pos_DMR <- row_info[rownames(Features_pos),]

# vs HBC
DMR.classes <- c(rep("One", 11), rep("Others", 14))
Des <- model.matrix(~0 + DMR.classes)
colnames(Des) <- levels(factor(DMR.classes))

Data_sub <- cbind(LFS_BC, HBC)
LimmaFit <- lmFit(Data_sub, Des)%>%
  contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
  eBayes(., trend = FALSE)%>%
  topTable(., number = nrow(Data_sub))

LimmaFit <- LimmaFit%>%.[order(.$t,decreasing = T),]
Features_hbc <- LimmaFit[1:100,]

# CV accuracy 0.80
BC_hbc_DMR <- row_info[rownames(Features_hbc),]

BC_DMR_list <- list()
BC_DMR_list[['HBC']] <- BC_hbc_DMR
BC_DMR_list[['pos']] <- BC_pos_DMR
BC_DMR_list[['neg']] <- BC_neg_DMR

saveRDS(BC_DMR_list, 'BC_DMR_list.rds')
BC_DMR_list <- readRDS('BC_DMR_list.rds')

save(Features, Features_pos, Features_hbc, file = "BC_DMRs_list.RData")

# 28 signatures related to BC
test <- semi_join(BC_DMR_list[['neg']], BC_DMR_list[['pos']], by=c('bin_chr', 'bin_start', 'bin_end'))
# 10 signatures
test2 <- semi_join(BC_DMR_list[['neg']], BC_DMR_list[['HBC']], by=c('bin_chr', 'bin_start', 'bin_end'))
# 0 signatures
test3 <- semi_join(BC_DMR_list[['pos']], BC_DMR_list[['HBC']], by=c('bin_chr', 'bin_start', 'bin_end'))

tmp <- full_join(test, test2, by=c('bin_chr', 'bin_start', 'bin_end'))


intersect(tmp$Features_vec,rownames(Features_pos))
# accuracy 0.85 for BC and HBC

write.table(test, 'BC_pos_neg_DMR.tsv', quote=F, row.names=F, sep='\t')
write.table(test2, 'BC_pos_neg_HBC_DMR.tsv', quote=F, row.names=F, sep='\t')
write.table(tmp, 'BC_tmp.tsv', quote=F, row.names=F, sep='\t')

########################### CV
mtx <- cbind(LFS_BC, HBC)[test3,]
mtx <- as.data.frame(t(mtx))

mtx$Type <- c(rep("One", 11), rep("Others", 14))
mtx[is.na(mtx)] <- 0

control <- trainControl(method='repeatedcv', number=10, repeats=1)
tunegrid <- expand.grid(.mtry=c(10))
rf_default <- train(Type~.,
                    data=mtx,
                    method='rf',
                    metric='Accuracy',
                    tuneGrid=tunegrid,
                    trControl=control)
#cat('beta:',thres_beta,',','cg:',thres_cg,'\n')
print(rf_default)
