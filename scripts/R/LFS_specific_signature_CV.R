# Extract pure negative and HBC samples
library(arrow)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(caret)
library(randomForest)
library(limma)

# ------------------- load data -------------------
setwd('/cluster/projects/pughlab/projects/CHARM/LFS/Ping_medremix')

data <- read_parquet('neg_HBC.parquet')

# ------------------- identify DMR -------------------
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

thres_beta_vec <- c(0.1,0.2,0.3)
thres_cg_vec <- c(1,2,3,4,5,6)
nDMR_vec <- c(150,300)
nFreq_vec <- c(30,50,70)
# N <- 14

DMR.classes <- c(rep("One", 45), rep("Others", 12))

for (nDMR in nDMR_vec) {
  for (nFreq in nFreq_vec) {
    for (thres_beta in thres_beta_vec) {
      for (thres_cg in thres_cg_vec) {
        pbl_del <- df_black_removed[df_black_removed$median < thres_beta,]
        pbl_del['median'] <- NULL

        dense_cpg <- df_bin_cpg[cpg_density$density >= thres_cg, c('bin_chr', 'bin_start', 'bin_end')]

        # combine blacklist removed bins to cpg dense bins
        result_bin <- semi_join(pbl_del, dense_cpg, by=c('bin_chr', 'bin_start', 'bin_end'))

        Data_mapped <- semi_join(data, result_bin, by=c('bin_chr', 'bin_start', 'bin_end'))

        # cross-validation
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
          Features[[k]] <- LimmaFit[1:nDMR,]
        }

        Features_vec <- c()
        for (k in 1:100) {
          Features_vec <- c(Features_vec, rownames(Features[[k]]))
        }
        re <- as.data.frame(table(Features_vec))

        # serious error here
        # re$Features_vec <- as.numeric(as.character(re$Features_vec))
        tmp <- re[re$Freq >= nFreq,]
        # CV
        mtx <- Data_mapped[tmp$Features_vec,]
        mtx <- as.data.frame(t(mtx))

        mtx$Type <- c(rep("One", 51), rep("Others", 14))
        mtx[is.na(mtx)] <- 0

        control <- trainControl(method='repeatedcv', number=10, repeats=1)
        tunegrid <- expand.grid(.mtry=c(10))
        rf_default <- train(Type~.,
                            data=mtx,
                            method='rf',
                            metric='Accuracy',
                            tuneGrid=tunegrid,
                            trControl=control)
        cat('beta:',thres_beta,',','cg:',thres_cg,'\n')
        cat('nDMR:',nDMR,',','nFreq:',nFreq,'\n')
        print(rf_default)
      }
    }
  }
}
