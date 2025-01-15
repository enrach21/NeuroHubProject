library("DiffBind")
library(ggplot2)

#### Location of output ####
dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG'

setwd(dir)

#### Read in Meta data ####
samples <- read.table('diffbind.NH.oRG.vRG.txt', header=T)

#### Get samples #### 
s1 <- unique(samples$Condition)[1]
s2 <- unique(samples$Condition)[2]



#### Dba analysis ####
dba <- dba(sampleSheet=samples)
dba_count <- dba.count(dba,minOverlap=1, score = DBA_SCORE_NORMALIZED, bScaleControl = FALSE, summits=FALSE)


pdf(file=paste0(s1,'_',s2,'.NH.diffbind_count_PCA_ATAC.filt.pdf'))
dba.plotPCA(dba_count, attributes=DBA_ID)
dev.off()

pdf(file=paste0(s1,'_',s2,'.NH.diffbind_count_ATAC_9_13_23.filt.pdf'))
plot(dba_count)
dev.off()

#### 8.22.23 ####
dba_contrast <- dba.contrast(dba_count, design="~Factor + Condition", minMembers=2)
dba_normalize <- dba.normalize(dba_contrast)
dba_analysis <- dba.analyze(dba_normalize,method=DBA_ALL_METHODS)
res <- dba.report(dba_analysis, contrast=1,th=0.01)

#### Write outputs

pdf(file=paste0(s1,'_',s2,'.NH.diffbind_analyis_MA_NonNorm_ATAC.pdf'))
dba.plotMA(dba_analysis, bNormalized=FALSE, sub="Non-Normalized")
dev.off()

pdf(file=paste0(s1,'_',s2,'.NH.diffbind_analyis_MA_Norm_ATAC.pdf'))
dba.plotMA(dba_analysis , method=DBA_DESEQ2, sub="DESeq2:lib:full")
dev.off()

pdf(file=paste0(s1,'_',s2,'.diffbind_analyis_PCA_ATAC.pdf'))
dba.plotPCA(dba_analysis,contrast=1, attributes=DBA_CONDITION)
dev.off()


write.csv(res,paste0(s1,'_',s2,'diffbind_analysis.csv'))
