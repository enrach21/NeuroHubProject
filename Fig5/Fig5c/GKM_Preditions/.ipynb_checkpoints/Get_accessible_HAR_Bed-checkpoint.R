#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)


setwd('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/cCRE_enrichment_2024')

# Read in unique ATAC-seq peaks
ReadCRE <- function(file){
    df<-read.table(file, header=F)
    colnames(df) <- c('chr','start','end')
    return(df)
}

HAR2 <- 'GSE180714_HARs.tsv'
HAR2_DF <- read.table(HAR2, header=T)
head(HAR2_DF)

ID <- 'HARsv2_2742'

HAR2_DF[HAR2_DF$HAR_ID == ID,]

test <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/oRG.4reps.overlap.optimal_peak.unique.bed'
test_V2 <- bt.intersect(a=HAR2_DF,b=test,wa=T, u=T)
colnames(test_V2) <- colnames(HAR2_DF)
test_V2

# write out bed file
write.table(test_V2[,c(1:4,8)], file='/shen/shenlabstore3/ijones1/dependencies/vcf2maf/HAR/oRG.ATAC.HARs.bed', col.names = F, row.names = F, quote = F, sep='\t')