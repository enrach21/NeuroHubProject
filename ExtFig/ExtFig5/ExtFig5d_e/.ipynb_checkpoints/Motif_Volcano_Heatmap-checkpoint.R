# Packages
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(mixOmics)
library(statmod)
library(biomaRt)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)

#### Read in data ####
df <- read.table('diffTFBS_MG_OPC_bindetect_results.txt', h=T)

hocomoco <- read.csv('~/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv', sep='\t')
dim(hocomoco)
head(hocomoco, n=5)

hocomoco_filt <- hocomoco[,c('Model','Transcription.factor')]
colnames(hocomoco_filt) <- c('motif_id','Transcription.factor')
head(hocomoco_filt, n=5)

# Merge the two dataframe together
df2 <- left_join(df, hocomoco_filt, multiple = "all")
dim(df2)