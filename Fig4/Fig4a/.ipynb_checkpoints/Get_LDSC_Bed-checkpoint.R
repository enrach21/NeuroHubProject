# Author: Ian Jones
# Date: 4.11.23
# Goal: Generate files for LDSC analysis

#### Read in packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')

#### Functions ####

# Read in unique ATAC-seq peaks
ReadATAC <- function(dir, file){
    df<-read.table(paste0(dir,file))
    colnames(df) <- c('chr','start','end')
    return(df)
}

# Read in unique methylation peaks
ReadMRs <- function(dir, file){
    df<-read.table(paste0(dir,file), header=T)
    return(df)
}


# Read in Data
# Location of called peaks
ATAC_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/'

# Read in Data
MG_ATAC <- ReadATAC(ATAC_dir,'MG.overlap.optimal_peak.unique.bed')
nrow(MG_ATAC)
OPC_ATAC <- ReadATAC(ATAC_dir,'OPC.overlap.optimal_peak.unique.bed')
nrow(OPC_ATAC)
oRG_ATAC <- ReadATAC(ATAC_dir,'oRG.4reps.overlap.optimal_peak.unique.bed')
nrow(oRG_ATAC)
vRG_ATAC <- ReadATAC(ATAC_dir,'vRG.4reps.overlap.optimal_peak.unique.bed')
nrow(vRG_ATAC)

# Combine oRG_ATAC and vRG_ATAC
RG_ATAC <- rbind(oRG_ATAC, vRG_ATAC)
RG_ATAC <- bt.sort(RG_ATAC)
RG_ATAC <- bt.merge(RG_ATAC)
nrow(RG_ATAC)

# Location of called peaks
MRs_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie/'

MG_MRs <- ReadMRs(MRs_dir,'MG_UMRsLMRs.tab')
nrow(MG_MRs)
OPC_MRs <- ReadMRs(MRs_dir,'OPC_UMRsLMRs.tab')
nrow(OPC_MRs)
oRG_MRs <- ReadMRs(MRs_dir,'oRG_UMRsLMRs.tab')
nrow(oRG_MRs)
vRG_MRs <- ReadMRs(MRs_dir,'vRG_UMRsLMRs.tab')
nrow(vRG_MRs)

# Combine oRG_ATAC and vRG_ATAC
RG_MRs <- rbind(oRG_MRs, vRG_MRs)
RG_MRs <- bt.sort(RG_MRs)
RG_MRs <- bt.merge(RG_MRs)
nrow(RG_MRs)

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC-seq/BEDPE/'

OC <- 'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

MG <- 'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'


# Read in 
OC_DF <- read.table(paste0(DIR,OC), fill=T, header=T)
length(unique(OC_DF$InteractionName))

MG_DF <- read.table(paste0(DIR,MG), fill=T, header=T)
length(unique(MG_DF$InteractionName))

vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)
length(unique(vRG_DF$InteractionName))

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)
length(unique(oRG_DF$InteractionName))

RG_DF <- full_join(vRG_DF, oRG_DF)
length(unique(oRG_DF$InteractionName))


#### Write out LDSC Files ####

# OPC
# XOR
DF <- OC_DF
ATAC <- OPC_ATAC
MRs <- OPC_MRs
Cell <- 'OC'

df_XOR <- DF[DF$type == 'XOR',]
df_AND <- DF[DF$type == 'AND',]


# Get XOR ATAC-seq peaks
df_XOR_ATAC <- bt.intersect(ATAC,df_XOR, wa=T, u=T)
nrow(df_XOR_ATAC)
write.table(df_XOR_ATAC, paste0(Cell,'.XOR.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_XOR_MRs <- bt.intersect( MRs,df_XOR, wa=T, u=T)
nrow(df_XOR_MRs)
write.table(df_XOR_MRs, paste0(Cell,'.XOR.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_XOR_ATAC_MRs <- bt.intersect(df_XOR_ATAC, df_XOR_MRs)
nrow(df_XOR_ATAC_MRs)
write.table(df_XOR_ATAC_MRs, paste0(Cell,'.XOR.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_AND_ATAC <- bt.intersect( ATAC,df_AND, wa=T, u=T)
nrow(df_AND_ATAC)
write.table(df_AND_ATAC, paste0(Cell,'.AND.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_AND_MRs <- bt.intersect(MRs,df_AND, wa=T, u=T)
nrow(df_AND_MRs)
write.table(df_AND_MRs, paste0(Cell,'.AND.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')



df_AND_ATAC_MRs <- bt.intersect(df_AND_ATAC, df_AND_MRs)
nrow(df_AND_ATAC_MRs)
write.table(df_AND_ATAC_MRs, paste0(Cell,'.AND.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


# MG
# XOR
DF <- MG_DF
ATAC <- MG_ATAC
MRs <- MG_MRs
Cell <- 'MG'

df_XOR <- DF[DF$type == 'XOR',]
df_AND <- DF[DF$type == 'AND',]


# Get XOR ATAC-seq peaks
df_XOR_ATAC <- bt.intersect(ATAC,df_XOR, wa=T, u=T)
nrow(df_XOR_ATAC)
write.table(df_XOR_ATAC, paste0(Cell,'.XOR.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_XOR_MRs <- bt.intersect( MRs,df_XOR, wa=T, u=T)
nrow(df_XOR_MRs)
write.table(df_XOR_MRs, paste0(Cell,'.XOR.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_XOR_ATAC_MRs <- bt.intersect(df_XOR_ATAC, df_XOR_MRs)
nrow(df_XOR_ATAC_MRs)
write.table(df_XOR_ATAC_MRs, paste0(Cell,'.XOR.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_AND_ATAC <- bt.intersect( ATAC,df_AND, wa=T, u=T)
nrow(df_AND_ATAC)
write.table(df_AND_ATAC, paste0(Cell,'.AND.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_AND_MRs <- bt.intersect(MRs,df_AND, wa=T, u=T)
nrow(df_AND_MRs)
write.table(df_AND_MRs, paste0(Cell,'.AND.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')



df_AND_ATAC_MRs <- bt.intersect(df_AND_ATAC, df_AND_MRs)
nrow(df_AND_ATAC_MRs)
write.table(df_AND_ATAC_MRs, paste0(Cell,'.AND.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

# oRG
# XOR
DF <- oRG_DF
ATAC <- oRG_ATAC
MRs <- oRG_MRs
Cell <- 'oRG'

df_XOR <- DF[DF$type == 'XOR',]
df_AND <- DF[DF$type == 'AND',]


# Get XOR ATAC-seq peaks
df_XOR_ATAC <- bt.intersect(ATAC,df_XOR, wa=T, u=T)
nrow(df_XOR_ATAC)
write.table(df_XOR_ATAC, paste0(Cell,'.XOR.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_XOR_MRs <- bt.intersect( MRs,df_XOR, wa=T, u=T)
nrow(df_XOR_MRs)
write.table(df_XOR_MRs, paste0(Cell,'.XOR.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_XOR_ATAC_MRs <- bt.intersect(df_XOR_ATAC, df_XOR_MRs)
nrow(df_XOR_ATAC_MRs)
write.table(df_XOR_ATAC_MRs, paste0(Cell,'.XOR.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_AND_ATAC <- bt.intersect( ATAC,df_AND, wa=T, u=T)
nrow(df_AND_ATAC)
write.table(df_AND_ATAC, paste0(Cell,'.AND.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_AND_MRs <- bt.intersect(MRs,df_AND, wa=T, u=T)
nrow(df_AND_MRs)
write.table(df_AND_MRs, paste0(Cell,'.AND.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')



df_AND_ATAC_MRs <- bt.intersect(df_AND_ATAC, df_AND_MRs)
nrow(df_AND_ATAC_MRs)
write.table(df_AND_ATAC_MRs, paste0(Cell,'.AND.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

# vRG 
# XOR
DF <- vRG_DF
ATAC <- vRG_ATAC
MRs <- vRG_MRs
Cell <- 'vRG'

df_XOR <- DF[DF$type == 'XOR',]
df_AND <- DF[DF$type == 'AND',]


# Get XOR ATAC-seq peaks
df_XOR_ATAC <- bt.intersect(ATAC,df_XOR, wa=T, u=T)
nrow(df_XOR_ATAC)
write.table(df_XOR_ATAC, paste0(Cell,'.XOR.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_XOR_MRs <- bt.intersect( MRs,df_XOR, wa=T, u=T)
nrow(df_XOR_MRs)
write.table(df_XOR_MRs, paste0(Cell,'.XOR.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_XOR_ATAC_MRs <- bt.intersect(df_XOR_ATAC, df_XOR_MRs)
nrow(df_XOR_ATAC_MRs)
write.table(df_XOR_ATAC_MRs, paste0(Cell,'.XOR.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')


df_AND_ATAC <- bt.intersect( ATAC,df_AND, wa=T, u=T)
nrow(df_AND_ATAC)
write.table(df_AND_ATAC, paste0(Cell,'.AND.ATAC.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

df_AND_MRs <- bt.intersect(MRs,df_AND, wa=T, u=T)
nrow(df_AND_MRs)
write.table(df_AND_MRs, paste0(Cell,'.AND.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')



df_AND_ATAC_MRs <- bt.intersect(df_AND_ATAC, df_AND_MRs)
nrow(df_AND_ATAC_MRs)
write.table(df_AND_ATAC_MRs, paste0(Cell,'.AND.ATAC.MRs.bed'), col.names = F, row.names = F, quote = F, sep = '\t')