
#### Packages ####

library('bedtoolsr')
library('dplyr')
# library('tidyr')
library(ggplot2)
library("UpSetR")
library("reshape2")
library('ggpointdensity')
library(viridis)

#### Functions ####
source('Count_TSS.R')
source('Plot_Corr.NH.R')
source('Count_TSS_ATAC.V3.R')


#### Read in files ####

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

MG_PLAC <- paste0(DIR,'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe')

OC_PLAC <-  paste0(DIR,'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe')

oRG_PLAC <-  paste0(DIR,'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe')

vRG_PLAC <-  paste0(DIR,'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe')


BinnedGenome <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/BinnedGenome/hg38.2kb.bed'

Tss <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed'

# Location of ATAC-seq overlapping with UMR/LMR
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/'

ATAC_MR_MG <- paste0(DIR,'MG.Both.peak.bed')

ATAC_MR_OC <- paste0(DIR,'OPC.Both.peak.bed')

ATAC_MR_oRG <- paste0(DIR,'oRG.Both.peak.bed')

ATAC_MR_vRG <- paste0(DIR,'vRG.Both.peak.bed')

anchor <- '/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/Merged/Anchor.Merged.9.7.23.final.bed'


#### Get interaction bretween genes ####

MG <- Count_TSS(Tss,BinnedGenome,MG_PLAC,'MG_H3K4')
OC <- Count_TSS(Tss,BinnedGenome,OC_PLAC,'OC_H3K4')
oRG <- Count_TSS(Tss,BinnedGenome,oRG_PLAC,'oRG_H3K4')
vRG <- Count_TSS(Tss,BinnedGenome,vRG_PLAC,'vRG_H3K4')


Plot_Corr(OC,MG,contrast=1)[2]
Plot_Corr(oRG,MG,contrast=2)[2]
Plot_Corr(vRG,MG,contrast=3)[2]
Plot_Corr(oRG,OC,contrast=4)[2]
Plot_Corr(vRG,OC,contrast=5)[2]
Plot_Corr(vRG,oRG,contrast=6)[2]

#### Correlation for only LMAR cCREs ####
MG_ATAC.MR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,MG_PLAC,anchor,'MG_H3K4',ATAC_MR_MG)
OC_ATAC.MR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,OC_PLAC,anchor,'OC_H3K4',ATAC_MR_OC)
vRG_ATAC.MR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,vRG_PLAC,anchor,'vRG_H3K4',ATAC_MR_vRG)
oRG_ATAC.MR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,oRG_PLAC,anchor,'oRG_H3K4',ATAC_MR_oRG)

Plot_Corr(OC_ATAC.MR_V3 ,MG_ATAC.MR_V3 ,contrast=1)[2]
ggsave('MG-OC.Dif.cCRE.H3K4me3.5.30.23.pdf', width = 1.75, height = 1.75)
Plot_Corr(oRG_ATAC.MR_V3 ,MG_ATAC.MR_V3 ,contrast=2)[2]
Plot_Corr(vRG_ATAC.MR_V3 ,MG_ATAC.MR_V3 ,contrast=3)[2]
Plot_Corr(oRG_ATAC.MR_V3 ,OC_ATAC.MR_V3 ,contrast=4)[2]
Plot_Corr(vRG_ATAC.MR_V3 ,OC_ATAC.MR_V3 ,contrast=5)[2]
Plot_Corr(vRG_ATAC.MR_V3 ,oRG_ATAC.MR_V3 ,contrast=6)[2]
ggsave('oRG-vRG.Dif.cCRE.H3K4me3.5.30.23.pdf', width = 1.75, height = 1.75)

#### Look at just AR cCRE ####
# Location of ATAC-seq overlapping with UMR/LMR
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/'

ATAC_MG <- paste0(DIR,'MG.ATAC.peak.bed')

ATAC_OC <- paste0(DIR,'OPC.ATAC.peak.bed')

ATAC_oRG <- paste0(DIR,'oRG.ATAC.peak.bed')

ATAC_vRG <- paste0(DIR,'vRG.ATAC.peak.bed')

MG_ATAC_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,MG_PLAC,anchor,'MG_H3K4',ATAC_MG)
OC_ATAC_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,OC_PLAC,anchor,'OC_H3K4',ATAC_OC)
vRG_ATAC_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,vRG_PLAC,anchor,'vRG_H3K4',ATAC_vRG)
oRG_ATAC_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,oRG_PLAC,anchor,'oRG_H3K4',ATAC_oRG)

Plot_Corr(OC_ATAC_V3 ,MG_ATAC_V3 ,contrast=1)[2]
Plot_Corr(oRG_ATAC_V3 ,MG_ATAC_V3 ,contrast=2)[2]
Plot_Corr(vRG_ATAC_V3 ,MG_ATAC_V3 ,contrast=3)[2]
Plot_Corr(oRG_ATAC_V3 ,OC_ATAC_V3 ,contrast=4)[2]
Plot_Corr(vRG_ATAC_V3 ,OC_ATAC_V3 ,contrast=5)[2]
Plot_Corr(vRG_ATAC_V3 ,oRG_ATAC_V3 ,contrast=6)[2]

#### Look at just LMR ####
# Location of ATAC-seq overlapping with UMR/LMR
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/'

LMR_MG <- paste0(DIR,'MG.MR.peak.bed')

LMR_OC <- paste0(DIR,'OPC.MR.peak.bed')

LMR_oRG <- paste0(DIR,'oRG.MR.peak.bed')

LMR_vRG <- paste0(DIR,'vRG.MR.peak.bed')

MG_LMR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,MG_PLAC,anchor,'MG_H3K4',LMR_MG)
OC_LMR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,OC_PLAC,anchor,'OC_H3K4',LMR_OC)
vRG_LMR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,vRG_PLAC,anchor,'vRG_H3K4',LMR_vRG)
oRG_LMR_V3 <- Count_TSS_ATAC.V3(Tss,BinnedGenome,oRG_PLAC,anchor,'oRG_H3K4',LMR_oRG)

Plot_Corr(OC_LMR_V3 ,MG_LMR_V3 ,contrast=1)[2]
Plot_Corr(oRG_LMR_V3 ,MG_LMR_V3 ,contrast=2)[2]
Plot_Corr(vRG_LMR_V3 ,MG_LMR_V3 ,contrast=3)[2]
Plot_Corr(oRG_LMR_V3 ,OC_LMR_V3 ,contrast=4)[2]
Plot_Corr(vRG_LMR_V3 ,OC_LMR_V3 ,contrast=5)[2]
Plot_Corr(vRG_LMR_V3 ,oRG_LMR_V3 ,contrast=6)[2]