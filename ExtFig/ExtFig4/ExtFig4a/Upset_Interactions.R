#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(ggplot2)
library("UpSetR")
library("reshape2")

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

# Paste Read one and Read two together
Make_Int_List <- function(DF) {
    TempList <- paste0(DF[,1],DF[,2],DF[,3],DF[,4],DF[,5],DF[,6])
    return(TempList)
}

# Get input for upset

OC_UPSET_DF <- read.table(paste0(DIR,'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'), fill=T, header=T)[,1:6]

MG_UPSET_DF <- read.table(paste0(DIR,'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'), fill=T, header=T)[,1:6]

vRG_UPSET_DF <- read.table(paste0(DIR,'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'), fill=T, header=T)[,1:6]

oRG_UPSET_DF <- read.table(paste0(DIR,'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'), fill=T, header=T)[,1:6]


OC_UPSET_LIST <- Make_Int_List(OC_UPSET_DF)
MG_UPSET_LIST <- Make_Int_List(MG_UPSET_DF)
vRG_UPSET_LIST <- Make_Int_List(vRG_UPSET_DF)
oRG_UPSET_LIST <- Make_Int_List(oRG_UPSET_DF)


listInput <- list(OCP = OC_UPSET_LIST , MG = MG_UPSET_LIST , vRG=vRG_UPSET_LIST, oRG=oRG_UPSET_LIST)


pdf('PLAC.Upset.5.31.24.pdf', height=3, width=4.5)
upset(fromList(listInput), order.by = "freq")
dev.off()