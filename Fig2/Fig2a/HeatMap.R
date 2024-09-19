# Author: Ian Jones
# Date: 7-13-23
# Goal: Generate heatmap looking into dynamic changes in PLAC-seq contacts


#### Packages ####
library(data.table)
library(dplyr)
library('preprocessCore')
library('bedtoolsr')
library(pheatmap)
library(ggplot2)
library(ggplotify)
library(ComplexHeatmap)
library(patchwork)
library('tidyr')
library(data.table)
library(circlize)

#### Functions ####

# Function to read in modified PLAC-seq and select certain columns
source('/shen/shenlabstore3/ijones1/ShenAnalysis/iNGN_TC/Scripts/PLAC-seq/HeatMap.2.24.23/Get.XOR.R')

# Function to obtain normalized contact for non-sig interactions
source('/shen/shenlabstore3/ijones1/ShenAnalysis/iNGN_TC/Scripts/PLAC-seq/HeatMap.2.24.23/Fill.NonSig.R')

#### Read in PLAC-seq files ####

# Location of PLAC-seq files

DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

OC <- 'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

MG <- 'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

# Get XOR interactions
OC_XOR <- Get.XOR(paste0(DIR,OC), 'OC')
MG_XOR <- Get.XOR(paste0(DIR,MG), 'MG')
vRG_XOR <- Get.XOR(paste0(DIR,vRG), 'vRG')
oRG_XOR <- Get.XOR(paste0(DIR,oRG), 'oRG')

# Combine files
Merged_XOR <- full_join(OC_XOR , MG_XOR)
Merged_XOR  <- full_join(Merged_XOR , oRG_XOR)
Merged_XOR  <- full_join(Merged_XOR , vRG_XOR)
dim(Merged_XOR)

# Add labels for interactions
Merged_XOR$lab <- 'Other'

# Label oRG specific
Merged_XOR$lab[is.na(Merged_XOR$MG_norm)] <- 'Not_MG'

# Label OC specific
Merged_XOR$lab[is.na(Merged_XOR$MG_norm) & is.na(Merged_XOR$vRG_norm) & is.na(Merged_XOR$oRG_norm)] <- 'OC'

 # Label MG specific
Merged_XOR$lab[is.na(Merged_XOR$OC_norm) & is.na(Merged_XOR$vRG_norm) & is.na(Merged_XOR$oRG_norm)] <- 'MG'

# Label oRG specific
Merged_XOR$lab[is.na(Merged_XOR$OC_norm) & is.na(Merged_XOR$MG_norm)] <- 'RG'

# Label vRG specific
Merged_XOR$lab[is.na(Merged_XOR$OC_norm) & is.na(Merged_XOR$oRG_norm) & is.na(Merged_XOR$MG_norm)] <- 'vRG'
                                                
                                                 
# Label oRG specific
Merged_XOR$lab[is.na(Merged_XOR$OC_norm) & is.na(Merged_XOR$vRG_norm) & is.na(Merged_XOR$MG_norm)] <- 'oRG'

table(Merged_XOR$lab)
head(Merged_XOR)

#### Fill in Missing Data ####

chr_list <- c('chr1','chr2','chr3','chr4','chr5','chr6',
             'chr7','chr8','chr9','chr10','chr11','chr12',
             'chr13','chr14','chr15','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22')

# filter to only be autosomes
dim(Merged_XOR)
Merged_XOR <- Merged_XOR[Merged_XOR$chr %in% chr_list,]
dim(Merged_XOR)


# Get non-sig contacts for OC
Merged_XOR <- Fill.NonSig(chr_list,'/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/MAPS_output/Oligo_100bp.fastp.60mil.FRIP_20230908_103229/reg_raw.',
           '.Oligo_100bp.fastp.60mil.FRIP.2k.xor.MAPS2_pospoisson','OC',2000,Merged_XOR)

# Get non-sig contacts for MG
Merged_XOR <- Fill.NonSig(chr_list,'/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/MAPS_output/Microglia_100bp.fastp.60mil.FRIP_20230908_103229/reg_raw.',
           '.Microglia_100bp.fastp.60mil.FRIP.2k.xor.MAPS2_pospoisson','MG',2000,Merged_XOR)

# Get non-sig contacts for vRG
Merged_XOR <- Fill.NonSig(chr_list,'/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/MAPS_output/vRG_100bp.fastp.60mil.FRIP_20230908_103229/reg_raw.',
           '.vRG_100bp.fastp.60mil.FRIP.2k.xor.MAPS2_pospoisson','vRG',2000,Merged_XOR)

# Get non-sig contacts for oRG
Merged_XOR <- Fill.NonSig(chr_list,'/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/MAPS_output/oRG_100bp.fastp.60mil.FRIP_20230908_103229/reg_raw.',
           '.oRG_100bp.fastp.60mil.FRIP.2k.xor.MAPS2_pospoisson','oRG',2000,Merged_XOR)

Merged_XOR[is.na(Merged_XOR)] <- 0
head(Merged_XOR)
dim(Merged_XOR)
table(Merged_XOR$lab)

write.table(Merged_XOR, 'Merged_XOR.11.21.23.txt')

#### Read in ATAC-seq data ####

OC_ATAC <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Heatmap/OPC.11.21.23.overlap.optimal_peak.count.bed')
colnames(OC_ATAC) <- c('chr','start','stop','peak','OC_mean')
dim(OC_ATAC)

MG_ATAC <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Heatmap/MG.11.21.23.overlap.optimal_peak.count.bed')
colnames(MG_ATAC) <- c('chr','start','stop','peak','MG_mean')
dim(MG_ATAC)

vRG_ATAC <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Heatmap/vRG.4reps.8.15.23.overlap.optimal_peak.count.bed')
colnames(vRG_ATAC) <- c('chr','start','stop','peak','vRG_mean')
dim(vRG_ATAC)

oRG_ATAC <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Heatmap/oRG.4reps.8.15.23.overlap.optimal_peak.count.bed')
colnames(oRG_ATAC) <- c('chr','start','stop','peak','oRG_mean')
dim(oRG_ATAC)

# combine ATAC
Merge_ATAC <- full_join(OC_ATAC, MG_ATAC)
Merge_ATAC <- full_join(Merge_ATAC, vRG_ATAC)
Merge_ATAC <- full_join(Merge_ATAC, oRG_ATAC)
dim(Merge_ATAC)
head(Merge_ATAC)

# Perform Quantile normalization
Merge_ATAC[,5:8] <- normalize.quantiles(as.matrix(Merge_ATAC[,5:8]))
head(Merge_ATAC)

# Overlap peaks with 5kb bin

# Read in 5kb genome
BinnedGenome <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/BinnedGenome/hg38.2kb.bed'
BinnedGenomeDF <- read.table(BinnedGenome)
colnames(BinnedGenomeDF) <- c('bin_chr','bin_start','bin_end')
# BinnedGenomeDF$bin_end <- BinnedGenomeDF$bin_end-1
head(BinnedGenomeDF)

# Overlap with ATAC-seq with at least 50% overlap
Binned_Merge_ATAC <- bt.intersect(Merge_ATAC, BinnedGenomeDF,wa=T,wb=T, f=0.50)
colnames(Binned_Merge_ATAC) <- c(colnames(Merge_ATAC),colnames(BinnedGenomeDF))


Binned_Merge_ATAC <- Binned_Merge_ATAC[,c('bin_chr','bin_start','bin_end',
                                          'OC_mean','MG_mean','vRG_mean','oRG_mean')]

# Get the sum within 2kb bins
Binned_Merge_ATAC_Sum <- Binned_Merge_ATAC %>% group_by(bin_chr,bin_start,bin_end) %>% 
    summarise(OC_ATAC_SUM=sum(OC_mean),
             MG_ATAC_SUM=sum(MG_mean),
            vRG_ATAC_SUM=sum(vRG_mean),
             oRG_ATAC_SUM=sum(oRG_mean))

# Scale from 0-1
Binned_Merge_ATAC_Sum$OC_ATAC_ratio <- Binned_Merge_ATAC_Sum$OC_ATAC_SUM / (Binned_Merge_ATAC_Sum$OC_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$MG_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$vRG_ATAC_SUM +
                                                                              Binned_Merge_ATAC_Sum$oRG_ATAC_SUM)

Binned_Merge_ATAC_Sum$MG_ATAC_ratio <- Binned_Merge_ATAC_Sum$MG_ATAC_SUM / (Binned_Merge_ATAC_Sum$OC_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$MG_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$vRG_ATAC_SUM +
                                                                              Binned_Merge_ATAC_Sum$oRG_ATAC_SUM)


Binned_Merge_ATAC_Sum$vRG_ATAC_ratio <- Binned_Merge_ATAC_Sum$vRG_ATAC_SUM / (Binned_Merge_ATAC_Sum$OC_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$MG_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$vRG_ATAC_SUM +
                                                                              Binned_Merge_ATAC_Sum$oRG_ATAC_SUM)

Binned_Merge_ATAC_Sum$oRG_ATAC_ratio <- Binned_Merge_ATAC_Sum$oRG_ATAC_SUM / (Binned_Merge_ATAC_Sum$OC_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$MG_ATAC_SUM +
                                                                             Binned_Merge_ATAC_Sum$vRG_ATAC_SUM +
                                                                              Binned_Merge_ATAC_Sum$oRG_ATAC_SUM)

# change col names so that it can merge with PLAC-seq data
colnames(Binned_Merge_ATAC_Sum)[1:3] <- c('chr','start','stop')

#### Combine ATAC-seq and PLAC-seq ####

combineDF <- inner_join(Merged_XOR, Binned_Merge_ATAC_Sum)
dim(Merged_XOR)
dim(combineDF)

write.table(combineDF, 'Merged_XOR.ATAC.11.23.23.txt')

#### Read in Methylation ####

# Read in UMR/LMR
Dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/'

MG_MR <- read.table(paste0(Dir,'MG_UMRsLMRs.tab'), header=T)
head(MG_MR)

OC_MR <- read.table(paste0(Dir,'OPC_UMRsLMRs.tab'), header=T)
head(OC_MR)

oRG_MR <- read.table(paste0(Dir,'oRG_UMRsLMRs.tab'), header=T)
head(oRG_MR)

vRG_MR <- read.table(paste0(Dir,'vRG_UMRsLMRs.tab'), header=T)
head(vRG_MR)


MR_df <- rbind(MG_MR,OC_MR)
MR_df <- rbind(MR_df,oRG_MR)
MR_df <- rbind(MR_df,vRG_MR)
dim(MR_df)


CpG_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Methylation'

# Read in MG
MG_CpG <- read.table(paste0(CpG_dir,'/MG/ATAC_binned2_MG_CpG.bismark.10.cov.bed'))
colnames(MG_CpG) <- c('chr','start','stop','MG_Methyl','MG_Unmethyl','MG_count','MG_CpG_Percent')
head(MG_CpG )
dim(MG_CpG)

# Read in MG
OC_CpG <- read.table(paste0(CpG_dir,'/OC/ATAC_binned2_OPC_CpG.bismark.10.cov.bed'))
colnames(OC_CpG) <-  c('chr','start','stop','OC_Methyl','OC_Unmethyl','OC_count','OC_CpG_Percent')
head(OC_CpG)
dim(OC_CpG)

# Read in MG
oRG_CpG <- read.table(paste0(CpG_dir,'/oRG/ATAC_binned2_oRG_CpG.bismark.10.cov.bed'))
colnames(oRG_CpG) <-  c('chr','start','stop','oRG_Methyl','oRG_Unmethyl','oRG_count','oRG_CpG_Percent')
head(oRG_CpG )
dim(oRG_CpG )

# Read in MG
vRG_CpG <- read.table(paste0(CpG_dir,'/vRG/ATAC_binned2_vRG_CpG.bismark.10.cov.bed'))
colnames(vRG_CpG) <-  c('chr','start','stop','vRG_Methyl','vRG_Unmethyl','vRG_count','vRG_CpG_Percent')
head(vRG_CpG)
dim(vRG_CpG)

# combine CpG
Merge_CpG <- inner_join(MG_CpG, OC_CpG )
Merge_CpG  <- inner_join(Merge_CpG , oRG_CpG )
Merge_CpG <- inner_join(Merge_CpG , vRG_CpG)
dim(Merge_CpG )
head(Merge_CpG )


# Overlap with UMRs/LMRs
Merge_CpG_MR <- bt.intersect(a= Merge_CpG, b =MR_df, wa=T, u=T)
colnames(Merge_CpG_MR) <- colnames(Merge_CpG)
head(Merge_CpG_MR)
dim(Merge_CpG_MR)


Merge_CpG <- Merge_CpG_MR[,c('chr','start','stop',
                          'MG_CpG_Percent','OC_CpG_Percent',
                         'oRG_CpG_Percent','vRG_CpG_Percent')]
head(Merge_CpG)
dim(Merge_CpG)

# Overlap peaks with 2kb bin

# Read in 2kb genome
BinnedGenome <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/BinnedGenome/hg38.2kb.bed'
BinnedGenomeDF <- read.table(BinnedGenome)
colnames(BinnedGenomeDF) <- c('bin_chr','bin_start','bin_end')
head(BinnedGenomeDF)

# Overlap with ATAC-seq with at least 50% overlap
Binned_Merge_CpG <- bt.intersect(Merge_CpG, BinnedGenomeDF,wa=T,wb=T, f=0.50)
colnames(Binned_Merge_CpG) <- c(colnames(Merge_CpG),colnames(BinnedGenomeDF))
dim(Binned_Merge_CpG)

Binned_Merge_CpG <- Binned_Merge_CpG[,c('bin_chr','bin_start','bin_end',
                                          'MG_CpG_Percent','OC_CpG_Percent',
                         'oRG_CpG_Percent','vRG_CpG_Percent')]
colnames(Binned_Merge_CpG)[1:3] <- c('chr','start','stop')
head(Binned_Merge_CpG)
dim(Binned_Merge_CpG)
dim(unique(Binned_Merge_CpG[,1:3]))

# Get average of groups
dim(Binned_Merge_CpG )
Binned_Merge_CpG <- Binned_Merge_CpG %>%
    group_by(chr,start,stop) %>%
    dplyr::summarize(MG_CpG_Percent = mean(MG_CpG_Percent, na.rm=TRUE),
                    OC_CpG_Percent = mean(OC_CpG_Percent, na.rm=TRUE),
                    oRG_CpG_Percent = mean(oRG_CpG_Percent, na.rm=TRUE),
                    vRG_CpG_Percent = mean(vRG_CpG_Percent, na.rm=TRUE))
dim(Binned_Merge_CpG )

combineDF <- inner_join(combineDF, Binned_Merge_CpG)
# dim(Merged_XOR)
dim(combineDF)


write.table(combineDF, 'Merged_XOR.ATAC.CpG.11.21.23.txt')


#### Read in RNA ####

# Get RNA
RNA_df <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv')[,-1]
RNA_df <- unique(RNA_df)
RNA_df <- unique(RNA_df[,c('oRG','vRG','MG','Oligo','hgnc_symbol')])
colnames(RNA_df ) <- c('oRG','vRG','MG','OC','TargetGene')
dim(RNA_df)
head(RNA_df)

# Get PLAC Target Genes
Merged_temp <- combineDF[,c('TargetAnchor','TargetGene')]
dim(Merged_temp)
head(Merged_temp)

# Split by gene name
Merged_temp_sep <- separate_rows(Merged_temp, 'TargetGene', sep='\\|')
dim(Merged_temp_sep )
head(Merged_temp_sep )

# merge with expression
Merged_RNA <- unique(inner_join(Merged_temp_sep,RNA_df,multiple = "first"))
dim(Merged_RNA)
head(Merged_RNA)

# Perform log2 transformation
Merged_RNA[,3:6] <- log2(Merged_RNA[,3:6]+1)
head(Merged_RNA)

# Sum the RNA expression
Merge_RNA_Sum <- Merged_RNA %>% group_by(TargetAnchor) %>% 
    summarise(sum_oRG_TMM_RPKM=sum(oRG),
             sum_vRG_TMM_RPKM=sum(vRG),
             sum_MG_TMM_RPKM=sum(MG),
             sum_OC_TMM_RPKM=sum(OC))
head(Merge_RNA_Sum)


# scale out of 0-1
Merge_RNA_Sum$oRG_RNA_Ratio <- Merge_RNA_Sum$sum_oRG_TMM_RPKM / (Merge_RNA_Sum$sum_oRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_vRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_MG_TMM_RPKM +
                                                                              Merge_RNA_Sum$sum_OC_TMM_RPKM)

Merge_RNA_Sum$vRG_RNA_Ratio <- Merge_RNA_Sum$sum_vRG_TMM_RPKM / (Merge_RNA_Sum$sum_oRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_vRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_MG_TMM_RPKM +
                                                                              Merge_RNA_Sum$sum_OC_TMM_RPKM)

Merge_RNA_Sum$MG_RNA_Ratio <- Merge_RNA_Sum$sum_MG_TMM_RPKM / (Merge_RNA_Sum$sum_oRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_vRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_MG_TMM_RPKM +
                                                                              Merge_RNA_Sum$sum_OC_TMM_RPKM)

Merge_RNA_Sum$OC_RNA_Ratio <- Merge_RNA_Sum$sum_OC_TMM_RPKM / (Merge_RNA_Sum$sum_oRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_vRG_TMM_RPKM +
                                                                             Merge_RNA_Sum$sum_MG_TMM_RPKM +
                                                                              Merge_RNA_Sum$sum_OC_TMM_RPKM)


# Combine with other data
combineDF2 <- inner_join(combineDF, Merge_RNA_Sum)


# Remove rows not targeting any gene
combineDF2 <- combineDF2[!is.na(combineDF2$oRG_RNA_Ratio),]

write.table(combineDF2, 'Merged_XOR.ATAC.CpG.RNA.5.30.24.txt')


dim(combineDF2)
# dim(combineDF)
table(combineDF2$lab)


#### Obtain unqiue interactions with LMAR cCRE ####

# vRG specfiic

vRG_DF3 <- bt.intersect(a=combineDF2[combineDF2$lab == 'vRG',], 
             b=list('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed'),
            wa=T,
            u=T,
            f=0.5)
colnames(vRG_DF3) <- colnames(combineDF2)
dim(vRG_DF3)


oRG_DF3 <- bt.intersect(a=combineDF2[combineDF2$lab == 'oRG',], 
             b=list('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed'),
            wa=T,
            u=T,
            f=0.5)
colnames(oRG_DF3) <- colnames(combineDF2)
dim(oRG_DF3)

MG_DF3 <- bt.intersect(a=combineDF2[combineDF2$lab == 'MG',], 
             b=list('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.Both.peak.bed'),
            wa=T,
            u=T,
            f=0.5)
colnames(MG_DF3) <- colnames(combineDF2)
dim(MG_DF3)

OC_DF3 <- bt.intersect(a=combineDF2[combineDF2$lab == 'OC',], 
             b=list('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.Both.peak.bed'),
            wa=T,
            u=T,
            f=0.5)
colnames(OC_DF3) <- colnames(combineDF2)
dim(OC_DF3)


#### Plot each cell types ####

plot_heatmap <- function(combineDF_filt, name, height) {
    
    # Change factor for the order
    combineDF_filt$lab <- factor(combineDF_filt$lab, levels = rev(c('MG','OC','oRG','vRG')))

        # Reorder
    combineDF_filt$mean <- rowMeans(combineDF_filt[,rev(c('MG_norm','OC_norm','oRG_norm','vRG_norm'))]
                                    , na.rm=TRUE)
    combineDF_sort <- combineDF_filt[(order(combineDF_filt$lab,combineDF_filt$vRG_norm, decreasing = TRUE)),]

    # Get mat for each condition

    # PLAC
    PLAC <- as.matrix(combineDF_sort[,rev(c('MG_norm','OC_norm','oRG_norm', 'vRG_norm'))])
    dim(PLAC)
    # ATAC
    ATAC <- as.matrix(combineDF_sort[,rev(c('MG_ATAC_ratio','OC_ATAC_ratio','oRG_ATAC_ratio','vRG_ATAC_ratio'))])
    dim(ATAC)
    # CpG
    CpG <- as.matrix(combineDF_sort[,rev(c('MG_CpG_Percent','OC_CpG_Percent','oRG_CpG_Percent','vRG_CpG_Percent'))])
    dim(CpG)
    # RNA
    RNA <- as.matrix(combineDF_sort[,rev(c('MG_RNA_Ratio','OC_RNA_Ratio','oRG_RNA_Ratio','vRG_RNA_Ratio'))])
    dim(RNA)

    #### Heatmap for PLAC ####
    mat <- PLAC

    col_fun = colorRamp2(c(0,2.5, 5), c("#2166AC","#F7F7F7", "#B2182B"))

    Mean_Sig <- colMeans(mat)

    column_ha = HeatmapAnnotation(MeanSig = anno_simple(Mean_Sig, height = unit(2, "mm"), col =  col_fun, border = TRUE),show_legend = FALSE)

    H1 <- Heatmap(mat, 
    # top_annotation = ha,
            # set orderof the rows
            row_order =  rownames(mat),
            # set order of the columns
            column_order = colnames(mat),
            # bolean to show rowname
            show_row_names=FALSE,
            # change the color format
            col = col_fun,
            show_column_names=FALSE,
            # Split rows based on interaction
            # row_split = combineDF_sort$lab, 

            # Split rows based on interaction for fitness regions
            row_split = combineDF_sort$lab, 

            # Increase gap size
            row_gap = unit(2.5, "mm"),
            # Change the label of the rows
            # row_title = c('n=5,261','n=5,295','n=22,078'),
            # change Row names 90
            row_title_rot = 0,
            # Change size and color of rows
            row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
            # Add title to describe the data
            column_title = 'Norm. Contact Freq.',
            column_title_gp = gpar(fontsize = 8),
            #Change Legend
            show_heatmap_legend = FALSE, 
            top_annotation = column_ha, border = TRUE   
            )

    lgd1 =  Legend(col_fun = col_fun, title = "Norm Contact",
            direction = "horizontal",
            legend_width = unit(3, "cm"))

    #### Heatmap for ATAC ####
    mat <- ATAC

    col_fun = colorRamp2(c(0,.3, 0.5), c("#2166AC","#F7F7F7", "#B2182B"))

    Mean_Sig <- colMeans(mat)

    column_ha = HeatmapAnnotation(MeanSig = anno_simple(Mean_Sig, height = unit(2, "mm"), col =  col_fun, border = TRUE),show_legend = FALSE)

    H2 <- Heatmap(mat,
    #         top_annotation = ha,
            # set order of the rows
            row_order =  rownames(mat),
            # set order of the columns
            column_order = colnames(mat),
            # bolean to show rowname
            show_row_names=FALSE,
            show_column_names=FALSE,
            # change the color format
            col = col_fun,
            # Split rows based on interaction
            # row_split = combineDF_sort$lab, 

            # Split rows based on interaction for fitness regions
            row_split = combineDF_sort$lab, 

            # Increase gap size
            row_gap = unit(2.5, "mm"),
            # Change the label of the rows
            # row_title = c('n=6,237','n=6,005','n=23,450'),
            # change Row names 90
            row_title_rot = 0,
            # Change size and color of rows
            row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
            # Change col name and location
            # column_names_side = "none" ,
            # column_names_rot = 0,
            # column_names_centered = TRUE,
            # column_names_gp = gpar(fill = c("#332288", "#117733","#882255"),
            #                        col = c("white", "white","white"),
            # fontsize = c(3,3,3)),
            # Add title to describe the data
            column_title = 'Relative Accessibility',
            column_title_gp = gpar(fontsize = 8),
            #Change Legend
            show_heatmap_legend = FALSE, 
            top_annotation = column_ha , border = TRUE  )

    lgd2 =  Legend(col_fun = col_fun, title = "% ATAC",
            direction = "horizontal",
            legend_width = unit(3, "cm"))


    #### Heatmap for ATAC ####
    mat <- CpG

    col_fun = colorRamp2(c(0,.5, 1), rev(c("#2166AC","#F7F7F7", "#B2182B")))

    Mean_Sig <- colMeans(mat)

    column_ha = HeatmapAnnotation(MeanSig = anno_simple(Mean_Sig, height = unit(2, "mm"), col =  col_fun, border = TRUE),show_legend = FALSE)

    H3 <- Heatmap(mat,
    #         top_annotation = ha,
            # set order of the rows
            row_order =  rownames(mat),
            # set order of the columns
            column_order = colnames(mat),
            # bolean to show rowname
            show_row_names=FALSE,
            show_column_names=FALSE,
            # change the color format
            col = col_fun,
            # Split rows based on interaction
            # row_split = combineDF_sort$lab, 

            # Split rows based on interaction for fitness regions
            row_split = combineDF_sort$lab, 

            # Increase gap size
            row_gap = unit(2.5, "mm"),
            # Change the label of the rows
            # row_title = c('n=6,237','n=6,005','n=23,450'),
            # change Row names 90
            row_title_rot = 0,
            # Change size and color of rows
            row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
            # Change col name and location
            # column_names_side = "none" ,
            # column_names_rot = 0,
            # column_names_centered = TRUE,
            # column_names_gp = gpar(fill = c("#332288", "#117733","#882255"),
            #                        col = c("white", "white","white"),
            # fontsize = c(3,3,3)),
            # Add title to describe the data
            column_title = '% of Methylated CpG',
            column_title_gp = gpar(fontsize = 8),
            #Change Legend
            show_heatmap_legend = FALSE, 
            top_annotation = column_ha, border = TRUE)

    lgd3 =  Legend(col_fun = col_fun, title = "% CpG",
            direction = "horizontal",
            legend_width = unit(3, "cm"))

    #### Heatmap for RNA ####
    mat <- RNA

    col_fun = colorRamp2(c(0,.25, .4), c("#2166AC","#F7F7F7", "#B2182B"))

    Mean_Sig <- colMeans(mat)

    column_ha = HeatmapAnnotation(MeanSig = anno_simple(Mean_Sig, height = unit(2, "mm"), col =  col_fun, border = TRUE),show_legend = FALSE)


    H4 <- Heatmap(mat, 
            # top_annotation = ha,
            # set order of the rows
            row_order =  rownames(mat),
            # set order of the columns
            column_order = colnames(mat),
            # bolean to show rowname
            show_row_names=FALSE,
            show_column_names=FALSE,
            # change the color format
            col = col_fun,
            # Split rows based on interaction
            # row_split = combineDF_sort$lab, 

            # Split rows based on interaction for fitness regions
            row_split = combineDF_sort$lab, 


            # Increase gap size
            row_gap = unit(2.5, "mm"),
            # Change the label of the rows
            # row_title = c('n=6,237','n=6,005','n=23,450'),
            # change Row names 90
            row_title_rot = 0,
            # Change size and color of rows
            row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
            # Change col name and location
            # column_names_side = "none" ,
            # column_names_rot = 0,
            # column_names_centered = TRUE,
            # column_names_gp = gpar(fill = c("#332288", "#117733","#882255"),
            #                        col = c("white", "white","white"),
           #  fontsize = c(3,3,3)),
            # Add title to describe the data
            column_title = 'Relative Expression',
            column_title_gp = gpar(fontsize = 8),
            #Change Legend
            show_heatmap_legend = FALSE, 
            top_annotation = column_ha  , border = TRUE)

    lgd4 =  Legend(col_fun = col_fun, title = "% RNA",
            direction = "horizontal",
            legend_width = unit(3, "cm"))

        # Combine plots
    ht_list = H1 + H2 + H3 + H4
    options(repr.plot.width=4.5, repr.plot.height=height)
    draw(ht_list)

    # Combine plots
    ht_list = H1 + H2 + H3 +H4

    pdf(paste0(name,'.LMAR.Heatmap.pdf'), height = height, width = 4.5)
    draw(ht_list)
    dev.off()

    pd = packLegend(lgd1, lgd2, lgd3, lgd4, direction = "horizontal")
    pdf('NH.Total.legend.pdf', height = 1, width = 6)
    draw(pd)
    dev.off()
}


plot_heatmap(vRG_DF3, 'vRG', 0.75)

plot_heatmap(oRG_DF3, 'oRG', 0.75)

plot_heatmap(OC_DF3, 'OPC', 0.75)

plot_heatmap(MG_DF3, 'MG', 1)

