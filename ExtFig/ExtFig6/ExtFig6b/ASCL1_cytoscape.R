library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(data.table)
library(readxl)

# Read in Motifs
df <- read_excel('ASCL1_HUMAN.H11MO.0.A_ASCL1_HUMAN.H11MO.0.A_overview.xlsx') 


vRG_motif <- df[df$vRG_bound==1 & df$oRG_vRG_log2fc <= -1 ,]
vRG_motif <- vRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','vRG_score','oRG_vRG_log2fc')]
dim(vRG_motif)
write.table(vRG_motif, 'vRG_ASCL1_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')

#### IN bash ####
# PLAC
# bedtools intersect -a vRG_ASCL1_motif.bed -b /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.V2.bedpe -wa -wb | cut -f 1-5,22 > vRG_ASCL1_motif_PLAC.bed

# TSS
# bedtools intersect -a vRG_ASCL1_motif.bed -b /shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed -wa -wb | cut -f 1-5,12 > vRG_ASCL1_motif_TSS.bed
df <- read.table('vRG_ASCL1_motif_PLAC.bed', fill=T)
colnames(df) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')

# Separate bins that had multiple gene targets
df <- separate_rows(df,'hgnc_symbol', sep='\\|')

# Remove rows without target genes
df <- unique(df[df$hgnc_symbol!='',])
dim(df)
head(df)

#### Read in RNA

RNA <- read.csv('../RG.DeSEQ2.DAR.12.1.23.csv')
RNA <- RNA[,c('hgnc_symbol','log2FoldChange','cat')]


#### Combine RNA with PLAC Targets

df_Target <- left_join(df,RNA)

# Remove targets without expression
df_Target <- unique(df_Target[!is.na(df_Target$cat),])
table(df_Target$cat)

df_Target$type <- 'distal'

dim(df_Target)
head(df_Target)

# Get distal Targets
TSS <- read.table('vRG_ASCL1_motif_TSS.bed', fill = T)

colnames(TSS) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')
TSS <- separate_rows(TSS,'hgnc_symbol', sep='\\|')
dim(TSS)
head(TSS)

#### Combine RNA with PLAC Targets

df_TSS <- left_join(df_TSS,RNA)

# Remove targets without expression
df_TSS <- unique(df_TSS[!is.na(df_TSS$cat),])
table(df_TSS$cat)

df_TSS$type <- 'TSS'

dim(df_TSS)
length(unique(df_TSS$hgnc_symbol))
head(df_TSS)

# Combine the two results
df_Final <- rbind(df_Target, df_TSS)
df_Final$TF <- 'ASCL1'

# Add RPKM values

RNA <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', row.names=1)
RNA <- RNA[,c('vRG','oRG','hgnc_symbol')]
RNA[RNA$hgnc_symbol == 'ASCL1',]
head(RNA)

df_Final <- left_join(df_Final, RNA)
head(df_Final)

df_Final <- df_Final[df_Final$cat != 'NOT_SIG',]
dim(df_Final)

df_Final$log2_RPKM <- log2(df_Final$vRG/df_Final$oRG) 


write.csv(df_Final,'ASCL1.cytoscape.csv', row.names = F)