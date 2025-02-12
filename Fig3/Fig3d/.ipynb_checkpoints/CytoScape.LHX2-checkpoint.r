#### Read in Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(data.table)
library(readxl)

#### Get motif's of interest ####
# Read in Motifs
df <- read_excel('LHX2_HUMAN.H11MO.0.A_LHX2_HUMAN.H11MO.0.A_overview.xlsx') 

# Filter for motifs bound in oRG and >1 log2FC
oRG_motif <- df[df$oRG_bound==1 & df$oRG_vRG_log2fc >=1 ,]
oRG_motif <- oRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','oRG_score','oRG_vRG_log2fc')]

write.table(oRG_motif, 'oRG_LHX2_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')

#### Read in output from bash file ####

#### Read in PLAC-seq output ####

df <- read.table('oRG_LHX2_motif_PLAC.bed', fill=T)
dim(df)
colnames(df) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')

# Separate bins that had multiple gene targets
df <- separate_rows(df,'hgnc_symbol', sep='\\|')

# Remove rows without target genes
df <- unique(df[df$hgnc_symbol!='',])
length(unique(df$hgnc_symbol))
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
length(unique(df_Target$hgnc_symbol))
head(df_Target)

#### Read in PLAC-seq output ####

df_TSS <- read.table('oRG_LHX2_motif_TSS.bed', fill=T)
dim(df)
colnames(df_TSS) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')
head(df_TSS)

#### Combine RNA with PLAC Targets

df_TSS <- left_join(df_TSS,RNA)

# Remove targets without expression
df_TSS <- unique(df_TSS[!is.na(df_TSS$cat),])
table(df_TSS$cat)

df_TSS$type <- 'TSS'

dim(df_TSS)
length(unique(df_TSS$hgnc_symbol))
head(df_TSS)

df_Final <- rbind(df_Target,df_TSS)
dim(df_Final)
head(df_Final)
table(unique(df_Final[,c('hgnc_symbol', 'cat')])$cat)

# Add RPKM values

RNA <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', row.names=1)
RNA <- RNA[,c('vRG','oRG','hgnc_symbol')]
RNA[RNA$hgnc_symbol == 'LHX2',]
head(RNA)

df_Final <- left_join(df_Final, RNA)
head(df_Final)

df_Final <- df_Final[df_Final$cat != 'NOT_SIG',]
dim(df_Final)

df_Final$log2_RPKM <- log2(df_Final$vRG/df_Final$oRG) 


write.csv(df_Final,'LHX2.cytoscape.csv', row.names = F)