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

# Read in PLAC-seq output for oRG first

df <- read.table('oRG_LHX2_motif_PLAC.bed', fill=T)
dim(df)
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

dim(df_Target)
head(df_Target )