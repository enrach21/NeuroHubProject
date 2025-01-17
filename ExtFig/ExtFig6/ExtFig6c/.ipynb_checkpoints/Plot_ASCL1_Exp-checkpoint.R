library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(data.table)
library(readxl)
library(ggsignif)
library(ggpubr)

# Read in Motifs
df <- read_excel('ASCL1_HUMAN.H11MO.0.A_ASCL1_HUMAN.H11MO.0.A_overview.xlsx') 


vRG_motif <- df[df$vRG_bound==1 & df$oRG_vRG_log2fc <= -1 ,]
vRG_motif <- vRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','vRG_score','oRG_vRG_log2fc')]
dim(vRG_motif)

#### Bash ####
# PLAC
# bedtools intersect -a vRG_ASCL1_motif.bed -b /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.V2.bedpe -wa -wb | cut -f 1-5,22 > vRG_ASCL1_motif_PLAC.bed

# TSS
# bedtools intersect -a vRG_ASCL1_motif.bed -b /shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed -wa -wb | cut -f 1-5,12 > vRG_ASCL1_motif_TSS.bed

# Get distal Targets
Target <- read.table('vRG_ASCL1_motif_PLAC.bed', fill = T)
colnames(Target) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')
Target <- separate_rows(Target,'hgnc_symbol', sep='\\|')
dim(Target)
head(Target)

# Get TSS Targets
TSS <- read.table('vRG_ASCL1_motif_TSS.bed', fill = T)

colnames(TSS) <- c('chr','start','end','Score','Log2FC','hgnc_symbol')
TSS <- separate_rows(TSS,'hgnc_symbol', sep='\\|')
dim(TSS)
head(TSS)


TF_DF <- rbind(Target,TSS)
head(TF_DF)
dim(TF_DF)


#### Read in RNA

RNA <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', row.names=1)
RNA <- RNA[,c('vRG','oRG','hgnc_symbol')]
head(RNA)

# Filter
RNA_filt <- RNA[RNA$hgnc_symbol %in% TF_DF$hgnc_symbol,]

# Remove targets without expression
RNA_filt <- unique(RNA_filt[RNA_filt$hgnc_symbol != '',])

RNA_filt <- RNA_filt[RNA_filt$vRG > 1 | RNA_filt$oRG > 1, ]

dim(RNA_filt )
head(RNA_filt)
tail(RNA_filt)

# Get DEGS

DEG <- read.csv('../RG.DeSEQ2.DAR.12.1.23.csv')
DEG <- DEG[,c('hgnc_symbol','log2FoldChange','cat')]
head(DEG)

RNA_filt2 <- left_join(RNA_filt, DEG)
dim(RNA_filt2 )
head(RNA_filt2 )
table(RNA_filt2$cat)

RNA_filt2[RNA_filt2$cat == 'vRG',]

summary(log2(RNA_filt$vRG + 1))
summary(log2(RNA_filt$oRG + 1))
t.test(log2(RNA_filt$vRG + 1), log2(RNA_filt$oRG + 1), paired = TRUE, alternative = "two.sided")

df <- melt(RNA_filt)
head(df)
tail(df)
ggplot(df, aes(x=as.factor(variable), y=log2(value +1)))+ 
  geom_violin(aes(fill=variable)) +
  geom_boxplot(width=0.2) +
 geom_signif(comparisons = list(c('vRG', 'oRG')), 
             map_signif_level=TRUE,
              y_position = c(11),
              color = 'black',
             test='t.test',
            test.args=list(alternative = "two.sided", paired=TRUE)) +
  theme_classic() + 
  theme( # plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=10),
        legend.position = "none")  +
    xlab('') +
     ylim (0,13) +
    scale_fill_manual(values=c("#2166AC",'#B2182B')) 
ggsave('ASLC1.Target.Expression.pdf', width=1.25, height = 2)