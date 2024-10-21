# Plot differential TFs between vRG and oRG


#### Packages ####
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

#### Read in output from BINDectect ####

df <- read.table('diffTFBS_MG_OPC_bindetect_results.txt', h=T)

#### Read in HOMOCO TFs ####
hocomoco <- read.csv('~/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv', sep='\t')
hocomoco_filt <- hocomoco[,c('Model','Transcription.factor')]
colnames(hocomoco_filt) <- c('motif_id','Transcription.factor')

#### Merge both Dataframe Together ####
df2 <- left_join(df, hocomoco_filt, multiple = "all")
dim(df2)


#### Read in Transcription ####
df_rna <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/AVG_TMM_RPKM_NoLog_exp_geneID_4.13.23.txt', row.names = 1)
df_rna2 <-  df_rna[,c('oRG','vRG','MG','Oligo','hgnc_symbol')]
colnames(df_rna2)[5] <- 'Transcription.factor'


#### Join expression with Motif data ####
df3 <- left_join(df2, df_rna2, multiple = "all")
dim(df3)


#### Remove any row with a non-expressed TF ####
df4 <- df3[!is.na(df3$Oligo),]
dim(df3)
dim(df4)


#### Only retain TFs that are highly expressed ####
# Down cell type
Down <- 'Oligo'
Up <- 'MG'
# Cut off for TF expression
RPKM <- 10

# Remove low expressed TF
df5 <- df4[!(df4[,'MG_OPC_change']< 0 & (df4[,Down] < RPKM)), ]
dim(df5)

df6 <- df5[!(df5[,'MG_OPC_change'] > 0 & (df5[,Up] < RPKM)), ]
dim(df6)

# Add color for significants
df6$cat <- 'NOT_SIG'
df6$cat[df6[,12] == 'True' & df6[,10] > 0] <- Up
df6$cat[df6[,12] == 'True' & df6[,10] < 0] <- Down
table(df6$cat)


# Add TF names for significatns elements
df$label <- NA
df6$label[df6[,12] == 'True' & df6[,10] > 0] <- df6$Transcription.factor[df6[,12] == 'True' & df6[,10] > 0]
df6$label[df6$Transcription.factor %in% c('SOX9','IRF3','SP1','SOX10','ELK1') ] <- df6$Transcription.factor[df6$Transcription.factor %in% c('SOX9','IRF3','SP1','SOX10','ELK1')]

write.table(df6,'Diff.MG_OPC_binding.txt')


# Final DF
df7 <- unique(df6[,c(1:13,18,19)])
dim(df6)
dim(df7)

# oRG and vRG
mycolors <- rev(c("#1B7837","#762A83","black"))
names(mycolors) <- c(Up, Down, 'NOT_SIG')



options(repr.plot.width=2.5, repr.plot.height=2.5)
require("ggrepel")
p <- ggplot(df7) +
        geom_point(aes(x=MG_OPC_change, y=-log10(MG_OPC_pvalue), colour=cat), size=0.5) +
        geom_label_repel(aes(x=MG_OPC_change, y=-log10(MG_OPC_pvalue),label=label),min.segment.length = unit(0, 'lines'), size=1, max.overlaps = 30) +
        ggtitle(paste0(Up, ' / ',Down, ' motif binding')) +
        xlab("log2 fold change") + 
        ylab("-log10 p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = 8, hjust = 1),
              axis.title = element_text(size = rel(1.25))) +
        theme_classic(base_size = 6) + scale_color_manual(values=mycolors) + 
        theme(plot.title = element_text(size=8, hjust = 0.5)) +
        theme(legend.title = element_blank(), legend.key.size = unit(0.3, 'cm')) + 
        theme(legend.text=element_text(size=6)) + 
        guides(colour = guide_legend(override.aes = list(size=0.5)))
p  + guides(fill="none")+ 
    theme(legend.position = "none")
ggsave('OPC.MG.Vol.pdf',width = 2.5, height = 2.5)