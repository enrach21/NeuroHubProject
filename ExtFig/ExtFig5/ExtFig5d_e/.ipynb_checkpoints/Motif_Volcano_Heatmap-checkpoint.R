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

# Read in RNA-seq
df_rna <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/AVG_TMM_RPKM_NoLog_exp_geneID_4.13.23.txt', row.names = 1)
df_rna2 <-  df_rna[,c('oRG','vRG','MG','Oligo','hgnc_symbol')]
colnames(df_rna2)[5] <- 'Transcription.factor'
head(df_rna2)

# Join with bin detect
df3 <- left_join(df2, df_rna2, multiple = "all")
dim(df3)

# Remove any row with a non-expressed TF
df4 <- df3[!is.na(df3$vRG),]
dim(df3)
dim(df4)
head(df4)

# Down cell type
Down <- 'Oligo'
Up <- 'MG'
# Cut off for TF expression
RPKM <- 10

# Remove low expressed TF
df5 <- df4[!(df4[,10]< 0 & (df4[,Down] < RPKM)), ]
dim(df5)

df6 <- df5[!(df5[,10] > 0 & (df5[,Up] < RPKM)), ]
dim(df6)

# Add color for significants
df6$cat <- 'NOT_SIG'
df6$cat[df6[,12] == 'True' & df6[,10] > 0] <- Up
df6$cat[df6[,12] == 'True' & df6[,10] < 0] <- Down
table(df6$cat)

# Add TF names for significatns elements
df$label <- NA
df6$label[df6[,12] == 'True' & df6[,10] > 0] <- df6$Transcription.factor[df6[,12] == 'True' & df6[,10] > 0]
df6$label[df6[,12] == 'True' & df6[,10] < 0] <- df6$Transcription.factor[df6[,12] == 'True' & df6[,10] < 0]

# Final DF
df7 <- unique(df6[,c(1:13,18,19)])
dim(df6)
dim(df7)# oRG and vRG
mycolors <- c("#762A83","#1B7837","black")
names(mycolors) <- c(Up, Down, 'NOT_SIG')


options(repr.plot.width=2, repr.plot.height=2)
require("ggrepel")
p <- ggplot(df7) +
        geom_point(aes(x=MG_OPC_change, y=-log10(MG_OPC_pvalue), colour=cat), size=0.5) +
        geom_label_repel(aes(x=MG_OPC_change, y=-log10(MG_OPC_pvalue),label=label),
                         min.segment.length = unit(0, 'lines'), size=2, max.overlaps = 14) +
        ggtitle(paste0(Down, ' / ',Up, ' motif binding')) +
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
ggsave('MG.OPC.Vol.pdf',width = 2, height = 2)

# Plot Heatmap
df8 <- unique(df6[df6[,12] == 'True' ,])
df8 <- df8[order(df8[,10]),]
dim(df8)
head(df8)

TF_df <- matrix(df8[,10])
rownames(TF_df) <- df8$Transcription.factor
colnames(TF_df) <- 'Motif_Dif'
dim(TF_df )

RNA_df <- matrix(log2(df8[,c(Up)] / df8[,c(Down)]))
rownames(RNA_df) <- df8$Transcription.factor
colnames(RNA_df) <- 'RPKM'
RNA_df

p <- cor.test(TF_df, RNA_df, method = 'pearson')
p

library(circlize)

# colfunc <- colorRampPalette(c("#93AACF","White","#789678"))

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#1B7837",'White',"#762A83"))

p1 <- Heatmap(TF_df, cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col=col_fun ,  
              heatmap_legend_param = list(
                title = 'Motif Differences',
               legend_direction = "vertical", 
               legend_height = unit(0.5, "cm"),
              labels_gp = gpar( font = 4) ), 
              row_names_gp = gpar(fontsize = 6))

# colfunc <- colorRampPalette(c("blue","red"))
col_fun = colorRamp2(c(-5, 0, 5), c("#1B7837",'White',"#762A83"))
p2 <- Heatmap(as.matrix(RNA_df), cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col= (col_fun) ,
          heatmap_legend_param = list(
                title = 'RPKM (MG/OPC)',
               legend_direction = "vertical", 
               legend_height = unit(0.5, "cm"),
              labels_gp = gpar( font = 6)))

options(repr.plot.width=4, repr.plot.height=5)
pdf('MG.OC.Motif.Exp.pdf', height = 5, width = 4)
print(p1 + p2)
dev.off()
              
              
              