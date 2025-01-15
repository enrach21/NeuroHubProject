# Plot heatmap of TF expression and binding difference as well as correlation


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
library(circlize)

#### Read in data.frame ####

df6 <- read.table('../ExtFig8_OCMG/Diff.MG_OPC_binding.txt')

df8 <- unique(df6[df6[,12] == 'True' ,])
df8 <- df8[order(df8[,10]),]
dim(df8)
head(df8)


#### Get motif difference ####

TF_df <- matrix(df8[,10])
rownames(TF_df) <- df8$Transcription.factor
colnames(TF_df) <- 'Motif_Dif'

#### Get Expression Difference ####

# Down cell type
Down <- 'Oligo'
Up <- 'MG'


RNA_df <- matrix(log2(df8[,c(Up)] / df8[,c(Down)]))
rownames(RNA_df) <- df8$Transcription.factor
colnames(RNA_df) <- 'RPKM'

p <- cor.test(TF_df, RNA_df, method = 'pearson')
p


#### Plot motif and expression

options(repr.plot.width=2, repr.plot.height=3)
colfunc <- colorRampPalette(c("#1B7837","White","#762A83"))
p1 <- Heatmap(TF_df, cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col=colfunc(20),  
              heatmap_legend_param = list(
                title = 'Motif Differences',
               legend_direction = "vertical", 
               legend_height = unit(0.5, "cm"),
              labels_gp = gpar( font = 4) ), 
              row_names_gp = gpar(fontsize = 6))

p1

options(repr.plot.width=2, repr.plot.height=3)
col_fun = colorRamp2(c(-1, 0, 1),c("#1B7837","White","#762A83"))


p2 <- Heatmap(as.matrix(RNA_df), cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col= col_fun ,
          heatmap_legend_param = list(
                title = 'RPKM (oRG/vRG)',
               legend_direction = "vertical", 
               legend_height = unit(0.5, "cm"),
              labels_gp = gpar( font = 4) ), 
              row_names_gp = gpar(fontsize = 6))
p2

options(repr.plot.width=2.5, repr.plot.height=4)
print(p1 + p2)

pdf('OPC.MG.Motif.Exp4.pdf', height = 3, width = 2.5)
print(p1 + p2)
dev.off()