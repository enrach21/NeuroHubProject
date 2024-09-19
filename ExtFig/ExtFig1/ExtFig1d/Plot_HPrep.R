library(data.table)
library(plyr)
library(ggplot2)
library(tibble)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)


# Read in output from HPRep

cor <- read.table('/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/HPREP/2023_10mil_usable/HPRep_output/stage3/NeuroHub.10million.9_8_23.Norm.results.txt')
# Takes the mean correlation between all the chromosomes
cor$meanVal <- apply(cor[3:24], 1, mean, na.rm = TRUE)
#
cor_df <- cor[,c(1,2,25)]
cor_df2 <- cor_df[,c(2,1,3)]
colnames(cor_df2) <- c('V1','V2','meanVal')
cor_final <- rbind(cor_df, cor_df2 )


# undo melt
unmelt_cor <- dcast(cor_final,V1~V2, value.var = 'meanVal')
# Replace NA values of same sampel cor with 1
unmelt_cor[is.na(unmelt_cor)] <- 1


# replace rownames with the first column
rownames(unmelt_cor) <- unmelt_cor[,1]
# Remove the first column
unmelt_cor <- unmelt_cor[,2:13]


rownames(unmelt_cor) <- c('OC_R2','oRG_R3','MG_R2','vRG_R2','oRG_R2','OC_R3','vRG_R3', 'OC_R1','MG_R1','oRG_R1','vRG_R1','MG_R3')
colnames(unmelt_cor) <- c('OC_R2','oRG_R3','MG_R2','vRG_R2','oRG_R2','OC_R3','vRG_R3', 'OC_R1','MG_R1','oRG_R1','vRG_R1','MG_R3')


#quantile_breaks
quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(unmelt_cor), n = 20)
my.colors <- colorRampPalette(colors = c('#3B4CC0', '#8db0fe','#dddddd','#f49a7b', '#b40426'))(length(mat_breaks))

pdf('HPrep.10mil.9.12.23.pdf', height=4.5, width=4.5)
pheatmap(unmelt_cor, display_numbers=F, breaks= mat_breaks, color = my.colors,  fontsize = 8)
dev.off()
