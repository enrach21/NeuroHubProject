# Code to make figure 1C

#### Packages ####
library('dplyr')
library('pheatmap')
library(RColorBrewer)
library(ggplot2)

# color scheme
my_colors <- c('#2166AC','#B2182B','#1B7837','#762A83')
names(my_colors) <- c('vRG','oRG','OPC','MG')
my_colors

args <- commandArgs(trailingOnly=TRUE)

# Read in variable
RNA <- args[1]
print(RNA)
MARKER <- args[2]
print(MARKER)

# Read in expression
df <- read.csv(RNA, row.names = 1)

# Read in Markers
Chosen_Markers <- read.csv(MARKER)
colnames(Chosen_Markers) <- c('Cell Type','hgnc_symbol')
Chosen_Markers


# Join with normalized expression
Chosen_expression <- inner_join(df, Chosen_Markers)
row.names(Chosen_expression) <- make.names(Chosen_expression$hgnc_symbol, unique=T)
head(Chosen_expression )

# Reorder based on Markers
Chosen_expression$'Cell Type' <- factor(Chosen_expression$'Cell Type', levels=rev(c('MG','OPC','oRG','vRG')))
Chosen_expression <- Chosen_expression[order(Chosen_expression$'Cell Type'),]

# Prepare plot for pheatmap
df <- Chosen_expression
markerDF <- data.frame(df$'Cell Type')
order <- unique(markerDF$df.marker)
rownames(markerDF) <- rownames(df)
colnames(markerDF) <- 'Cell Type'
markerDF$'Cell Type'<- as.factor(markerDF$'Cell Type')

# Change the color of the Markers
annotation_color <- list('Cell Type'=my_colors)
annotation_color

# Get color
RdBu <- rev(brewer.pal(9,"RdBu"))
length(RdBu)

p <- pheatmap(t(scale(t(log2(Chosen_expression[,c(2,1,4,3)]+1)))), cluster_rows=F, 
            annotation_row=markerDF, annotation_names_row = FALSE, 
              annotation_colors = annotation_color, fontsize=8, main='Z-Score (RPKM)',
             color = colorRampPalette(RdBu)(100),
             cluster_cols=F)


pdf( file='Chosen_markers_log2_RPKM_Zscore.pdf', width=2.5,height = 2.5)
p
dev.off()
