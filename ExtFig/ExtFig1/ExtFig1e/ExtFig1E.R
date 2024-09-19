# Code to make figure 1C

#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library("UpSetR")
library("ComplexHeatmap")
library(gprofiler2)
library("ggrepel")

library(RColorBrewer)
library(viridis)

# color scheme
my_colors <- c('#2166AC','#B2182B','#1B7837','#762A83')
names(my_colors) <- c('vRG','oRG','OPC','MG')
my_colors

args <- commandArgs(trailingOnly=TRUE)

# Read in variable
CB.Output <- args[1]
print(CB.Output)
RNA.lab <- args[2]
print(RNA.lab)

# Read in Data 
df <- read.csv(CB.Output)
df
df <- df[c(1:14),]
df

Lab <- read.csv(RNA.lab)
Lab

# combine label and output
df_lab <- inner_join(df, Lab)
df_lab


# Reorder based on cell types
df_lab$Cell_Type <- factor(df_lab$Cell_Type,
                           levels = rev(c('MG','OPC','oRG','vRG')))
df_lab_sort <- df_lab[order(df_lab$Cell_Type),]
head(df_lab_sort)

# Change row names to sample ID's
rownames(df_lab_sort) <- df_lab_sort$Mixture

# Make annotation DF
anno_df <- df_lab_sort[,c('GW','Cell_Type')]

ann_colors = list(
    # Donor_ID = c('ARK-2021-002'="goldenrod1",'ARK-2021-005'="goldenrod2",'ARK-2021-008'="goldenrod3",'ARK-2021-009'="darkgoldenrod3",'ARK-2021-010'="goldenrod4",'ARK-2021-014'="darkgoldenrod4"),
    GW = c(GW15='#FEE391', GW18='#FB9A29',GW22='#CC4C02',GW23='#662506'),
    Cell_Type = c(MG = '#762A83', OPC = '#1B7837', oRG = '#B2182B', vRG = '#2166AC')
)

# Make matrix for pheatmap
plot_df <- df_lab_sort[,c('vRG','tRG','oRG','Astrocyte','OPC','Microglia','IPC','EN')]
head(plot_df)

# Get color
RdBu <- rev(brewer.pal(9,"RdBu"))
length(RdBu)

my.colors <- colorRampPalette(RdBu)(100)
length(my.colors)

mat_breaks <- seq(0,0.5, length.out=100)
length(mat_breaks)

p <- pheatmap(as.matrix(plot_df), row_title = NULL,
 cluster_rows=TRUE, cluster_cols=FALSE, annotation_row = anno_df,  show_rownames = F, border_color = NA, breaks=mat_breaks, color = my.colors,
              display_numbers = TRUE, annotation_colors = ann_colors, fontsize=8)


pdf("CiberSortHeatMap.pdf", width=4.5,height = 4)
p
dev.off()