#### Packages ####
library(data.table)
library(dplyr)
library('preprocessCore')
library('bedtoolsr')
library(pheatmap)
library(ggplot2)
library(ggplotify)
library(ComplexHeatmap)
library(patchwork)
library('tidyr')
library(data.table)
library(circlize)

#### Read in all interactions ####
df <- read.csv('PLAC_RNA_Cor - All_int.csv', row.names=1)


df <- read.csv('PLAC_RNA_Cor - All_int.csv', row.names=1)

options(repr.plot.width=2.5, repr.plot.height=2.5)

col_fun = colorRamp2(c(0,.25, .5), c("#2166AC","#F7F7F7", "#B2182B"))


h <- Heatmap(df, cluster_rows = FALSE,, cluster_columns = FALSE,
               # change Row names 90
        row_title_rot = 0,
        # Change size and color of rows
        row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
       column_title = 'Correlation between PLAC and RNA',
               column_title_gp = gpar(fontsize = 8),
       column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(
        fontsize = c(6,6,6,6)),
       row_names_rot = 0,
        row_names_centered = TRUE,
        row_names_gp = gpar(
        fontsize = c(6,6,6,6)),
        col = col_fun,
        #Change Legend
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
}
       
       )
        # Add title to describe the data)

pdf('All.Int.Cor.PLAC.RNA.7.17.24.pdf', height = 2.5, width = 2.5)
draw(h)
dev.off()


#### LMAR cCREs ####

df <- read.csv('PLAC_RNA_Cor - cCRE_int.csv', row.names=1)

options(repr.plot.width=2.5, repr.plot.height=2.5)

col_fun = colorRamp2(c(0,.25, .5), c("#2166AC","#F7F7F7", "#B2182B"))


h <- Heatmap(df, cluster_rows = FALSE,, cluster_columns = FALSE,
               # change Row names 90
        row_title_rot = 0,
        # Change size and color of rows
        row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
       column_title = 'Correlation between PLAC and RNA',
               column_title_gp = gpar(fontsize = 8),
       column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(
        fontsize = c(6,6,6,6)),
       row_names_rot = 0,
        row_names_centered = TRUE,
        row_names_gp = gpar(
        fontsize = c(6,6,6,6)),
        col = col_fun,
        #Change Legend
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
}
       
       )
        # Add title to describe the data)

pdf('cCRE.Int.Cor.PLAC.RNA.7.17.24.pdf', height = 2.5, width = 2.5)
draw(h)
dev.off()


#### AR ####

df <- read.csv('PLAC_RNA_Cor - ATAC_int.csv', row.names=1)

options(repr.plot.width=2.5, repr.plot.height=2.5)

col_fun = colorRamp2(c(0,.25, .5), c("#2166AC","#F7F7F7", "#B2182B"))


h <- Heatmap(df, cluster_rows = FALSE,, cluster_columns = FALSE,
               # change Row names 90
        row_title_rot = 0,
        # Change size and color of rows
        row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
       column_title = 'Correlation between PLAC and RNA',
               column_title_gp = gpar(fontsize = 8),
       column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(
        fontsize = c(6,6,6,6)),
       row_names_rot = 0,
        row_names_centered = TRUE,
        row_names_gp = gpar(
        fontsize = c(6,6,6,6)),
        col = col_fun,
        #Change Legend
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
}
       
       )
        # Add title to describe the data)

pdf('ATAC.Int.Cor.PLAC.RNA.7.17.24.pdf', height = 2.5, width = 2.5)
draw(h)
dev.off()

#### LMR ####

df <- read.csv('PLAC_RNA_Cor - LMR_int.csv', row.names=1)

options(repr.plot.width=2.5, repr.plot.height=2.5)

col_fun = colorRamp2(c(0,.25, .5), c("#2166AC","#F7F7F7", "#B2182B"))


h <- Heatmap(df, cluster_rows = FALSE,, cluster_columns = FALSE,
               # change Row names 90
        row_title_rot = 0,
        # Change size and color of rows
        row_title_gp = gpar(col = c("#332288","#882255"), fontsize = c(6,6)),
       column_title = 'Correlation between PLAC and RNA',
               column_title_gp = gpar(fontsize = 8),
       column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(
        fontsize = c(6,6,6,6)),
       row_names_rot = 0,
        row_names_centered = TRUE,
        row_names_gp = gpar(
        fontsize = c(6,6,6,6)),
        col = col_fun,
        #Change Legend
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
}
       
       )
        # Add title to describe the data)

pdf('LMR.Int.Cor.PLAC.RNA.7.17.24.pdf', height = 2.5, width = 2.5)
draw(h)
dev.off()