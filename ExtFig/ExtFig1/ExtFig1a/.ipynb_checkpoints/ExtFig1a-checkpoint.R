
# Load libraries
library(edgeR)
library("DESeq2")
library("tximport")
library("ggplot2")
library("pheatmap")
library(dplyr)
library(tidyr)
library(biomaRt)
library(reshape2)
library("limma")
library("apeglm")
library(genefilter)
library(ggplot2)
library(ggrepel)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

my_colors <- c('#93AACF','#789678','#E3803B','#6B4E30')
names(my_colors) <- c('MG','OPC','oRG','vRG')
my_colors

#### All samples ####
file <- 'All.replicates.txt'
name <- 'All.replicates'

dir <- '.' # location of sample.txt
samples <- read.table(file.path(dir,file), header=TRUE)
samples

#### Updated ####

# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM
# Useful site

dir = '/wynton/group/shen/NeuroHub.7.23.21/RNA-seq/STAR_RSEM_9_7_23/genes.results/' # location of RSEM output
files <- file.path(dir, paste0(samples$samples, ".genes.results"))
names(files) <- samples$samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)


#### Testing changing to using replicate # of batch effect
txi.rsem$length[txi.rsem$length == 0] <- 1 # Make it so genes length cannot be zero
ddsTxi <- DESeqDataSetFromTximport(txi.rsem,
                                    colData = samples,
                                   design = ~date + condition)


# Perform pre-filtering
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]


# Perfrom differential expression
dds <- DESeq(dds)

rld <- rlog(dds, blind=FALSE)

rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)
mm <- model.matrix(~condition, colData(rld))
mat <- limma::removeBatchEffect(mat, batch=rld$date, design=mm)
assay(rld) <- mat
head(assay(rld), 3)



sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(rld$condition,'_R', rld$replicate)
# rownames(sampleDistMatrix) <- paste(vsd$samples, 'GW', vsd$GW, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(c('#3B4CC0', '#8db0fe','#dddddd','#f49a7b', '#b40426')) )(255)
p <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

save_pheatmap_pdf(p, paste0(name,'.sampleDistMatrix.pdf'), width=4.5, height=4)