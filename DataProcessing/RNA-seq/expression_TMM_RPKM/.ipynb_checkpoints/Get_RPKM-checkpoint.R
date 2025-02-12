#### Packages
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library("tximport")
library(mixOmics)
library(statmod)
library(biomaRt)
library(dplyr)
library(pheatmap)
library(tidyr)

#### Read in Data

# Make sample table
dir <- '.' # location of sample.txt
samples <- read.table(file.path(dir,"All.replicates.txt"), header=TRUE)
samples

# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM
# Useful site

dir = '/wynton/group/shen/NeuroHub.7.23.21/RNA-seq/STAR_RSEM_9_7_23/genes.results/' # location of RSEM output
files <- file.path(dir, paste0(samples$samples, ".genes.results"))
names(files) <- samples$samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

y <- DGEList(txi.rsem$counts, group=samples$condition, genes=data.frame(Length=txi.rsem$length))
y$samples

# Calculate normalized factors
y <- calcNormFactors(y)
y$sample

# Gety$counts$PKM
norm_rpkm_count <- rpkm(y, log=F, gene.length = y$genes, normalized.lib.size=T)
dim(norm_rpkm_count)
head(norm_rpkm_count)

# Remove NA
norm_rpkm_count <- na.omit(norm_rpkm_count)
norm_rpkm_count <- norm_rpkm_count[!is.infinite(rowSums(norm_rpkm_count)),]
dim(norm_rpkm_count)
head(norm_rpkm_count)

# Get the mean expression
norm_rpkm_count<-data.frame(norm_rpkm_count)
oRG_mean <- (norm_rpkm_count$Length.IJ345 + norm_rpkm_count$Length.JJ027 + norm_rpkm_count$Length.JJ101 +  norm_rpkm_count$Length.JJ184)/4
vRG_mean <- (norm_rpkm_count$Length.IJ344 + norm_rpkm_count$Length.JJ028 + norm_rpkm_count$Length.JJ102 +  norm_rpkm_count$Length.JJ185)/4
MG_mean <- (norm_rpkm_count$Length.IJ346 + norm_rpkm_count$Length.JJ029 + norm_rpkm_count$Length.JJ096)/3
Oligo_mean <- (norm_rpkm_count$Length.IJ347 + norm_rpkm_count$Length.JJ031 + norm_rpkm_count$Length.JJ097)/3

# Make new data frame with the average expression
norm_rpkm_count_avg <- data.frame(row.names = rownames(norm_rpkm_count), 
                                   'oRG'=oRG_mean,'vRG'=vRG_mean,
                                  'MG'=MG_mean,'Oligo'=Oligo_mean)
dim(norm_rpkm_count_avg )
head(norm_rpkm_count_avg)

# Combine the two data.frames
df <- data.frame(norm_rpkm_count_avg,norm_rpkm_count)
dim(df)
head(df)

write.csv(df,'AVG_TMM_RPKM_NoLog_exp_9.20.23.csv' )

#### Convert ENSG to geneID

df$egids <- rownames(df)

df <- df %>% separate(egids, c('Tmp1', 'Tmp2'), sep='\\.')

egids <- df$Tmp1

df$ensembl_gene_id <- egids

df <- df[,!colnames(df) %in% c('Tmp1', 'Tmp2')]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

translated <- getBM(filters="ensembl_gene_id", 
                    attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                    values=egids, mart=ensembl)

df_new <- unique(inner_join(df, translated[,1:2], multiple = "all"))
head(df_new, n=10)
dim(df_new )
dim(df)

write.csv(df_new,'AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv' )