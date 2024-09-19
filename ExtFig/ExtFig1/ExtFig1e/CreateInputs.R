# Make sample table
dir <- '/shen/shenlabstore3/ijones1/NeuroHub/RNA-seq/' # location of sample.txt
samples <- read.table(file.path(dir,"samples4.txt"), header=TRUE)
samples

# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM
# Useful site

dir = '/wynton/group/shen/NeuroHub.7.23.21/RNA-seq/GeneResults/' # location of RSEM output
files <- file.path(dir, paste0(samples$samples, ".genes.results"))
names(files) <- samples$samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

# Make input for CIBERSORT
# Make sure to read txi.rsem before changing it
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
row.names(txi.rsem$counts) <- make.names( unlist(strsplit(as.character(row.names(txi.rsem$counts)), split="\\."))[(1:length(row.names(txi.rsem$counts)))*2-1],unique=T)
CIBER_BULK_REF <-as.matrix(df)
write.table(CIBER_BULK_REF, '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Bisque/Intermediate_Files/CIBER_BULK_Output_NH_10_14_22.tsv',quote=F, sep='\t')


#### Single-cell input ####
#### Read in Packages ####
library(BisqueRNA)
library(Biobase)
library(dplyr)
library("tximport")
library(annotables)

#### Get single cell data ####

# Read in raw count matrix from Tom's 2017 2nd trimester paper
sc.counts <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Bisque/Intermediate_Files/TOM_N_COUNT_DATA/count_matrix.tsv')

# Get cell name id's and put into a vector and then make into a df with column called Cell
sample.ids <- colnames(sc.counts)
sample_df <- data.frame('Cell'=sample.ids)

# Read in meta data for each cell
meta <- read.table('//shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Bisque/Intermediate_Files/TOM_N_COUNT_DATA/meta.tsv', sep='\t', header=T)

# Get the labels for each cell cluster
cell.type.labels <- meta$WGCNAcluster

# Merge the cells identified above with their meta data
meta_sub <- inner_join(sample_df, meta)

# Combine all eN clusters into one
meta_sub$WGCNAcluster[(meta_sub$WGCNAcluster == 'EN-PFC1' | 
                      meta_sub$WGCNAcluster == 'EN-PFC2' |
                      meta_sub$WGCNAcluster == 'EN-PFC3' |
                      meta_sub$WGCNAcluster == 'nEN-early1' |
                      meta_sub$WGCNAcluster == 'nEN-early2' |
                      meta_sub$WGCNAcluster == 'EN-late')
                      ] <- 'EN'

# Combine all iN clusters into one
meta_sub$WGCNAcluster[(meta_sub$WGCNAcluster == 'IN-CTX-CGE1' | 
        meta_sub$WGCNAcluster == 'IN-CTX-CGE2' |
        meta_sub$WGCNAcluster == 'IN-CTX-MGE1' |
        meta_sub$WGCNAcluster == 'IN-CTX-MGE2' |
        meta_sub$WGCNAcluster == 'IN-STR' |
        meta_sub$WGCNAcluster == 'InIN1' |
        meta_sub$WGCNAcluster == 'InIN2' |
        meta_sub$WGCNAcluster == 'InIN3' |
        meta_sub$WGCNAcluster == 'InIN4' |                     
        meta_sub$WGCNAcluster == 'InIN5')] <- 'IN'

# Combine all iPC clusters into one
meta_sub$WGCNAcluster[(meta_sub$WGCNAcluster == 'IPC-nEN1' | 
        meta_sub$WGCNAcluster == 'IPC-nEN2' |
        meta_sub$WGCNAcluster == 'IPC-nEN3')] <- 'IPC'

# Select the cluster of single cells wanted to be used in the final comparison
meta_sub2 <- meta_sub[meta_sub$WGCNAcluster == 'IPC' | meta_sub$WGCNAcluster == 'Astrocyte' | meta_sub$WGCNAcluster == 'EN' | meta_sub$WGCNAcluster == 'OPC' | meta_sub$WGCNAcluster == 'Microglia' | meta_sub$WGCNAcluster == 'tRG' | meta_sub$WGCNAcluster == 'oRG' | meta_sub$WGCNAcluster == 'vRG',]

# Remove expression from celltypes of non interest
sc.counts.sort <- sc.counts[,colnames(sc.counts) %in% meta_sub2$Cell]
sc.counts.matrix <- as.matrix(sc.counts.sort)

# Write out single cell reference sample file
CIBER_SC_REF <- sc.counts.matrix
colnames(CIBER_SC_REF) <- meta_sub2$WGCNAcluster
head(CIBER_SC_REF)
write.table(CIBER_SC_REF, '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Bisque/Intermediate_Files/CIBER_SC_Output_Tom.tsv',quote=F, sep='\t')