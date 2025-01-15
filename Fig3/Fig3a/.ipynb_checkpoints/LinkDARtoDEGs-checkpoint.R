# Read in packages
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(ggsignif)
library(data.table)
library(gprofiler2)
library("ggrepel")
library(forcats)

# Read in RNA
RNA <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/DeSeq2.RG.Genes/RG.DeSEQ2.9.13.23.csv')
dim(RNA)

# Get significant Genes
RNA.sig <- RNA[RNA$cat == 'oRG' | RNA$cat == 'vRG',]
dim(RNA.sig)
table(RNA.sig$cat)

# Get genes with only gene.id
RNA.sig2 <- RNA.sig[RNA.sig$hgnc_symbol != '',]
dim(RNA.sig2)

# Read in expression
RPKM <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Marker_expression/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv')[,c('oRG','vRG','hgnc_symbol')]
head(RPKM)

RNA.Final <- inner_join(RNA.sig2,RPKM)

# Get Genes
vRG.genes <- RNA.Final[RNA.Final$cat == 'vRG','hgnc_symbol']
length(vRG.genes )

oRG.genes <- RNA.Final[RNA.Final$cat == 'oRG','hgnc_symbol']
length(oRG.genes )

# Read in TSS  
Tss <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed'
TssDF <- read.table(Tss, header=FALSE)
colnames(TssDF) <- c('chr','start','end', 'strand','ENST','ENSG','gene')
print('Dimensions of TSS file...')
dim(TssDF)
print('View TSS file...')
head(TssDF)

# get vRG TSS
vRG_TSS <- TssDF[TssDF$gene %in% vRG.genes,]
dim(vRG_TSS)
head(vRG_TSS)

oRG_TSS <- TssDF[TssDF$gene %in% oRG.genes,]
dim(oRG_TSS)
head(oRG_TSS)

#### Read in PLAC-seq ####

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'



vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

# Read in 
vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)

vRG_DF2 <- separate_rows(data=vRG_DF, TargetGene,TargetENSG,sep="\\|")

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)

oRG_DF2 <- separate_rows(data=oRG_DF, TargetGene,TargetENSG,sep="\\|")

#### vRG ####

# Read in DARs
vRG_DAR <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG/vRG.DAR.bed')
head(vRG_DAR)
dim(vRG_DAR)

# Filter for DEGs from PLAC contacts

vRG.DEG <- vRG_DF2[vRG_DF2$TargetGene %in% vRG.genes,]
head(vRG.DEG)

vRG_TSS_DAR <- unique(bt.intersect( vRG_DAR, vRG_TSS,wa=T, wb=T)[,c(1:3,11)])
colnames(vRG_TSS_DAR) <- c('chr','start','end','TargetGene')
table(vRG_TSS_DAR$TargetGene)
vRG_TSS_DAR

vRG_DAR_PLAC <- unique(bt.intersect( vRG_DAR, vRG.DEG,wa=T, wb=T)[,c(1:3,21)])
dim(vRG_DAR_PLAC) 
colnames(vRG_DAR_PLAC) <- c('chr','start','end','TargetGene')
table(vRG_DAR_PLAC$TargetGene)
vRG_DAR_PLAC

vRG_final <- rbind(vRG_TSS_DAR, vRG_DAR_PLAC)

length(table(vRG_final$TargetGene))
table(vRG_final$TargetGene)

# Write number of overlap per gene
write.table(data.frame(table(vRG_final$TargetGene)),'vRG.DARs_linked_DEGs.txt')


# Write out gene name of genes overlappign with DARs
write.table(unique(vRG_final$TargetGene),'vRG.DAR.genes.txt', quote = F, row.names = F, col.names = F)

#### oRG ####

# Read in DAR update 11.14.24
oRG_DAR <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG/oRG.DAR.bed')
head(oRG_DAR)
dim(oRG_DAR)

# Filter for DEGs from PLAC contacts

oRG.DEG <- oRG_DF2[oRG_DF2$TargetGene %in% oRG.genes,]

oRG_TSS_DAR <- unique(bt.intersect( oRG_DAR, oRG_TSS,wa=T, wb=T)[,c(1:3,11)])
colnames(oRG_TSS_DAR) <- c('chr','start','end','TargetGene')
table(oRG_TSS_DAR$TargetGene)
oRG_TSS_DAR

# Filter for DEGs from PLAC contacts

oRG.DEG <- oRG_DF2[oRG_DF2$TargetGene %in% oRG.genes,]

oRG_TSS_DAR <- unique(bt.intersect( oRG_DAR, oRG_TSS,wa=T, wb=T)[,c(1:3,11)])
colnames(oRG_TSS_DAR) <- c('chr','start','end','TargetGene')
table(oRG_TSS_DAR$TargetGene)
oRG_TSS_DAR

oRG_DAR_PLAC <- unique(bt.intersect( oRG_DAR, oRG.DEG,wa=T, wb=T)[,c(1:3,21)])
dim(oRG_DAR_PLAC) 
colnames(oRG_DAR_PLAC) <- c('chr','start','end','TargetGene')
table(oRG_DAR_PLAC$TargetGene)
oRG_DAR_PLAC

oRG_final <- rbind(oRG_TSS_DAR, oRG_DAR_PLAC)

length(table(oRG_final$TargetGene))
table(oRG_final$TargetGene)

# Write out overlap per differentially expressed genes
write.table(data.frame(table(oRG_final$TargetGene)),'oRG.DARs_linked_DEGs.txt')

# Write out gene names themselves
write.table(unique(oRG_final$TargetGene),'oRG.DAR.genes.txt', quote = F, row.names = F, col.names = F)