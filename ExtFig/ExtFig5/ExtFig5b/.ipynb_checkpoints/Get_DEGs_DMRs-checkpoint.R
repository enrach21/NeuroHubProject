# Read in packages
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(ggsignif)
library(data.table)

#### Read in RNA
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


# Get vRG Genes
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


#### Read in PLAC-seq

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'



vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

# Read in 
vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)

vRG_DF2 <- separate_rows(data=vRG_DF, TargetGene,TargetENSG,sep="\\|")

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)

oRG_DF2 <- separate_rows(data=oRG_DF, TargetGene,TargetENSG,sep="\\|")


#### vRG
df <-read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/DMRs/DMR_vRGvsoRG.txt', header=T)
df$chr <- paste0('chr',df$chr)
table(df$dir)

vRG_DMR <- df[df$dir == 'hypo',]
oRG_DMR <- df[df$dir == 'hyper',]

# Filter for DEGs from PLAC contacts

vRG.DEG <- vRG_DF2[vRG_DF2$TargetGene %in% vRG.genes,]

vRG_TSS_DMR <- unique(bt.intersect( vRG_DMR, vRG_TSS,wa=T, wb=T)[,c(1:3,15)])
colnames(vRG_TSS_DMR) <- c('chr','start','end','TargetGene')
table(vRG_TSS_DMR$TargetGene)
vRG_TSS_DMR

vRG_DMR_PLAC <- unique(bt.intersect( vRG_DMR, vRG.DEG,wa=T, wb=T)[,c(1:3,25)])
dim(vRG_DMR_PLAC) 
colnames(vRG_DMR_PLAC) <- c('chr','start','end','TargetGene')
table(vRG_DMR_PLAC$TargetGene)
vRG_DMR_PLAC

vRG_final <- rbind(vRG_TSS_DMR, vRG_D<R_PLAC)

length(table(vRG_final$TargetGene))
table(vRG_final$TargetGene)

write.table(data.frame(table(vRG_final$TargetGene)),'vRG.DMRs_linked_DEGs.txt')

write.table(unique(vRG_final$TargetGene),'vRG.DMR.genes.txt', quote = F, row.names = F, col.names = F)

#### oRG
# Filter for DEGs from PLAC contacts

oRG.DEG <- oRG_DF2[oRG_DF2$TargetGene %in% oRG.genes,]

oRG_TSS_DMR <- unique(bt.intersect( oRG_DMR, oRG_TSS,wa=T, wb=T)[,c(1:3,15)])
colnames(oRG_TSS_DMR) <- c('chr','start','end','TargetGene')
table(oRG_TSS_DMR$TargetGene)
oRG_TSS_DMR

oRG_DAR_PLAC <- unique(bt.intersect( oRG_DMR, oRG.DEG,wa=T, wb=T)[,c(1:3,25)])
dim(oRG_DMR_PLAC) 
colnames(oRG_DMR_PLAC) <- c('chr','start','end','TargetGene')
table(oRG_DMR_PLAC$TargetGene)
oRG_DMR_PLAC

oRG_final <- rbind(oRG_TSS_DMR, oRG_DMR_PLAC)

length(table(oRG_final$TargetGene))
table(oRG_final$TargetGene)

write.table(data.frame(table(oRG_final$TargetGene)),'oRG.DMRs_linked_DEGs.txt')

write.table(unique(oRG_final$TargetGene),'oRG.DMR.genes.txt', quote = F, row.names = F, col.names = F)