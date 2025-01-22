#### Packages 

#### Read in packages #####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library("UpSetR")
library("ComplexHeatmap")
# library(gprofiler2)
library("ggrepel")
library(biomaRt)
library(dplyr)
library(corrplot)
library('seqinr')
library(stringr)

#### Read in snps
AD_SNPs_DF <- read.csv('41588_2023_1506_MOESM3_ESM.csv')

AD_SNPs <- unique(na.omit(AD_SNPs_DF$rsid))
head(AD_SNPs)
length(AD_SNPs)

# Used snps as input into UCSC table browser to get genomic locations
write.table(AD_SNPs, 'AD_SNPs.txt', quote = F, col.names = F, row.names = F)

# Read in input from table browser
snps <- read.csv('AD.SNPs.1.4.24.TableBrowser.csv')
head(snps)
dim(snps)


ATAC1 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS//vRG.4reps.overlap.optimal_peak.narrowPeak.gz'
ATAC2 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS//oRG.4reps.overlap.optimal_peak.narrowPeak.gz'
ATAC3 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS//OPC.overlap.optimal_peak.narrowPeak.gz'
ATAC4 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS//MG.overlap.optimal_peak.narrowPeak.gz'


MR1 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/vRG_UMRsLMRs.tab'
MR2 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/oRG_UMRsLMRs.tab'
MR3 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/OPC_UMRsLMRs.tab'
MR4 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/MG_UMRsLMRs.tab'

PLAC1 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'
PLAC2 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'
PLAC3 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'
PLAC4 <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'


Anchor <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Anchor/Merged/Anchor.Merged.7.27.22.final.bed'


#### Overlap with ATAC

tmp<- bt.intersect(snps,ATAC1, C=T)
colnames(tmp) <- c(colnames(snps),'vRG_ATAC')
snps <- tmp

tmp<- bt.intersect(snps,ATAC2, C=T)
colnames(tmp) <- c(colnames(snps),'oRG_ATAC')
snps <- tmp

tmp<- bt.intersect(snps,ATAC3, C=T)
colnames(tmp) <- c(colnames(snps),'OC_ATAC')
snps <- tmp

tmp<- bt.intersect(snps,ATAC4, C=T)
colnames(tmp) <- c(colnames(snps),'MG_ATAC')
snps <- tmp

# Overlap with methyl

MRseq <- read.table(MR1,header=T)

tmp<- bt.intersect(snps,MRseq, C=T)
colnames(tmp) <- c(colnames(snps),'vRG_MR')
snps <- tmp

MRseq <- read.table(MR2,header=T)

tmp<- bt.intersect(snps,MRseq, C=T)
colnames(tmp) <- c(colnames(snps),'oRG_MR')
snps <- tmp

MRseq <- read.table(MR3,header=T)

tmp<- bt.intersect(snps,MRseq, C=T)
colnames(tmp) <- c(colnames(snps),'OC_MR')
snps <- tmp

MRseq <- read.table(MR4,header=T)

tmp<- bt.intersect(snps,MRseq, C=T)
colnames(tmp) <- c(colnames(snps),'MG_MR')
snps <- tmp

# Overlap PLAC-seq

PLACseq <- fread(PLAC1)

tmp<- bt.intersect(snps,PLACseq, C=T)
colnames(tmp) <- c(colnames(snps),'vRG_PLAC')
snps <- tmp

PLACseq <- fread(PLAC2)

tmp<- bt.intersect(snps,PLACseq, C=T)
colnames(tmp) <- c(colnames(snps),'oRG_PLAC')
snps <- tmp

PLACseq <- fread(PLAC3)

tmp<- bt.intersect(snps,PLACseq, C=T)
colnames(tmp) <- c(colnames(snps),'OC_PLAC')
snps <- tmp

PLACseq <- fread(PLAC4)

tmp<- bt.intersect(snps,PLACseq, C=T)
colnames(tmp) <- c(colnames(snps),'MG_PLAC')
snps <- tmp


# Overlap with Anchors
tmp<- bt.intersect(snps,Anchor , C=T)
colnames(tmp) <- c(colnames(snps),'Anchor')
snps <- tmp


sum(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 )
sum((snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 ) & snps$Anchor == 0)

filt_snps <-snps[(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 )  & snps$Anchor ==0,]
head(filt_snps)
nrow(filt_snps)

# write.table(snps[(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0),], 'Accessible.RG.scz.8.2.23.txt')
write.csv(snps[(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 ),], 'Accessible.AD.8.1.24.csv')


#### Write out accessible snps

df_final <- snps[(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 ),]
dim(df_final)
df_final <- df_final[df_final$ref != '',]
dim(df_final )
head(df_final)

# Split alternative aleles
# Get multiple alternative alleles
df_final <- separate_rows(df_final, alts, sep = ',')
# Remove empty rows
df_final <- df_final[df_final$alts != '',]

dim(df_final)
head(df_final)

write.csv(df_final, 'Accessible.AD.8.1.24.csv')

write.table(df_final$name, 'AD_SNPs.txt', quote = F, col.names = F, row.names = F)

# change size to 200bp
df_final$chromStart <- df_final$chromStart-99
df_final$chromEnd <- df_final$chromStart+200
summary(df_final$chromEnd-df_final$chromStart)

write.table(df_final$name, 'AD_SNPs.txt', quote = F, col.names = F, row.names = F)

# Get alternate allele
alt <- df_final$alts
alt
length(alt)


write.table(df_final, 'MG.AD.Filt.allATAC.SNPs.200bp.bed',col.names=F, row.names = F,sep='\t', quote=F)

# IN BASH

# bedtools getfasta -fi /shen/shenlabstore3/shared/reference_genome/hg38/rsem_star/GRCh38.p10.genome.fa -bed ./MG.AD.Filt.allATAC.SNPs.200bp.bed > MG.AD.Filt.allATAC.SNPs.200bp.fa

# Read in the fasta
fasta <- read.fasta('MG.AD.Filt.allATAC.SNPs.200bp.fa')

# Make the fasta have unique names
write.fasta(fasta,make.unique(names(fasta)), 'MG.AD.Filt.allATAC.SNPs.200bp.fa')
fasta<- read.fasta('MG.AD.Filt.allATAC.SNPs.200bp.fa')

### Replace variant with Alt at pos 100
for (x in 1:length(fasta)){
    fasta[[x]][100:(100 + nchar(alt[x]) - 1)] <- str_split(alt[x], pattern = "")[[1]]
}

write.fasta(fasta,make.unique(names(fasta)), 'MG.AD.Filt.allATAC.ALT.SNPs.200bp.fa')

#### Make shuffled controls 
#### Run shuffle.200bp.sh

# Get alternate allele
ref <- df_final$ref
ref
length(ref)

# Get alternate allele
alt <- df_final$alts
alt
length(alt)

shuf <- read.fasta('Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa')

fasta_alt <- shuf
fasta_ref <- shuf
num <- 3
for (x in 1:length(ref)){
    for (y in 1:num) {
        z <- ((x-1)*num)+y

        fasta_alt[[z]][100:(100 + nchar(alt[x]) - 1)] <- str_split(alt[x], pattern = "")[[1]]
    
        fasta_ref[[z]][100:(100 + nchar(ref[x]) - 1)] <- str_split(ref[x], pattern = "")[[1]]
    
        
    }
}

write.fasta(fasta_alt,make.unique(names(fasta_alt)), 'Alt_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa')
write.fasta(fasta_ref,make.unique(names(fasta_ref)), 'Ref_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa')

#### Make fasta 1 line in bash
# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' Alt_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa > Alt_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.fa

# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' Ref_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa > Ref_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.fa

# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' MG.AD.Filt.allATAC.SNPs.200bp.fa > MG.AD.Filt.allATAC.SNPs.200bp.V2.fa

# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' MG.AD.Filt.allATAC.ALT.SNPs.200bp.fa > MG.AD.Filt.allATAC.ALT.SNPs.200bp.V2.fa