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
library('stringr')

#### Read in SNPs
snps <- read.csv('TableBrowser-var_list_scz3GWAS_2more_sig_deepGWASsig_hg38.csv')
head(snps)

dim(snps)

snps <- snps[snps$X.chrom %in% paste0('chr',c(1:22)),]
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

#### ATAC

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

#### Methyl

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

#### PLAC
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


df_final <- snps[(snps$vRG_ATAC > 0 | snps$oRG_ATAC > 0 | snps$OC_ATAC > 0  | snps$MG_ATAC > 0 ),]
dim(df_final)


# Get multiple alternative alleles
df_final <- separate_rows(df_final, alts, sep = ',')
# Remove empty rows
df_final <- df_final[df_final$alts != '',]


# change size to 200bp
df_final$chromStart <- df_final$chromStart-99
df_final$chromEnd <- df_final$chromStart+200
summary(df_final$chromEnd-df_final$chromStart)

write.csv(df_final,'Scz.Filt.allATAC.SNPs.200bp.5.1.24.csv')

# Get alternate allele
alt <- df_final$alts
alt
length(alt)



#200bp
write.table(df_final, 'Scz.Filt.allATAC.SNPs.200bp.5.1.24.bed',col.names=F, row.names = F,sep='\t', quote=F)

bedtools getfasta -fi /shen/shenlabstore3/shared/reference_genome/hg38/rsem_star/GRCh38.p10.genome.fa -bed ./Scz.Filt.allATAC.SNPs.200bp.5.1.24.bed > Scz.Filt.allATAC.SNPs.200bp.5.1.24.fa

bedtools getfasta -fi /shen/shenlabstore3/shared/reference_genome/hg38/rsem_star/GRCh38.p10.genome.fa -bed ./Scz.Filt.allATAC.SNPs.200bp.5.1.24.bed > Scz.Filt.allATAC.SNPs.200bp.5.1.24.fa

write.fasta(fasta,make.unique(names(fasta)), 'Scz.Filt.allATAC.SNPs.200bp.5.1.24.fa')
fasta<- read.fasta('Scz.Filt.allATAC.SNPs.200bp.5.1.24.fa')


#### Replace pos 100 with alt
for (x in 1:length(fasta)){
    fasta[[x]][100:(100 + nchar(alt[x]) - 1)] <- str_split(alt[x], pattern = "")[[1]]
}

write.fasta(fasta,make.unique(names(fasta)), 'Scz.Filt.allATAC.ALT.SNPs.200bp.5.1.24.fa')

#### Make shuffled regions
df <- read.csv('Scz.Filt.allATAC.SNPs.200bp.5.1.24.csv', row.names=1)
head(df)
dim(df)

# Get ref allele
ref <- df$ref
ref
length(ref)

# Get alt allele
alt <- df$alts
alt
length(alt)

shuf <- read.fasta('Scz.Shuffled.SNPs.200bp.5.1.24.fa')
length(shuf)

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

write.fasta(fasta_alt,make.unique(names(fasta_alt)), 'Scz_Alt_Shuffled.SNPs.200bp.5.1.24.fa')
write.fasta(fasta_ref,make.unique(names(fasta_ref)), 'Scz_Ref_Shuffled.SNPs.200bp.5.1.24.fa')


#### Make on line one

# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' Scz_Alt_Shuffled.SNPs.200bp.5.1.24.fa > Scz_Alt_Shuffled.SNPs.200bp.5.1.24.V2.fa

# awk '!/^>/ { printf "%s", $0; n = "\n" } 
# /^>/ { print n $0; n = "" }
# END { printf "%s", n }
# ' Scz_Ref_Shuffled.SNPs.200bp.5.1.24.fa > Scz_Ref_Shuffled.SNPs.200bp.5.1.24.V2.fa
