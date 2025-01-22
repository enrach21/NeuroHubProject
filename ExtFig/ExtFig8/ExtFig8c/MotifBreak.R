# Date: 12.16.24
# Goal: Evaluate SNP

#### Follow this tutorial ####
# http://www.bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html

#### Read in packages ####
library(dplyr)
library(data.table)
library(motifbreakR)
library(BSgenome)


#### Read in SNPs from rsID: ####

# Read in snp data base
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# read in genome
library(BSgenome.Hsapiens.UCSC.hg38)  
# Read in snps

setwd('/shen/shenlabstore3/ijones1/GKM_explain_test/MG/AD_scores_update')

# obtain list of 192 snipes
snps_ref <- read.table('MG.AD.Sig.Ref.txt', header=0)$X
snps_alt <- read.table('MG.AD.Sig.Alt.txt', header=0)$X
snps <- c(snps_ref, snps_alt)

# Read into motifbreakR
snps.mb <- snps.from.rsid(rsid = 'rs636317',
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)
# Look at variants

snps.mb

# Get the number of variants (180), 12 variants were lost
length(unique(snps.mb$SNP_id ))


#### Find Broken Motifs ####

library(MotifDb)

### Here we can see which organisms are availible under which sources
### in MotifDb
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource)

# leveraged the MotifList introduced by MotifDb to include an additional set of motifs that have been gathered to this package:
data(motifbreakR_motif)
motifbreakR_motif


results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

#### Filter the restuls ####

### Filter Results for strong effects
results_v1 <- results[results$effect == 'strong',]

results_v2 <- calculatePvalue(results_v1, granularity = 1e-6)


RNA <- unique(read.csv('//shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', header=T)[,c('MG','hgnc_symbol')])
colnames(RNA) <- c('MG','Transcription.factor')

### read in hocomoco
hocomoco <- read.csv('~/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv', sep='\t')
dim(hocomoco)
head(hocomoco, n=5)

hocomoco_filt <- hocomoco[,c('Model','Transcription.factor')]
head(hocomoco_filt, n=5)

hocomoco_filt_RNA <- inner_join(hocomoco_filt, RNA)
dim(hocomoco_filt)
dim(hocomoco_filt_RNA )

# Remove if oRG and vRG lowly expressed
hocomoco_filt_RNA2 <- hocomoco_filt_RNA[hocomoco_filt_RNA$MG > 1,]


# Keep only TF greater than 5 RPKM
results_v3 <- results_v2[results_v2$providerName %in% hocomoco_filt_RNA2$Model]


# Write out results
fwrite(data.frame(results_v3@elementMetadata), 'rs636317_SNP_Output.12.16.24.csv')

### Plot 
pdf(file = "rs636317.TF.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4)
plotMB(results = results_v3, rsid = "rs636317", effect = "strong", altAllele = "T")
dev.off()