#### Read in Package ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library(gprofiler2)
library("ggrepel")
library(forcats)

# Get HAR regions
# Read in unique ATAC-seq peaks
ReadCRE <- function(file){
    df<-read.table(file, header=F)[,1:3]
    colnames(df) <- c('chr','start','end')
    return(df)
}

HAR <- '../GSE180714_HARs.bed'
HAR2 <- '../GSE180714_HARs.tsv'
HAR_DF <- ReadCRE(HAR)
HAR_DF <- HAR_DF[-1,]
dim(HAR_DF)
head(HAR_DF)

HAR2_DF <- read.table(HAR2, header=T)
head(HAR2_DF)

# Read in GKM signifcant hits
GKM <- read.csv('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/HARs_update/HAR.8.12.24.csv', row.names = 1)
GKM <- GKM[GKM$ISM_pval <= 0.05 & GKM$delta_p_val <= 0.05 & GKM$explain_pval <= 0.05,]
dim(GKM)
head(GKM)

file <- GKM
HAR <- bt.intersect(a= HAR2_DF, b=file, wa=T, u=T)
dim(unique(HAR))
head(HAR)

# Read in PLAC

# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'


oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)

oRG_DF <- separate_rows(oRG_DF, c('TargetENSG', 'TargetGene'), sep='\\|')

# oRG_DF <- oRG_DF[oRG_DF$type == 'XOR',]
head(oRG_DF)

# PLAC_overlap with HARs
HAR_Target <- bt.intersect(a=oRG_DF, b=HAR, wa=T, u=T)$V17
length(HAR_Target)
head(HAR_Target)

# Overlap with TSS
TSS <- unique(bt.intersect( a='/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed',b=HAR , u=T)$V7)
length(TSS)
head(TSS)

# Combine Target and TSS to get total genes 
total <- unique(c(HAR_Target,TSS))
length(total)
total


# Filter for exprssion
RNA <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', header=T)[,c(2,21)]

RNA_filt <- RNA[(RNA$hgnc_symbol %in% total) & (RNA$hgnc_symbol != '') & (RNA$oRG > 1),]

length(RNA_filt$hgnc_symbol)

RNA_filt

write.table(RNA_filt$hgnc_symbol,'GKM_Sig_Target_genes.txt',quote = F, col.names = F, row.names = F)

Gene_1 <- RNA_filt$hgnc_symbol

g.term3 <- gost(query = list("HAR"=Gene_1), organism = "hsapiens", ordered_query = F, significant = T, correction_method = "g_SCS", sources = "GO:BP", evcodes = TRUE)


test2 <- test[test$term_id %in% c('GO:0008283',
'GO:0006357',
'GO:0048731',
'GO:0007417'),]
test2

fwrite(test, 'oRG.GKM.HAR.Targets.go.csv')

# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis
options(repr.plot.width=5, repr.plot.height=1.5)
ggplot(test2, aes(x=query, y= fct_reorder(term_name, -log10(p_value)))) +
    geom_point(aes(size=precision, color = -log10(p_value))) +
    theme_classic() +
    scale_colour_gradient(lim=c(0, 4), high='#B2182B', low='white',na.value = '#B2182B') +
    labs(title="GO Biological Process Enrichment",x="Timepoint", y = "") +
    labs(size="Precision") +
    theme(legend.title = element_text(color = "black", size = 8),
            legend.text = element_text(color = "black"),
            legend.key.size = unit(.5, "cm"),
            axis.text.x = element_text(face="plain", color="black", 
                           size=6, angle=45, hjust=1),
            axis.text.y = element_text(face="plain", color="black", 
                           size=6, angle=0),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            strip.text.x = element_blank(),
             plot.title = element_text(face="plain", color="black", 
                              size = 10, 
                              hjust = .5), legend.box = "horizontal") 
ggsave('oRG.GKM.HAR.Targets.go.pdf',width = 5, height = 1.5)