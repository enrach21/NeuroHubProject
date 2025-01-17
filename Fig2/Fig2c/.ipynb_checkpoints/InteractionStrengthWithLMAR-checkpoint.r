


#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(ggsignif)


#### Functions ####

# Read in unique ATAC-seq peaks
ReadATAC <- function(dir, file){
    df<-read.table(paste0(dir,file))
    colnames(df) <- c('chr','start','end')
    return(df)
}


# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

OC <- 'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

MG <- 'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

# Read in 
OC_DF <- read.table(paste0(DIR,OC), fill=T, header=T)

MG_DF <- read.table(paste0(DIR,MG), fill=T, header=T)

vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)

# cCRE Overlapping

# Location of called peaks
cCRE_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/'

# Read in Data
MG_cCRE <- ReadATAC(cCRE_dir,'MG.Both.peak.bed')
nrow(MG_cCRE)
OPC_cCRE  <- ReadATAC(cCRE_dir,'OPC.Both.peak.bed')
nrow(OPC_cCRE )
oRG_cCRE  <- ReadATAC(cCRE_dir,'oRG.Both.peak.bed')
nrow(oRG_cCRE )
vRG_cCRE <- ReadATAC(cCRE_dir,'vRG.Both.peak.bed')
nrow(vRG_cCRE)

# Read in Data
MG_ATAC <- ReadATAC(cCRE_dir,'MG.ATAC.peak.bed')
nrow(MG_ATAC)
OPC_ATAC  <- ReadATAC(cCRE_dir,'OPC.ATAC.peak.bed')
nrow(OPC_ATAC )
oRG_ATAC  <- ReadATAC(cCRE_dir,'oRG.ATAC.peak.bed')
nrow(oRG_ATAC )
vRG_ATAC <- ReadATAC(cCRE_dir,'vRG.ATAC.peak.bed')
nrow(vRG_ATAC)

# Read in Data
MG_MR <- ReadATAC(cCRE_dir,'MG.MR.peak.bed')
nrow(MG_MR)
OPC_MR  <- ReadATAC(cCRE_dir,'OPC.MR.peak.bed')
nrow(OPC_MR )
oRG_MR  <- ReadATAC(cCRE_dir,'oRG.MR.peak.bed')
nrow(oRG_MR )
vRG_MR <- ReadATAC(cCRE_dir,'vRG.MR.peak.bed')
nrow(vRG_MR)


classify <- function(PLAC,cCRE, ATAC, LMR) {
    
    # Recover only the distal interactions
    PLAC_XOR <- PLAC[PLAC$type=='XOR',]
    print(paste0('The total number of XOR interactions is ',dim(PLAC_XOR)[1] ))
    
    # Overlap XOR interactions with ATAC-seq bed file
    PLAC_XOR_cCRE<- bt.intersect(a=PLAC_XOR, b=cCRE, wa=TRUE, C=T )
    # Add column names to new data.frame
    colnames(PLAC_XOR_cCRE) <- c(colnames(PLAC_XOR), 'cCRE')
    # Confirm the length is the same as before
    print(paste0('The total number of XOR interactions is ',dim(PLAC_XOR_cCRE)[1] ))
    
    # Overlap XOR interactions with ATAC-seq bed file
    PLAC_XOR_ATAC<- bt.intersect(a=PLAC_XOR_cCRE, b=ATAC, wa=TRUE, C=T )
    # Add column names to new data.frame
    colnames(PLAC_XOR_ATAC) <- c(colnames(PLAC_XOR_cCRE), 'ATAC')
    # Confirm the length is the same as before
    print(paste0('The total number of XOR interactions is ',dim(PLAC_XOR_ATAC)[1] ))
    
    # Overlap XOR interactions with ATAC-seq bed file
    PLAC_XOR_LMR<- bt.intersect(a=PLAC_XOR_ATAC, b=LMR, wa=TRUE, C=T )
    # Add column names to new data.frame
    colnames(PLAC_XOR_LMR) <- c(colnames(PLAC_XOR_ATAC), 'LMR')
    # Confirm the length is the same as before
    print(paste0('The total number of XOR interactions is ',dim(PLAC_XOR_LMR)[1] ))

    # Get Cluster
    cluster <- PLAC_XOR_LMR[,c('ClusterLabel','ClusterSize','ClusterType','ClusterNegLog10P','TargetGene',
                                       'cCRE','ATAC','LMR')]
    head(cluster)
    
    cluster2 <- cluster %>% group_by(ClusterLabel,ClusterSize,ClusterType,ClusterNegLog10P) %>% 
      summarise(cCRE=sum(cCRE),
                ATAC=sum(ATAC),
                LMR=sum(LMR),
                .groups = 'drop')
    dim(cluster)
    head(cluster2)
    dim(cluster2)          

    df <- PLAC_XOR_LMR
    df$lab <- 'Other'
    df$lab[df$cCRE > 0  ] <- 'cCRE'

    df$lab <- factor(df$lab, levels=c('cCRE','Other'))
    
    
    df2 <- cluster2
    df2$lab <- 'Other'
    df2$lab[df2$cCRE > 0 ] <- 'cCRE'
    df2$lab <- factor(df2$lab, levels=c('cCRE','Other'))
    
    return(list(df,df2))
}




#### Obtain results ####

OC_res <- classify(OC_DF,OPC_cCRE,OPC_ATAC, OPC_MR )
MG_res <- classify(MG_DF,MG_cCRE,MG_ATAC,MG_MR)
oRG_res <- classify(oRG_DF,oRG_cCRE, oRG_ATAC, oRG_MR)
vRG_res <- classify(vRG_DF,vRG_cCRE, vRG_ATAC, vRG_MR)

# check counts for XOR contacts
table(OC_res[[1]]$lab)
table(MG_res[[1]]$lab)
table(oRG_res[[1]]$lab)
table(vRG_res[[1]]$lab)

#### Plot percentage interactions Strengths
df <- (oRG_res[[1]])
table(df$lab)
tapply(log2(df$count/df$expected), df$lab, summary)
df$lab <- factor(df$lab, levels = c('cCRE','Other'))
options(repr.plot.width=2, repr.plot.height=3)
ggplot(df, aes(x=lab, y=log2(count/expected)))+ 
  geom_violin(aes(fill=lab)) +
  geom_boxplot(width=0.2) +
  # geom_hline(yintercept=1.919, linetype="dashed", color = "gray15") +

  geom_signif(comparisons = list(c('cCRE', 'Other')), 
             map_signif_level=TRUE,
              y_position = c(7),
              color = 'black',
             test='t.test',
            test.args=list(alternative = "two.sided", paired=FALSE)) +
  theme_classic() + 
  theme( # plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=10),
        legend.position = "none")  +
    xlab('') +
    ylim(c(0,9)) +
        scale_x_discrete(labels=c('cCRE\nn=26,913',
                              'Other\nn=85,267'))    + 
    scale_fill_manual(values=c("#B2182B",'grey85')) 
ggsave('oRG.Interaction_score.pdf', width=2, height = 3)

df <- (MG_res[[1]])
table(df$lab)
tapply(log2(df$count/df$expected), df$lab, summary)
df$lab <- factor(df$lab, levels =  c('cCRE','Other'))
options(repr.plot.width=2, repr.plot.height=3)
ggplot(df, aes(x=lab, y=log2(count/expected)))+ 
  geom_violin(aes(fill=lab)) +
  geom_boxplot(width=0.2) +
  # geom_hline(yintercept=1.868, linetype="dashed", color = "gray15") +

  geom_signif(comparisons = list(c('cCRE', 'Other')), 
             map_signif_level=TRUE,
              y_position = c(7),
              color = 'black',
             test='t.test',
            test.args=list(alternative = "two.sided", paired=FALSE)) +
  theme_classic() + 
  theme( # plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=10),
        legend.position = "none")  +
    xlab('') +
    ylim(c(0,9)) +
        scale_x_discrete(labels=c('cCRE\nn=27,209',
                              'Closed\nn=75,841'))      + 
    scale_fill_manual(values=c("#762A83",'grey85')) 
ggsave('MG.Interaction_score.pdf', width=2, height = 3)

options(repr.plot.width=2, repr.plot.height=3)
df <- (OC_res[[1]])
table(df$lab)
tapply(log2(df$count/df$expected), df$lab, summary)
df$lab <- factor(df$lab, levels = c('cCRE','Other'))
table(df$lab)
ggplot(df, aes(x=lab, y=log2(count/expected)))+ 
  geom_violin(aes(fill=lab)) +
  geom_boxplot(width=0.2) +
  # geom_hline(yintercept=1.949, linetype="dashed", color = "red") +

  geom_signif(comparisons = list(c('cCRE', 'Other')), 
             map_signif_level=TRUE,
              y_position = c(7),
              color = 'black',
             test='t.test',
            test.args=list(alternative = "two.sided", paired=FALSE)) +
  theme_classic() + 
  theme( # plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=10),
        legend.position = "none")  +
    xlab('') +
    ylim(c(0,9)) +
        scale_x_discrete(labels=c('cCRE\nn=23,296',
                              'Other\nn=84,964'))       + 
    scale_fill_manual(values=c("#228833",'grey85')) 
ggsave('OC.Interaction_score.pdf', width=2, height = 3)

