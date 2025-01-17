# Read in packages
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(ggsignif)

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

# This time consider both XOR and AND interaction
classify2 <- function(PLAC,cCRE) {
    
    # Recover only the distal interactions
    PLAC <- PLAC
    
    # Overlap XOR interactions with ATAC-seq bed file
    PLAC_cCRE<- bt.intersect(a=PLAC, b=cCRE, wa=TRUE, C=T )
    # Add column names to new data.frame
    colnames(PLAC_cCRE) <- c(colnames(PLAC), 'cCRE')
    # Confirm the length is the same as before
    print(paste0('The total number of interactions is ',dim(PLAC_cCRE)[1] ))

    # Get Cluster
    cluster <- PLAC_cCRE[,c('ClusterLabel','ClusterSize','ClusterType','ClusterNegLog10P','TargetGene',
                                       'cCRE')]
    head(cluster)
    
    cluster2 <- cluster %>% group_by(ClusterLabel,ClusterSize,ClusterType,ClusterNegLog10P) %>% 
      summarise(cCRE=sum(cCRE),
                .groups = 'drop')
    dim(cluster)
    head(cluster2)
    dim(cluster2)          

    df <- PLAC_cCRE
    df$lab <- 'Other'
    df$lab[df$cCRE > 0] <- 'cCRE'

    df$lab <- factor(df$lab, levels=c('cCRE','Other'))
    
    
    df2 <- cluster2
    df2$lab <- 'Other'
    df2$lab[df2$cCRE > 0 ] <- 'cCRE'
    df2$lab <- factor(df2$lab, levels=c('cCRE','Other'))
    
    return(list(df,df2))
}

OC_res <- classify2(OC_DF,OPC_cCRE)
MG_res <- classify2(MG_DF,MG_cCRE)
oRG_res <- classify2(oRG_DF,oRG_cCRE)
vRG_res <- classify2(vRG_DF,vRG_cCRE)

output2 <- data.frame()
list1 <-  c(OC_res,MG_res,oRG_res,vRG_res)
list2 <-  c('OPC','OPC','MG','MG','oRG','oRG','vRG','vRG')
for (x in c(2,4,6,8)) {
    print(x)
    df <- list1[x]
    name <- list2[x]
    tab <- table(df[[1]]$lab)
    print(tab)
    temp_df <- data.frame(tab/sum(tab))
    temp_df$counts <- tab
    temp_df$cell <- name
    print(temp_df)
    output2 <- rbind(output2, temp_df)
}

output2$color <- paste0(output2$cell,'_',output2$Var1)
output2$color <- factor(output2$color, levels = c(
        'MG_cCRE','MG_Other',
    'OPC_cCRE','OPC_Other',
    'oRG_cCRE','oRG_Other',
    'vRG_cCRE','vRG_Other'))
output2

options(repr.plot.width=5, repr.plot.height=1.75)
ggplot(data=output2, aes(x=cell, y=Freq, fill=color, color=Var1)) +
  geom_bar(stat="identity")  +
        xlab(NULL) + 
        ylab(NULL) +
         theme_classic() +
        theme(
            legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=30, hjust=1)) +
   scale_fill_manual(values=c('#762A83', 'gray50',
                              '#1B7837','gray50',
                              '#B2182B','gray50',
                    '#2166AC','gray50'),guide = "none") + 
    scale_color_manual(values=c("black", "black", "black","black"))  +
    coord_flip() +
    geom_text(aes(label=counts),color="white",size=4,position=position_stack(vjust=0.5))



ggsave('cCRE_Cluster.Overlap.pdf',width=5, height=1.75)



