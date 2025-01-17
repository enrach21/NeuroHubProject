


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

#### Plot percentage of LMAR interactions
my_colors <- rev(c('#762A83','#1B7837','#B2182B','#2166AC'))
names(my_colors) <- rev(c('MG','OPC','oRG','vRG'))
my_colors

output <- data.frame()
list1 <-  c(OC_res,MG_res,oRG_res,vRG_res)
list2 <-  c('OPC','OPC','MG','MG','oRG','oRG','vRG','vRG')
for (x in c(1,3,5,7)) {
    print(x)
    df <- list1[x]
    name <- list2[x]
    tab <- table(df[[1]]$lab)
    print(tab)
    temp_df <- data.frame(tab/sum(tab))
    temp_df$counts <- tab
    temp_df$cell <- name
    print(temp_df)
    output <- rbind(output, temp_df)
}



output <- output[output$Var1 == 'cCRE',]
output$color <- paste0(output$cell,'_',output$Var1)
output$color <- factor(output$color, levels = rev(c(
        'MG_cCRE',
    'OPC_cCRE',
    'oRG_cCRE',
    'vRG_cCRE')))
output$cell <- factor(output$cell, levels = rev(c(
        'MG',
    'OPC',
    'oRG',
    'vRG')))
output


options(repr.plot.width=1.75, repr.plot.height=1.75)
ggplot(data=output, aes(x=cell, y=Freq * 100, fill=color, color=Var1)) +
  geom_bar(stat="identity", width=0.8)  +
        xlab(NULL) + 
        ylab('Percentage of XOR 
interations with cCRE') +
         theme_classic() +
        theme(
            legend.position = "none",
       axis.text.y = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=8)) +
    scale_fill_manual(values=c('#2166AC',
                              '#B2182B',
                              '#1B7837',
                    '#762A83'),guide = "none") + 
    scale_color_manual(values=c("black", "black", "black","black")) +
    ylim(0,30)

ggsave('cCRE.Pecentage.pdf',width=1.75, height=1.75)