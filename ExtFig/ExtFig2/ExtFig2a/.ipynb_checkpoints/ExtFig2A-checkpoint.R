#### Packages ####
library('dplyr')
library('ggplot2')
library('bedtoolsr')
library(reshape2) 


# color scheme
my_colors <- c('#2166AC','#B2182B','#1B7837','#762A83')
names(my_colors) <- c('vRG','oRG','OPC','MG')
my_colors

#### Functions ####

# Read in unique ATAC-seq peaks
ReadATAC <- function(dir, file){
    df<-read.table(paste0(dir,file))
    colnames(df) <- c('chr','start','end')
    return(df)
}

# Read in unique methylation peaks
ReadMRs <- function(dir, file){
    df<-read.table(paste0(dir,file), header=T)
    return(df)
}

#### ATAC-seq ####
# ATAC-seq Overlapping

# Location of called peaks
ATAC_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/'

# Read in Data
MG_ATAC <- ReadATAC(ATAC_dir,'MG.overlap.optimal_peak.unique.bed')
nrow(MG_ATAC)
OPC_ATAC <- ReadATAC(ATAC_dir,'OPC.overlap.optimal_peak.unique.bed')
nrow(OPC_ATAC)
oRG_ATAC <- ReadATAC(ATAC_dir,'oRG.4reps.overlap.optimal_peak.unique.bed')
nrow(oRG_ATAC)
vRG_ATAC <- ReadATAC(ATAC_dir,'vRG.4reps.overlap.optimal_peak.unique.bed')
nrow(vRG_ATAC)

#### LMR ####

# Location of called peaks
MRs_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/'

MG_MRs <- ReadMRs(MRs_dir,'MG_UMRsLMRs.tab')
nrow(MG_MRs)
OPC_MRs <- ReadMRs(MRs_dir,'OPC_UMRsLMRs.tab')
nrow(OPC_MRs)
oRG_MRs <- ReadMRs(MRs_dir,'oRG_UMRsLMRs.tab')
nrow(oRG_MRs)
vRG_MRs <- ReadMRs(MRs_dir,'vRG_UMRsLMRs.tab')
nrow(vRG_MRs)


ATAC <- list(MG_ATAC, OPC_ATAC, oRG_ATAC, vRG_ATAC)
MRs <- list(MG_MRs, OPC_MRs, oRG_MRs,vRG_MRs)
Cell <- list('MG','OPC','oRG','vRG')
df <- data.frame()
for (i in c(1:4)) {
    a <- ATAC[[i]]
    m <- MRs[[i]]
    c <- Cell[i]
    
    x <- nrow(bt.intersect(a,m, wa=T, u=T))
    # print(x)
    y <- nrow(bt.intersect(m,a, wa=T, u=T))
    # print(y)
    z <- nrow(bt.intersect(a,m, v=T))
    # print(z)
    z2 <- nrow(bt.intersect(m,a, v=T))
    # print(z2)
    
    write.table(bt.intersect(a,m, wa=T, u=T), paste0(c,'.Both.peak.bed'), col.names = F, row.names = F, quote = F, sep='\t')
    write.table(bt.intersect(a,m, v=T), paste0(c,'.ATAC.peak.bed'), col.names = F, row.names = F, quote = F, sep='\t')
    write.table(bt.intersect(m,a, v=T), paste0(c,'.MR.peak.bed'), col.names = F, row.names = F, quote = F, sep='\t')
    
    df <- rbind(df, c(max(x,y),z,z2, c))
    colnames(df) <- c('LMARs','ARs','LMRs','Cell')
}
df2 <- melt(df)


# Add a variable for the color
df2$color <- paste0(df2$Cell,'_',df2$variable)
df2$color <- factor(df2$color, levels = c(
        'MG_LMRs','MG_ARs','MG_LMARs',
    'OPC_LMRs','OPC_ARs', 'OPC_LMARs',
    'oRG_LMRs','oRG_ARs','oRG_LMARs',
    'vRG_LMRs','vRG_ARs','vRG_LMARs'))
df2


# Change the order of the variable
df2$variable <- factor(df2$variable, levels = c('LMARs', 'ARs','LMRs'))


options(repr.plot.width=3, repr.plot.height=2)
ggplot(data=df2, aes(x=Cell, y=value, fill=color, color=variable)) +
  geom_bar(stat="identity", position=position_dodge(width=.9), size=.25)  +
        xlab(NULL) + 
        ylab(NULL) +
         theme_classic() +
        theme(
            #legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=30, hjust=1),
              axis.text.y = element_text(face="plain", color="black", 
                           size=8, angle=0, hjust=1),
              legend.text = element_text(size=6),
              legend.title=element_blank(),
            legend.key.width = unit(.4, "cm"),
             legend.key.height = unit(.4, "cm")
        ) +
   scale_fill_manual(values=rev(c("#2166AC",'#4393C3','#92C5DE',
                   '#B2182B','#D6604D','#F4A582',
                    '#1B7837','#5AAE61','#ACD39E',
                   '#762A83','#9970AB','#C2A5CF')),guide = "none") + 
    scale_color_manual(values=c("#262626", "#7F7F7F", "#D9D9D9"))  +
    coord_flip() +
    geom_text(aes(label=value), vjust=0.5, hjust=1.23, color=c(rep("white",8),rep("black",4)),
            position = position_dodge(0.9), size=2)


ggsave('ATAC.MR.Overlap.pdf',width=3, height=2)