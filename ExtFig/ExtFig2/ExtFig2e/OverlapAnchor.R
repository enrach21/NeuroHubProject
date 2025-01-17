
#### Load Packages ####
library('dplyr')
library('ggplot2')
library('bedtoolsr')
library(reshape2) 


#### Read in LMAR cCREs ####
vRG <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed')
oRG <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed')
OPC <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.Both.peak.bed')
MG <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.Both.peak.bed')


#### Read in Anchors ####
anchor <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Anchor/Merged/Anchor.Merged.9.7.23.final.bed')

head(anchor)
dim(anchor)

bin <- read.table('/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/BinnedGenome/hg38.2kb.bed')
head(bin)

Anchor_2kb <- bt.intersect(a=bin, b=anchor, wa=T, u=T)
dim(Anchor_2kb)
head(Anchor_2kb)

#### Get counts ####

# Example of how the counts were obtained and then entered into csv
dim(bt.intersect(vRG, Anchor_2kb, wa=T, u=T))
dim(bt.intersect(vRG, Anchor_2kb, v=T))


#### LMAR ####

df <- read.csv('cCRE_Overlap - cCRE.csv')
df

options(repr.plot.width=3, repr.plot.height=1.25)
ggplot(data=df, aes(x=Cell, y=Value, fill=Color, color=Variable)) +
  geom_bar(stat="identity")  +
        xlab(NULL) + 
        ylab(NULL) +
         theme_classic() +
        theme(
            legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=30, hjust=1)) +
   scale_fill_manual(values=c( '#762A83','#BBBBBB',
                              '#1B7837','#BBBBBB',
                             '#B2182B','#BBBBBB',
                    "#2166AC",'#BBBBBB'),guide = "none") + 
    scale_color_manual(values=c("black", "black"))  +
    coord_flip()  +
    geom_text(aes(label=Value), vjust=0.5, hjust=1.23, color=c(rep("black",8)), size=3)


ggsave('LMAR.H3K4me3.Overlap.pdf',width=3, height=1.25)


#### AR ####

df <- read.csv('cCRE_Overlap - ATAC.csv')
df

options(repr.plot.width=3, repr.plot.height=1.25)
ggplot(data=df, aes(x=Cell, y=Value, fill=Color, color=Variable)) +
  geom_bar(stat="identity")  +
        xlab(NULL) + 
        ylab(NULL) +
         theme_classic() +
        theme(
            legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=30, hjust=1)) +
   scale_fill_manual(values=c( '#762A83','#BBBBBB',
                              '#1B7837','#BBBBBB',
                             '#B2182B','#BBBBBB',
                    "#2166AC",'#BBBBBB'),guide = "none") + 
    scale_color_manual(values=c("black", "black"))  +
    coord_flip()  
    # geom_text(aes(label=Value), vjust=0.5, hjust=1.23, color=c(rep("black",8)), size=3)


ggsave('AR.H3K4me3.Overlap.pdf',width=3, height=1.25)



#### LMR ####

df <- read.csv('cCRE_Overlap - LMR.csv')
df

options(repr.plot.width=3, repr.plot.height=1.5)
ggplot(data=df, aes(x=Cell, y=Value, fill=Color, color=Variable)) +
  geom_bar(stat="identity")  +
        xlab(NULL) + 
        ylab(NULL) +
         theme_classic() +
        theme(
            legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=30, hjust=1)) +
   scale_fill_manual(values=c( '#762A83','#BBBBBB',
                              '#1B7837','#BBBBBB',
                             '#B2182B','#BBBBBB',
                    "#2166AC",'#BBBBBB'),guide = "none") + 
    scale_color_manual(values=c("black", "black"))  +
    coord_flip()  
    # geom_text(aes(label=Value), vjust=0.5, hjust=1.23, color=c(rep("black",8)), size=3)


ggsave('LMR.H3K4me3.Overlap.pdf',width=3, height=1.25)