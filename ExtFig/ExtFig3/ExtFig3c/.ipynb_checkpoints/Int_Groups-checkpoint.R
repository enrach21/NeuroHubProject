# Packages 
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(ggplot2)
library("UpSetR")
library("reshape2")

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

# Input modified PLAC-seq file and return number of AND and XOR Interacitons
Get.ClusterType <- function(DF,NAME){
    TEMP_DF <- data.frame(table(unique(DF[,c('InteractionName','ClusterType')])[,'ClusterType']))
    TEMP_DF$SAMPLE <- NAME
    return(TEMP_DF)
                                    
}

# Cluster Types
OC_TYPE_DF <- Get.ClusterType(OC_DF, 'OPC')
MG_TYPE_DF <- Get.ClusterType(MG_DF, 'MG')
vRG_TYPE_DF <- Get.ClusterType(vRG_DF, 'vRG')
oRG_TYPE_DF <- Get.ClusterType(oRG_DF, 'oRG')

# Merge
PLOT_DF <- rbind(OC_TYPE_DF, MG_TYPE_DF, vRG_TYPE_DF, oRG_TYPE_DF)
head(PLOT_DF)

PLOT_DF$SAMPLE <- factor(PLOT_DF$SAMPLE , levels = c('vRG','oRG','OPC','MG'))

options(repr.plot.width=2, repr.plot.height=2)
ggplot(data=PLOT_DF, aes(x=Var1,y=Freq,fill=SAMPLE)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    theme(legend.title = element_text(color = "blue", size = 10),
            legend.text = element_text(color = "black"),
            legend.key.size = unit(.4, "cm"),
            axis.text.x = element_text(face="plain", color="black", 
                           size=6, angle=45,hjust=1),
            axis.text.y = element_text(face="plain", color="black", 
                           size=6, angle=0),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            strip.text.x = element_blank()) + 
    labs(fill = "",
         title = "",
         subtitle = '') + 
    scale_fill_manual(values=c("#2166AC","#B2182B","#1B7837", "#762A83"))
ggsave('Peak_Type_Barplot.pdf', height = 2, width = 2)