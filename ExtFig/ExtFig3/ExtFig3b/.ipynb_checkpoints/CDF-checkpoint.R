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

# Return the distance of each interaction
Get.Distance <- function(DF,NAME){
    TEMP_DF <- data.frame(unique(DF[,c('InteractionName','Distance')])[,'Distance'])
    colnames(TEMP_DF) <- 'DIST'
    TEMP_DF$SAMPLE <- NAME
    return(TEMP_DF)
}

# Cluster Types
OC_DIST_DF <- Get.Distance(OC_DF, 'OPC')
MG_DIST_DF <- Get.Distance(MG_DF, 'MG')
vRG_DIST_DF <- Get.Distance(vRG_DF, 'vRG')
oRG_DIST_DF <- Get.Distance(oRG_DF, 'oRG')

# Merge
PLOT_DF <- rbind(OC_DIST_DF, MG_DIST_DF, vRG_DIST_DF, oRG_DIST_DF)
head(PLOT_DF)

PLOT_DF$SAMPLE <- factor(PLOT_DF$SAMPLE , levels = c('vRG','oRG','OPC','MG'))

vline_df <- data.frame(SAMPLE = levels(as.factor(PLOT_DF$SAMPLE)),
                       Mean = tapply(X = PLOT_DF$DIST, INDEX = PLOT_DF$SAMPLE,
                                        FUN = mean))
vline_df 

options(repr.plot.width=2.5, repr.plot.height=2.5)
ggplot(data=PLOT_DF, aes(x=DIST, color=SAMPLE)) +
    stat_ecdf(linetype='solid',size=.2, alpha=2) + 
    scale_colour_manual(values=c("#2166AC","#B2182B","#1B7837", "#762A83")) +
    theme_classic() +
    theme(legend.title = element_blank(),
            legend.text = element_text(color = "black"),
            legend.key.size = unit(.4, "cm"),
            axis.text.x = element_text(face="plain", color="black", 
                           size=6, angle=45, hjust=1),
            axis.text.y = element_text(face="plain", color="black", 
                           size=6, angle=0),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8)) +

    geom_text(
      data    = vline_df[1,],
      mapping = aes(x = 500000, y = 0.5, label = paste0('Mean = ', round(Mean, 0), 'bp')),
      hjust   = -0.1,
      vjust   = -1,
      size =2) +

    geom_text(
          data    = vline_df[2,],
          mapping = aes(x = 500000, y = 0.4, label = paste0('Mean = ', round(Mean, 0), 'bp')),
          hjust   = -0.1,
          vjust   = -1,
          size =2) +


    geom_text(
          data    = vline_df[3,],
          mapping = aes(x = 500000, y = 0.3, label = paste0('Mean = ', round(Mean, 0), 'bp')),
          hjust   = -0.1,
          vjust   = -1,
          size =2) +


    geom_text(
          data    = vline_df[4,],
          mapping = aes(x = 500000, y = 0.2, label = paste0('Mean = ', round(Mean, 0), 'bp')),
          hjust   = -0.1,
          vjust   = -1,
          size =2) +

    labs(fill = "",
         title = "",
         subtitle = '') +
    xlab("Interaction Distance (kb)") + ylab("Cumulative Distribution") + 
geom_vline(data = vline_df, aes(xintercept = Mean), linetype = "dashed",
               colour = c("#2166AC","#B2182B","#1B7837", "#762A83"))
ggsave('PLAC_CDF.pdf', height = 2.5, width = 2.5)