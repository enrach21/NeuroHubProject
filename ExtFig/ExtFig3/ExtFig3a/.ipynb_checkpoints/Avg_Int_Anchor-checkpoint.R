### Packages ####
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


# Get the number of interactions for Anchor
Get.INT.Per.Anchor <- function(DF,NAME){
    TEMP_DF <- unique(DF[,c('InteractionName','TargetAnchor')])
    TEMP_DF <- TEMP_DF %>% count(TargetAnchor)
    TEMP_DF$SAMPLE <- NAME
    return(TEMP_DF)
}


# Cluster Types
OC_ANCHOR_COUNT_DF <- Get.INT.Per.Anchor(OC_DF, 'OPC')
MG_ANCHOR_COUNT_DF <- Get.INT.Per.Anchor(MG_DF, 'MG')
vRG_ANCHOR_COUNT_DF <- Get.INT.Per.Anchor(vRG_DF, 'vRG')
oRG_ANCHOR_COUNT_DF <- Get.INT.Per.Anchor(oRG_DF, 'oRG')

# Merge
PLOT_DF <- rbind(OC_ANCHOR_COUNT_DF, MG_ANCHOR_COUNT_DF, vRG_ANCHOR_COUNT_DF, oRG_ANCHOR_COUNT_DF)
head(PLOT_DF)

PLOT_DF$SAMPLE <- factor(PLOT_DF$SAMPLE , levels = c('vRG','oRG','OPC','MG'))

vline_df <- data.frame(SAMPLE = levels(as.factor(PLOT_DF$SAMPLE)),
                       Mean = tapply(X = PLOT_DF$n, INDEX = PLOT_DF$SAMPLE,
                                        FUN = mean))
vline_df 

options(repr.plot.width=4, repr.plot.height=4)
p <- ggplot(PLOT_DF, aes(x=n, fill=SAMPLE)) + 
      geom_histogram(alpha=1, position="identity")+ 
      scale_fill_manual(values=c("#2166AC","#B2182B","#1B7837", "#762A83")) +
      theme_classic() +
      xlim(c(0, 40)) +
      facet_wrap(~ SAMPLE) +
    labs(fill = "",
         title = "",
         subtitle = '') + geom_vline(data = vline_df, aes(xintercept = Mean), linetype = "dashed",
               colour = "black") +  
    geom_text(
      data    = vline_df,
      mapping = aes(x = 10, y = 5000, label = paste0('Mean = ', round(Mean, 2))),
      hjust   = -0.1,
      vjust   = -1,
      size =2) +
     theme(legend.title = element_text(color = "blue", size = 10),
            legend.text = element_text(color = "black"),
            legend.key.size = unit(.4, "cm"),
            axis.text.x = element_text(face="plain", color="black", 
                           size=6, angle=0),
            axis.text.y = element_text(face="plain", color="black", 
                           size=6, angle=0),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            strip.text.x = element_blank()) +
    labs(fill = "",
         title = "",
         subtitle = '') +
    xlab("# of Interactions per Anchor") + ylab("# of Anchors")
p
ggsave('Interactions_per_anchor.pdf', height = 4, width = 4)