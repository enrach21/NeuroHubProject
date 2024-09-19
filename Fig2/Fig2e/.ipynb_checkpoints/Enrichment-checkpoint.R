#### Library #####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)


#### Negative vista overlap ####

df <- data.frame()
cCRE_list <- c('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed',
              '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.ATAC.peak.bed',
            '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.MR.peak.bed')
PLAC_list <- c('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
'/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
              '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe')

Cell_list <- c('MG','OC','oRG','vRG','MG_ATAC','OC_ATAC','oRG_ATAC','vRG_ATAC','MG_MR','OC_MR','oRG_MR','vRG_MR')

for (i in (1:length(PLAC_list))) {
    res <- c(enrichement_test('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Vista/file_8_negative.hg38.bed',cCRE_list[i],PLAC_list[i], head=F), Cell_list[i])
    
    df <- rbind(df,res)

    

}

colnames(df) <- c('cCRE_Overlap','Control_Overlap', 'cCRE_Non_Overlap','Control_Non_Overlap',
                      'OR','Low','High','p.val','Bed','cCRE','PLAC', 'Celltype')
df$p.val <-  as.numeric(df$p.val)
df$sig <- ''
df$sig[df$p.val <= 0.05] <- '*'
df$sig[df$p.val <= 0.01] <- '**'
df$sig[df$p.val <= 0.001] <- '***'

#### Postive Vista Overlap ####

df <- data.frame()
cCRE_list <- c('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed',
              '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.ATAC.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.ATAC.peak.bed',
            '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/MG.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/OPC.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.MR.peak.bed',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.MR.peak.bed')
PLAC_list <- c('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
'/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
              '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe',
         '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe')

Cell_list <- c('MG','OC','oRG','vRG','MG_ATAC','OC_ATAC','oRG_ATAC','vRG_ATAC','MG_MR','OC_MR','oRG_MR','vRG_MR')

for (i in (1:length(PLAC_list))) {
    res <- c(enrichement_test('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Vista/file_7_vistarPosNeu.hg38.bed',cCRE_list[i],PLAC_list[i], head=F), Cell_list[i])
    
    df <- rbind(df,res)

}

colnames(df) <- c('cCRE_Overlap','Control_Overlap', 'cCRE_Non_Overlap','Control_Non_Overlap',
                      'OR','Low','High','p.val','Bed','cCRE','PLAC', 'Celltype')
df$p.val <-  as.numeric(df$p.val)
df$sig <- ''
df$sig[df$p.val <= 0.05] <- '*'
df$sig[df$p.val <= 0.01] <- '**'
df$sig[df$p.val <= 0.001] <- '***'


#### Plot enrichment ####

df <- read.csv('Vista Enrichment - Sheet1.csv')
df$shape <-c(rep('cCRE',4), rep('ATAC',4), rep('LMR',4))

df$shape <- factor(df$shape, levels = c('cCRE','ATAC','LMR'))

df$Conditon <- factor(df$Conditon, levels = rev(c('vRG_cCRE','vRG_ATAC','vRG_LMR',
                                             'oRG_cCRE','oRG_ATAC','oRG_LMR',
                                             'OPC_cCRE','OPC_ATAC','OPC_LMR',
                                             'MG_cCRE','MG_ATAC','MG_LMR')))

colors2 <- rev(c('#2166AC', 
            '#4393C3', 
           '#92C5DE', 
           '#B2182B', 
           '#D6604D', 
           '#F4A582', 
           '#1B7837', 
           '#5AAE61', 
           '#ACD39E',
           '#762A83', 
           '#9970AB', 
           '#C2A5CF'))

options(repr.plot.width=2, repr.plot.height=2)
ggplot(data=df, aes(y=Conditon, x=log2(OR), xmin=log2(Low), xmax=log2(High),color=Conditon )) +
  geom_point(size=2,position=position_dodge(width = .75), aes(shape=shape)) + 
  # geom_linerange(size=1,position=position_dodge(width = 1), show.legend = FALSE) +
  geom_errorbarh(position=position_dodge(width = .75),height=.5) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
 scale_color_manual(values=colors2)  + 
  ggtitle('Vista Enhancer Enrichment') +
  ylab('')+
   xlim(-3.5, 3.5) +
  xlab('log2(Odds Ratio)') +
  theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0),
       axis.text = element_text(size = 8),
       axis.title = element_text(size = 8),
       legend.position = "none") +
 scale_y_discrete(labels = rep(c('LMR','ATAC','cCRE'),4))
ggsave('Vista.FOrest.NH.5.31.24.pdf', width=2, height = 2)


#### Extended figure for with the Song et al data ####
df <- read.csv('Vista Enrichment - Song.csv')

options(repr.plot.width=2.5, repr.plot.height=1.25)
ggplot(data=df, aes(y=Conditon, x=log2(OR), xmin=log2(Low), xmax=log2(High),color=Conditon )) +
  geom_point(size=2,position=position_dodge(width = .75)) + 
  # geom_linerange(size=1,position=position_dodge(width = 1), show.legend = FALSE) +
  geom_errorbarh(position=position_dodge(width = .75),height=.5) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  ggtitle('Vista Enhancer Enrichment') +
  ylab('')+
   xlim(-3.5, 3.5) +
  xlab('log2(Odds Ratio)') +
  theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 6),
       axis.title = element_text(size = 6),
       legend.position = "none") 
ggsave('Vista.Forest.Song.NH.6.3.24.pdf', width=2.5, height = 1.25)