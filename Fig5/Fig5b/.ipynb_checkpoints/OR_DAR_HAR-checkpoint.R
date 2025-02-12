#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)


HAR <- read.table('GSE180714_HARs.bed', header=T)
head(HAR)

oRG_DAR <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG/oRG.DAR.bed')

vRG_DAR <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG/oRG.DAR.bed')

# 72 overlap
bt.intersect(oRG_DAR, HAR, wa=T)

# 13 overlap
bt.intersect(vRG_DAR, HAR, wa=T)


df <- read.csv('Vista Enrichment - HAR_DAR.csv')
head(df)



options(repr.plot.width=2, repr.plot.height=.9)
ggplot(data=df, aes(y=Conditon, x=log2(OR), xmin=log2(Low), xmax=log2(High) )) +
  geom_point(size=2,position=position_dodge(width = .75)) + 
  # geom_linerange(size=1,position=position_dodge(width = 1), show.legend = FALSE) +
  geom_errorbarh(position=position_dodge(width = .75),height=.5) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  ggtitle('HAR Enrichment') +
  ylab('')+
   xlim(-3.5, 3.5) +
  xlab('log2(Odds Ratio)') +
  theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0),
       axis.text = element_text(size = 8),
       axis.title = element_text(size = 8),
       legend.position = "none") +
 scale_y_discrete(labels = '')
ggsave('HAR.Forest.NH.8.20.24.pdf', width=2, height = .75)