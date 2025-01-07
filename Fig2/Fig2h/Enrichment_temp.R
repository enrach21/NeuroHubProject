#### Library #####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)




#### Plot enrichment ####

df <- read.csv('Vista Enrichment - Sheet1.csv')
df$shape <-c(rep('cCRE',4), rep('ATAC',4), rep('LMR',4))


df$shape <- factor(df$shape, levels = c('cCRE','ATAC','LMR'))

df$Conditon <- factor(df$Conditon, levels = rev(c('vRG_cCRE','vRG_ATAC','vRG_LMR',
                                             'oRG_cCRE','oRG_ATAC','oRG_LMR',
                                             'OPC_cCRE','OPC_ATAC','OPC_LMR',
                                             'MG_cCRE','MG_ATAC','MG_LMR')))

df <- df[df$shape == 'cCRE',]

colors2 <- rev(c('#2166AC',  
           '#B2182B', 
           '#1B7837',
           '#762A83'))

options(repr.plot.width=1, repr.plot.height=1)
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
 scale_y_discrete(labels = rep(rev(c('vRG','oRG','OPC','MG'))))
ggsave('Vista.FOrest.NH.11.26.24.pdf', width=2, height = 1.2)

