# load libraries
library('bedtoolsr')
library('dplyr')
library('ggplot2')
library('data.table')
library('pheatmap')



df <- read.csv('LDSC_figure - Sheet1 (2).csv')


# Add significant
# Add significants
df$sig <- ''
df$sig[df$Enrichment_p <= 0.05] <- '*'
df$sig[df$Enrichment_p <= 0.01] <- '**'
df$sig[df$Enrichment_p <= 0.001] <- '***'

df$Celltype <- factor(df$Celltype, levels = c('vRG','oRG','OC','MG'))

df$Trait <- factor(df$Trait, levels = rev(c('AD','ADHD','ASD','DS','MDD','NEU','PD','SCZ')))

head(df)


my_colors <- c('#762A83','#1B7837','#B2182B','#2166AC')
names(my_colors) <- c('MG','OPC','oRG','vRG')
my_colors


#quantile_breaks
quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(df[df$Annot=='ATAC.MRs','EnrichmentScore'], n = 6)
my.colors <- colorRampPalette(colors = c("#2166AC","#F7F7F7", "#B2182B"))(length(mat_breaks))



options(repr.plot.width=3, repr.plot.height=2)
ggplot(df[df$Annot=='ATAC.MRs',], aes(y = Trait, x = Celltype, fill = (EnrichmentScore))) +
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  theme_minimal() +
  guides(fill = guide_colourbar(title = "Log(Score)")) +
  geom_text(aes(label = sig), color = "black", size = 4)  +
    theme( plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
       axis.text.x = element_text( angle=0, hjust = 0.5, size=8),
       axis.text = element_text(size = 8),
       axis.title=element_text(size=8),
    legend.key.size = unit(0.25, 'cm')) +
    xlab('') +
    ylab('') +
#    scale_fill_gradient2(  low = ("blue"),
#                           mid = "white",
#                           high = ("red"),
#                           midpoint = 0,
#                           space = "Lab",
#                          na.value = "grey50")
  scale_fill_gradientn(colors = c('#1965B0','white','#DC050C'), 
                       breaks=c(-25,0,25),
                       labels=c(-25,0,25),
                       limits=c(-25,25),
                      na.value = '#DC050C')
ggsave('ATAC.MRs.LDSC.pdf',width=3,height=2)

