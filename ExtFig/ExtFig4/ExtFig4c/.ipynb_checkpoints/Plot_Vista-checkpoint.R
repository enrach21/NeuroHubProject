library(ggplot2)

# Read in Michael's Data
df <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/Michael_tested/2022-10-21_for_Ian.tsv', header=T)
df1 <- df[,c('chr','start','end','element','cell','tissue')]
df1$outcome <- 'Pos'
df1$outcome[df1$tissue == 'NEGATIVE'] <- 'Neg'
df1

df2 <- df1[,c('cell','outcome')]
df3 <- data.frame(table(df2))
df3$color <- paste0(df3$cell,'.',df3$outcome)
df3

df3$cell <- factor(df3$cell, levels = c('RG','IPC','eN','iN'))

options(repr.plot.width=2.5, repr.plot.height=2.5)
ggplot(data=df3, aes(x=cell, y=Freq, fill=outcome)) +
  geom_bar(stat="identity", width = 0.8)  +
        xlab(NULL) + 
        ylab('Vista Enhancers') +
         theme_classic() +
        theme(
            #legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=6, angle=30, hjust=1), 
              axis.title.y = element_text(face="plain", color="black", 
                           size=8)) +
   scale_fill_manual(values=c('blue','red')) +
    ylim(0,10)
ggsave('Vista.Michael.Song.5.31.24.pdf', width=2.5, height = 2.5)