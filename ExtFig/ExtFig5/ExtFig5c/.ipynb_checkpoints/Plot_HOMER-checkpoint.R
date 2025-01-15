# Read in packages
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library('reshape2')
library('limma')

Read.Homer <- function(file, cell) {
    df <- read.table(file, header=1,row.names=NULL)
    colnames(df)<- c('Motif','Consensus','P.value','Log.P.Value','q.value',
                 'num_of_Target','percent_of_Target','num_of_Background','percent_of_Background')
    df <- df[,c('Motif','P.value')]
    colnames(df)[2] <- cell
    return(df)
}

vRG <- Read.Homer('vRG_genome.Background_10_5_23/knownResults.txt','vRG')
oRG <- Read.Homer('oRG_genome.Background_10_5_23/knownResults.txt','oRG')

df <- full_join(vRG, oRG, by = 'Motif', multiple = 'all')

df$allias <- sapply(strsplit(df$Motif,"\\("), `[`, 1)
df$allias  <- toupper(df$allias )

df$gene_id <- alias2SymbolTable(df$allias , species = "Hs")
length(df$gene_id)

df_rna <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', row.names = 1)
df_rna2 <-  df_rna[,c('oRG','vRG','MG','Oligo','hgnc_symbol')]
colnames(df_rna2)[5] <- 'gene_id'
head(df_rna2)

#### Join together expression and HOMER results ####

df2 <- inner_join(melt(df), melt(df_rna2), by = c('gene_id','variable'), suffix = c('.homer','.RPKM'))

vRG_genes <- (df2 %>% filter((variable == 'vRG') & value.RPKM >10) %>% arrange(value.homer) %>% select(gene_id))[1:20,]
oRG_genes <- (df2 %>% filter((variable == 'oRG') & value.RPKM >10) %>% arrange(value.homer) %>% select(gene_id))[1:20,]

vRG_genes2 <- c('TCF12','TCF4','NEUROG2','ASCL1','FOS')
oRG_genes2 <- c('LHX2','PAX6','ETV1')

df3 <- df2[df2$gene_id %in% c(vRG_genes2,oRG_genes2),]
df3
df3$gene_id <- factor(df3$gene_id , levels = c(vRG_genes2,oRG_genes2))
df3$variable <- factor(df3$variable , levels = c('MG','Oligo','oRG','vRG'))
df3[df3 == 0] <- 1e-322


# oRG color
options(repr.plot.width=3.5, repr.plot.height=2)
ggplot(df3, aes(variable, gene_id)) +
  geom_point(aes(size = -log10(value.homer), colour=log2(value.RPKM+1))) +
  theme_classic()  +
  ylab('') +
  xlab('') +
        theme(
            #legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=0, hjust=0.5),
              axis.text.y = element_text(face="plain", color="black", 
                           size=8, angle=0, hjust=0.5),
        legend.key.size = unit(.2, "cm")) +
   scale_colour_gradient(
  low = "white",
  high = "black",
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "colour"
)
ggsave('Black.Homer.ATAC_MR.pdf', width = 3.5, height = 2)