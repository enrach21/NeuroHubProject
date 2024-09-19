#### Read in packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library('reshape2')
library('limma')
library(RColorBrewer)

#### Function ####
Read.Homer <- function(file, cell) {
    df <- read.table(file, header=1,row.names=NULL)
    colnames(df)<- c('Motif','Consensus','P.value','Log.P.Value','q.value',
                 'num_of_Target','percent_of_Target','num_of_Background','percent_of_Background')
    df <- df[,c('Motif','P.value')]
    colnames(df)[2] <- cell
    return(df)
}

#### Workflow ####

# Read in results
vRG <- Read.Homer('../Homer/vRG_5_17_24/knownResults.txt','vRG')
oRG <- Read.Homer('../Homer/oRG_5_17_24/knownResults.txt','oRG')
OPC <- Read.Homer('../Homer/OPC_5_17_24/knownResults.txt','Oligo')
MG <- Read.Homer('../Homer/MG_5_17_24/knownResults.txt','MG')


df <- full_join(vRG, OPC, by = 'Motif', multiple = 'all')

df <- full_join(df, MG, by = 'Motif', multiple = 'all')

df <- full_join(df, oRG, by = 'Motif', multiple = 'all')


df$allias <- sapply(strsplit(df$Motif,"\\("), `[`, 1)
df$allias  <- toupper(df$allias )


df$gene_id <- alias2SymbolTable(df$allias , species = "Hs")
length(df$gene_id)

#### Read in RNA ####

df_rna <- read.csv('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv', row.names = 1)
df_rna2 <-  df_rna[,c('oRG','vRG','MG','Oligo','hgnc_symbol')]
colnames(df_rna2)[5] <- 'gene_id'
head(df_rna2)


#### Combine the datasets ####
df2 <- inner_join(melt(df), melt(df_rna2), by = c('gene_id','variable'), suffix = c('.homer','.RPKM'))


#### Filter for TF expressed greater than 10 ####
vRG_genes <- (df2 %>% filter((variable == 'vRG') & value.RPKM >10) %>% arrange(value.homer) %>% select(gene_id))[1:10,]
oRG_genes <- (df2 %>% filter((variable == 'oRG') & value.RPKM >10) %>% arrange(value.homer) %>% select(gene_id))[1:10,]
MG_genes <- (df2 %>% filter((variable == 'MG') & value.RPKM >10)  %>% arrange(value.homer)  %>% select(gene_id))[1:10,]
OPC_genes <- (df2 %>% filter((variable == 'Oligo') & value.RPKM >10)  %>% arrange(value.homer) %>% select(gene_id))[1:10,]


# Follow up by hand selecting TFs of interest
RG_genes2 <- c('LHX2','RFX2','POU3F3','ZIC2')
MG_genes2 <- c('IRF8','IRF3','SPI1','ELF1')
OPC_genes2 <-  c('TCF12','OLIG2','ASCL1','SOX10')
All_genes2 <- c('CTCF')

# Make new DF for plotting
df3 <- df2[df2$gene_id %in% c(RG_genes2, MG_genes2,OPC_genes2,All_genes2),]
df3
df3$gene_id <- factor(df3$gene_id , levels = rev(c(RG_genes2,OPC_genes2, MG_genes2,All_genes2)))
df3$variable <- factor(df3$variable , levels = rev(c('MG','Oligo','oRG','vRG')))
df3[df3 == 0] <- 1e-322


# Color
RdBu_seq <- seq(0,1, length.out=9)
length(RdBu_seq)
RdBu <-rev(brewer.pal(9,"RdBu"))
length(RdBu)


#### Plot vertical ####
options(repr.plot.width=1.5, repr.plot.height=3.75)
ggplot(df3, aes(variable, gene_id)) +
  geom_point(aes(size = -log10(value.homer), colour=log2(value.RPKM+1))) +
  theme_classic()  +
  ylab('') +
  xlab('') +
        theme(
            #legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=45, hjust=.5),
              axis.text.y = element_text(face="plain", color="black", 
                           size=8, angle=0, hjust=0.5),
             legend.position="bottom", legend.box = "vertical",
        legend.text = element_text(size=6),
        legend.title = element_text(size=7),
        legend.key.size = unit(.2, 'cm')) +
   scale_colour_gradientn(values=YlOrBr_seq, colours = YlOrBr) + 
  scale_x_discrete(labels= c('vRG
(n=8,343)'
                             ,'oRG
(n=5,279)'
                             ,'OPC
(n=11,697)'
                             ,'MG
(n=32,811)'))+ 
    guides(colour = guide_colourbar(title.position = "top"), 
          size = guide_legend(title.position = "top")) +
 scale_size_continuous(range = c(0, 2.5))

ggsave('Final.Homer.ATAC_MR.pdf', width = 1.5, height = 3.75)

#### Plot Horizontal ####
options(repr.plot.width=4, repr.plot.height=1.25)
ggplot(df3, aes( gene_id,variable)) +
  geom_point(aes(size = -log10(value.homer), colour=log2(value.RPKM+1))) +
  theme_classic()  +
  ylab('') +
  xlab('') +
        theme(
            #legend.position = "none",

              axis.title = element_text(size = rel(1.25)),
              axis.text.x = element_text(face="plain", color="black", 
                           size=8, angle=45, hjust=1),
              axis.text.y = element_text(face="plain", color="black", 
                           size=8, angle=0, hjust=0.5),
        legend.box = "horizontal",
        legend.direction = "vertical",
        legend.position = "right",
        #legend.title=element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=7),
        legend.key.size = unit(.2, 'cm')) + 
  scale_y_discrete(labels= rev(c('vRG'
                             ,'oRG'
                             ,'OPC'
                             ,'MG'))) +
   scale_colour_gradientn(values=YlOrBr_seq, colours = YlOrBr)  +
 scale_size_continuous(range = c(0, 2.5)) + 
    guides(colour = guide_colourbar(title= 'log2(RPKM)', title.position = "top"), 
          size = guide_legend(title= 'log10(p.val)',title.position = "top"))


ggsave('Transposed.Final.Homer.ATAC_MR.pdf', width = 4, height = 1.25)