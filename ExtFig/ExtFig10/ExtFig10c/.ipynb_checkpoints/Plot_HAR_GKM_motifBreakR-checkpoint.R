#### Read in packages #####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library("UpSetR")
library("ComplexHeatmap")
# library(gprofiler2)
library("ggrepel")
library(biomaRt)
library(dplyr)
library(corrplot)
library('seqinr')
library('stringr')

df <- fread('/shen/shenlabstore3/ijones1/dependencies/vcf2maf/HAR_MAF_V2/HAR.VAR.vcf', fill=T)[,1:5]
colnames(df) <- c('X.chrom','chromEnd','name','ref','alt')
df$chromStart <- df$chromEnd - 1
df <- df[, c('X.chrom','chromStart','chromEnd','name','ref','alt')]
dim(df)
head(df)

delta <- read.csv('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/HARs_update/HAR.8.12.24.csv', row.names=1)

df2 <- cbind(df, delta[,7:39])
df2$SNP_id <- paste0(df2$X.chrom,':',df2$chromEnd,':',df2$ref,':',df2$alt)

head(df2)
summary(df2$V2)

MotifBreakDf <- read.csv('HAR.Accessible_1e-4_.All.RPKM_gt_5.7.22.24.csv')
head(MotifBreakDf)
dim(MotifBreakDf)

# Get number of elements wih motif
length(unique(MotifBreakDf$SNP_id))

final_df <- left_join(df2, MotifBreakDf, by=c('SNP_id'))
head(final_df)
dim(final_df)

# Filter for only significant variants
final_df$Sig <- 'No'
final_df[final_df$explain_pval <= 0.05 & final_df$delta_p_val <= 0.05& final_df$ISM_pval <= 0.05,'Sig'] <- 'Yes'

dim(final_df)
head(final_df)

# Remove elements without motif
final_df2<- final_df[!is.na(final_df$alleleEffectSize),]

# Only keep one predicted motif
final_df3 <- final_df2 %>% group_by(name) %>% top_n(1, abs(alleleEffectSize))

final_df4 <-  final_df3
final_df4$lab <- NA
final_df4$lab[abs(final_df4$explain_score) > 0.8] <- final_df4$geneSymbol[abs(final_df4$explain_score) > 0.8]



ggplot(data=final_df4, aes(x=explain_score, y=alleleEffectSize, color=Sig)) + 
    geom_point(size=0.5) + theme_classic() +
    geom_label_repel(aes(x=explain_score, y=alleleEffectSize,label=lab), max.overlaps=30, size = 2.5, box.padding = 0.5) +
    geom_hline(yintercept=0, linetype="dashed", color = "red") +
    geom_vline(xintercept=0, linetype="dashed", color = "red") +
    scale_color_manual(values=c("grey", "#B2182B")) + 
    theme(legend.position = "none",
           axis.title = element_text(size = rel(1.25)),
                  axis.text.x = element_text(face="plain", color="black", 
                               size=6, angle=0, hjust=1),
                axis.text.y = element_text(face="plain", color="black", 
                   size=6, angle=0, hjust=1),
                
                  axis.title.x = element_text(face="plain", color="black", 
                               size=8),
                  axis.title.y = element_text(face="plain", color="black", 
                               size=8)) 
ggsave('Test.HAR.PLOT.8.20.24.pdf',height=3, width=3)