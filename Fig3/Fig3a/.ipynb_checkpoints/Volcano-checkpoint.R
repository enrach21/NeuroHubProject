# Make Volcano plot


# Load libraries
library(edgeR)
library("DESeq2")
library("tximport")
library("ggplot2")
library("pheatmap")
library(dplyr)
library(tidyr)
library(biomaRt)
library(reshape2)
library("limma")

# Sample info
file <- 'RG.txt'
name <- 'RG'

# Color info
my_colors <- c('#762A83','#1B7837','#B2182B','#2166AC', 'grey','grey')
names(my_colors) <- c('MG','OPC','oRG','vRG', 'NOT_DAR','NOT_SIG')
my_colors

# Read in samples
dir <- '.' # location of sample.txt
samples <- read.table(file.path(dir,file), header=TRUE)
samples


# Read in Counts
dir = '/wynton/group/shen/NeuroHub.7.23.21/RNA-seq/STAR_RSEM_9_7_23/genes.results/' # location of RSEM output
files <- file.path(dir, paste0(samples$samples, ".genes.results"))
names(files) <- samples$samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

# Filter low quality reads and create design
txi.rsem$length[txi.rsem$length == 0] <- 1 # Make it so genes length cannot be zero
ddsTxi <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = samples,
                                   design = ~ date + condition
                                  )

# Perform pre-filtering
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

# Perfrom differential expression
dds <- DESeq(dds)

# Get DEGs
sample1 <- 'vRG'
sample2 <- 'oRG'
res <- results(dds, contrast=c("condition",sample1,sample2), alpha=0.05)
head(res)
summary(res)

resOrdered <- res[order(res$pvalue),]
resOrdered$cat <- "NOT_SIG"
resOrdered$cat[resOrdered$log2FoldChange >= 0.5 & resOrdered$padj < 0.05]  <- sample1
resOrdered$cat[resOrdered$log2FoldChange <= -0.5 & resOrdered$padj < 0.05]  <- sample2
sum(resOrdered$cat == sample2)

write.csv(resOrdered, 'RG.filt.DeSEQ2.Total.6.5.24.csv')


# Get Gene Symbols
df <- data.frame(resOrdered)
df$ENSMBL_ID <- row.names(df)

egids <- unlist(strsplit(as.character(df$ENSMBL_ID), split="\\."))[(1:length(df$ENSMBL_ID))*2-1]
head(egids)
df$ensembl_gene_id <- egids

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl",  host = "https://useast.ensembl.org/")

translated <- getBM(filters="ensembl_gene_id", 
                    attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                    values=egids, mart=ensembl)

df_new<- unique(inner_join(df, translated[,1:2], multiple='all'))
df_new <- df_new[!is.na(df_new$padj),]
dim(df_new)
head(df_new, n=20)

write.csv(df_new, 'RG.DeSEQ2.12.1.23.csv')

df_new <- read.csv('RG.DeSEQ2.12.1.23.csv')
head(df_new)

#### Add DARs ####
vRG.DARs <- read.table('../../Link_DARs_to_DEGs//vRG.DARs_linked_DEGs.txt')
oRG.DARs <- read.table('../../Link_DARs_to_DEGs//oRG.DARs_linked_DEGs.txt')
DARs <- rbind(vRG.DARs, oRG.DARs)
colnames(DARs)  <- c('hgnc_symbol','DARs')

df_new2 <- left_join(df_new, DARs)

# empty column
df_new2$delabel <- NA
# Actually adding the labels
genes <- c('MOXD1','LIFR','FBN2','PDGFC','FAM107A','SLCO1C1','HOPX','SPRY2','ETV5',
          'CRYAB','CCN2','FOXJ1','FBX032','INSM1','FOSB','EGFR','EGR1','NR4A1')
df_new2$delabel[df_new2$hgnc_symbol %in% genes] <- df_new2[df_new2$hgnc_symbol %in% genes,'hgnc_symbol']

require("ggrepel")
sample1 <- 'vRG'
sample2 <- 'oRG'
options(repr.plot.width = 4, repr.plot.height = 3)
df_new2[df_new2$DARs == 0.5, 'DARs'] <- 0
p <- ggplot(df_new2) +
        geom_point(aes(x=log2FoldChange, y=-log10(padj), fill=cat, size=DARs, stroke = 0.2), colour="black", pch=21) +
        geom_label_repel(aes(x=log2FoldChange, y=-log10(padj),label=delabel), max.overlaps=30, size = 2, box.padding = 0.5) +
        ggtitle(paste0(sample1, ' expression over ', sample2)) +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = 8, hjust = 1),
              axis.title = element_text(size = rel(1))) +
        theme_classic() + scale_fill_manual(values=my_colors) + 
        theme(plot.title = element_text(size=8, hjust = 0.5)) +
        theme(legend.title = element_blank(), legend.key.size = unit(.5, 'cm')) + 
        theme(legend.text=element_text(size=8)) + 
        guides(colour = guide_legend(override.aes = list(size=2))) +
        geom_vline(xintercept = 0.5, linetype="dotted", 
                color = "red", linewidth=.3) +
        geom_vline(xintercept = -0.5, linetype="dotted", 
                color = "red", linewidth=.3) +
        geom_hline(yintercept = 1.3, linetype="dotted", 
                color = "red", linewidth=.3)+
      scale_size_continuous(range  = c(0.5,4), 
                        limits = c(0,16), 
                        breaks = c(0, 1, 2, 4))
p

p
ggsave('DEG_with_DARs.pdf',width = 4, height=3)