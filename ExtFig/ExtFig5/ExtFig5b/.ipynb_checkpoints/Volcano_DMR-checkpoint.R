df_new <- read.csv('RG.DeSEQ2.12.1.23.csv')
head(df_new)

#### Add DARs ####
vRG.DMRs <- read.table('../../Link_DARs_to_DEGs//vRG.DMRs_linked_DEGs.txt')
oRG.DMRs <- read.table('../../Link_DARs_to_DEGs//oRG.DMRs_linked_DEGs.txt')
DMRs <- rbind(vRG.DMRs, oRG.DMRs)
colnames(DMRs)  <- c('hgnc_symbol','DMRs')


df_new2 <- left_join(df_new, DMRs)

df_new2$DMRs[is.na(df_new2$DMRs)] <- 0

# empty column
df_new2$delabel <- NA
# Actually adding the labels
df_new2$delabel[df_new2$DMRs >= 1] <- df_new2[df_new2$DMRs >= 1,'hgnc_symbol']

require("ggrepel")
sample1 <- 'vRG'
sample2 <- 'oRG'
options(repr.plot.width = 4, repr.plot.height = 3)
p <- ggplot(df_new2) +
        geom_point(aes(x=log2FoldChange, y=-log10(padj), fill=cat, size=DMRs, stroke = 0.2), colour="black", pch=21) +
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
      scale_size_continuous(range  = c(0.5,3), 
                        limits = c(0,16), 
                        breaks = c(0, 1, 2, 4))
p


p
ggsave('DEG_with_DARs.pdf',width = 4, height=3)