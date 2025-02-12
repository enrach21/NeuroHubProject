#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)

#### Functions ####

# Read in unique ATAC-seq peaks
ReadCRE <- function(file){
    df<-read.table(file, header=F)[,1:3]
    colnames(df) <- c('chr','start','end')
    return(df)
}

Plot_enrichment_V2 <- function(BED, POP, HAR, HAR2, name, color, n=100) {
  
    ####  Read in HARs ####
    
    # All HARs
    HAR_DF <- ReadCRE(HAR)
    HAR_DF <- HAR_DF[-1,]
    dim(HAR_DF)
    head(HAR_DF)
    
    # HARs with activity
    HAR2_DF <- read.table(HAR2, header=T)
    head(HAR2_DF)
    HAR2_DF  <- HAR2_DF[HAR2_DF$caMPRA == 1,]
    dim(HAR2_DF)
    head(HAR2_DF)
    
    # Human specific HARs
    HAR3_DF  <- HAR2_DF[HAR2_DF$caMPRA_up == 'human',]
    dim(HAR3_DF)
    head(HAR3_DF)

    # Chimp Specific HARs
    HAR4_DF  <- HAR2_DF[HAR2_DF$caMPRA_up == 'chimp',]
    dim(HAR4_DF)
    head(HAR4_DF)
    
    # Zoo Hars (https://www.science.org/doi/10.1126/science.abm1696#supplementary-materials)
    Zoo_DF <- read.csv('science.abm1696_table_s1.csv')[,1:4]
    dim(Zoo_DF)
    head(Zoo_DF)
    
    BED_DF <- ReadCRE(BED)
    print('Number of peaks...')
    print(dim(BED_DF))
    head(BED_DF)
    
    POP_DF <- ReadCRE(POP)
    print('Number of peaks to pick from...')
    print(dim(POP_DF))
    head(POP_DF)
    
    set.seed(1) 
    sample_dist <- c()
    for (x in 1:n) {
        # Get random peaks equal to the number of DARs
        sample_DF <- POP_DF[sample(nrow(POP_DF), size=nrow(BED_DF)), ]
        # Overlap random peaks with HARs
        overlap <- nrow(bt.intersect(a=HAR_DF,b=sample_DF, wa=T, u=T))
        # Append to vecto
        sample_dist <- c(sample_dist, overlap)
    }
    
  print('All 3k HARs')
    
    print(paste0('Total number of overlaps: ', nrow(bt.intersect(a=HAR_DF,b=BED_DF, wa=T, u=T))))
    
    score <- (nrow(bt.intersect(a=HAR_DF,b=BED_DF, wa=T, u=T)) - mean(sample_dist)) / (sd(sample_dist))
    print(paste0('Total z-score: ', score))

    pval <- sum(abs(score) < abs(scale(sample_dist))) / length(sample_dist)
    print(paste0('P-val: ', pval))

    
    options(repr.plot.width=3, repr.plot.height=3)
    ggplot(data.frame(scale(sample_dist)), aes(x=scale.sample_dist.)) + 
    #   geom_histogram(binwidth=.2,color="black", fill="#B2182B", alpha=.25, bins = 500) +
     geom_density(alpha=0.25, fill=color) +
      geom_vline(aes(xintercept=((nrow(bt.intersect(a=HAR_DF,b=BED_DF, wa=T, u=T)) - mean(sample_dist)) / (sd(sample_dist)))),
                color=color, linetype="dashed", size=1) +
      theme_classic() +
            theme(
                #legend.position = "none",

                  axis.title = element_text(size = rel(1.25)),
                  axis.text.x = element_text(face="plain", color="black", 
                               size=6, angle=0, hjust=1),
                
                  axis.title.x = element_text(face="plain", color="black", 
                               size=8),
                  axis.title.y = element_text(face="plain", color="black", 
                               size=8)) +
     xlab('Z-score') +
    xlim(-11,11)

    ggsave(paste0(name,'.HAR.cCRE.pdf'),height = 1.5, width = 1.5)
}
    
    
#### Get enrichment for unique LMAR ####
Plot_enrichment_V2('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Intervene_Both/sets/0100_oRG.bed',
               '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Total.both.peak.merge.bed',
               'GSE180714_HARs.bed',
               'GSE180714_HARs.tsv',
               'oRG_cCRE',
                "#B2182B",
                1000
               )

Plot_enrichment_V2('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Intervene_Both/sets/1000_vRG.bed',
               '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Total.both.peak.merge.bed',
               'GSE180714_HARs.bed',
               'GSE180714_HARs.tsv',
               'vRG_cCRE',
                "#2166AC",
                1000
               )

Plot_enrichment_V2('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Intervene_Both/sets/0010_OPC.bed',
               '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Total.both.peak.merge.bed',
               'GSE180714_HARs.bed',
               'GSE180714_HARs.tsv',
               'OPC_cCRE',
                "#1B7837",
                1000
               )

Plot_enrichment_V2('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Intervene_Both/sets/0001_MG.bed',
               '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Total.both.peak.merge.bed',
               'GSE180714_HARs.bed',
               'GSE180714_HARs.tsv',
               'MG_cCRE',
                "#762A83",
                1000
               )