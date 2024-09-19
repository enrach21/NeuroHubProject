#### Meta ####
# Author: Ian Jones
# Date: 2.2.23
# Email: Ian.Jones3@ucsf.edu

#### input ####
# Tss <- location of Tss file, formart ('chr','start','end', 'strand','ENST','ENSG','gene')
# BinnedGenome <- location of binned genome fromat ('chr','start','end')
# PLAC <- location of plac-seq bedpe file
# name <- character of column name in final output
# ATAC <- location of ATAC-seq peaks to filter PLAC-seq interaction between accessible regions

#### output ####
# dataframe, col1 is gene name, col2 is total interactions with gene

Count_TSS_ATAC.V3 <- function(Tss='/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed',
                      BinnedGenome='/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/BinnedGenome/hg38.2kb.bed',
                     PLAC,
                     Anchor,
                     name,
                     ATAC){
    
    # Read in binned genome
    BinnedGenomeDF <- read.table(BinnedGenome, header=FALSE)
    colnames(BinnedGenomeDF) <- c('chr','start','end')

    # Read in TSS  
    TssDF <- read.table(Tss, header=FALSE)
    colnames(TssDF) <- c('chr','start','end', 'strand','ENST','ENSG','gene')
    
    # Overlap 5kb with TSS
    BinnedGenomeTssDF <- unique(bt.intersect(a=BinnedGenomeDF, b=TssDF, wa=T, wb=T)[,c(1:3,9,10)])
    colnames(BinnedGenomeTssDF) <- c('chr','start','end','ENSG','gene')
    
    
    # Overlap TSS with anchors
    AnchorDF <- read.table(Anchor)
    BinnedGenomeTssDF <- unique(bt.intersect(a=BinnedGenomeTssDF, b=AnchorDF, wa=T, u=T))
    colnames(BinnedGenomeTssDF) <- c('chr','start','end','ENSG','gene')
    
    # Read in PLAC-seq file and only retain region 1 and region 2
    PLAC_DF <- read.table(PLAC, fill=T, header=T)[,1:6]
    print(paste0('Total PLAC-seq interactions ... ', dim(PLAC_DF)[1]))
    
    # Read in ATAC-seq
    ATAC_df <- read.table(ATAC)
    
    # Filter for peaks overlapping anchors
    ATAC_df <- bt.intersect(ATAC_df, AnchorDF, v=T) 
    
    # Use bed PairToBed to only get interactions between accessible regions
    PLAC_DF_ATAC <- unique(bt.pairtobed(PLAC_DF, ATAC_df,type='either')[1:6])
    print(paste0('Total PLAC-seq interactions between accessible reigons ... ', dim(PLAC_DF_ATAC)[1]))
    
    # Overlap TSS with regions 1
    temp1 <- bt.intersect(a=BinnedGenomeTssDF, b=PLAC_DF_ATAC[,1:3], wa=T )
    # Overlap TSS with regions 2
    temp2 <- bt.intersect(a=BinnedGenomeTssDF, b=PLAC_DF_ATAC[,4:6], wa=T )
    # Combine overlapping regions
    temp <- rbind(temp1, temp2)
    
    # Get dataframe of interaction counts
    df <- data.frame(table(temp$V5))
    df <- df[rev(order(df$Freq)),]
    colnames(df) <- c('hgnc_symbol', name)
    
    return(df)
}