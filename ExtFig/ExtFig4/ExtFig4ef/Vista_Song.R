#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)

# Read in unique ATAC-seq peaks
ReadCRE <- function(file){
    df<-read.table(file, header=F)
    colnames(df) <- c('chr','start','end')
    return(df)
}

# Perform fisher enrichment for cCREs and controls
enrichement_test <- function(bed, cCRE, PLAC, head=TRUE) {
    
    # Read in PLAC-seq
    PLAC_DF <- fread(PLAC)
    dim(PLAC_DF)
    head(PLAC_DF)

    # filter for only XOR interactions
    PLAC_DF <- PLAC_DF[PLAC_DF$type == 'XOR',]
    dim(PLAC_DF)
    head(PLAC_DF)
    
    # Read in cCREs
    cCRE_DF <- ReadCRE(cCRE)
    head(cCRE_DF)
    dim(cCRE_DF)

    # Retain only XOR interactions in with cCRE in distal region
    PLAC_cCRE_DF <- bt.intersect(a=PLAC_DF, b=cCRE_DF)
    dim(PLAC_cCRE_DF)
    head(PLAC_cCRE_DF)

    # Make test bed file
    test_df <- PLAC_cCRE_DF[,1:3]
    print('Test...')
    print(dim( test_df))
    print(head(test_df))
    
    # Make controls bed file  
    temp_df <- PLAC_cCRE_DF %>% separate('V4',into=c('chr','loc'),sep=':' ) %>% 
    separate('loc',into=c('start','end'),sep='-')
    # Get distance from anchor to distal element
    d <- as.numeric(temp_df$start) - test_df$V3
    d1 <- as.numeric(temp_df$start) - test_df$V2
    new_start <- (as.integer(temp_df$start)+d)
    new_end <- (as.integer(temp_df$start)+d1)
    control_df <- data.frame('chr'=temp_df$chr, 'start'=as.integer(new_start) , 'end'=as.integer(new_end))
    head(control_df)
    dim(control_df)

    # Remove any less than zero
    control_df<- control_df[!(control_df$start < 0),]
    head(control_df)
    dim(control_df)

    # Remove any controls also in test
    control_df <- bt.intersect(control_df, test_df , v=T)
    head(control_df)
    
    print('Control...')
    print(dim(control_df))
    print(head(control_df))
    
    # Make sure control and test don't overlap
    print('Make sure no overlap...')
    print(dim(bt.intersect(test_df, control_df,u=T)))
    dim(bt.intersect(test_df, test_df, u=T))
    dim(bt.intersect(control_df, control_df, u=T))
    
    # read in bed file
    bed_df <- read.table(bed, header=head)
    head(bed_df)
    
    # Get overlap with test
    a <- dim(bt.intersect(bed_df[,1:3],test_df,wa=T, u=T))[1]
    a
    b <- dim(bt.intersect(bed_df[,1:3],test_df,v=T))[1]
    b
    c <- dim(bt.intersect(bed_df[,1:3],control_df,wa=T, u=T))[1]
    c
    d <- dim(bt.intersect(bed_df[,1:3],control_df,v=T))[1]
    d
    
    final_df <- data.frame( 'Bed_Y'=c(a,c),
                      'Bed_N'=c(b,d))
    row.names(final_df) <- c('cCRE','Control')
    print(final_df)
    
    res <- fisher.test(final_df)
    
    return(c(a,c, b,d,  res$estimate,res$conf.int[1],res$conf.int[2],res$p.value, bed, cCRE, PLAC))
    
}


#### Positive vista overlap with Song et al.
df <- data.frame()
cCRE_list <- c('/shen/shenlabstore3/ijones1/Song2020/ATAC/RG.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/IPC.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/iN.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/eN.ATAC-seq.narrowPeak')
PLAC_list <- c(         '/shen/shenlabstore3/ijones1/Song2020/PLAC/RG.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/IPC.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/iN.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/eN.MAPS.peaks.mod.txt')

Cell_list <- c('RG_song','IPC_Song','iN_Song','eN_Song')

for (i in (1:length(PLAC_list))) {
    res <- c(enrichement_test('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Vista/file_7_vistarPosNeu.hg38.bed',cCRE_list[i],PLAC_list[i], head=F), Cell_list[i])
    
    df <- rbind(df,res)

}

colnames(df) <- c('cCRE_Overlap','Control_Overlap', 'cCRE_Non_Overlap','Control_Non_Overlap',
                      'OR','Low','High','p.val','Bed','cCRE','PLAC', 'Celltype')
df$p.val <-  as.numeric(df$p.val)
df$sig <- ''
df$sig[df$p.val <= 0.05] <- '*'
df$sig[df$p.val <= 0.01] <- '**'
df$sig[df$p.val <= 0.001] <- '***'

#### Negative vista overlap with Song et al.
df <- data.frame()
cCRE_list <- c('/shen/shenlabstore3/ijones1/Song2020/ATAC/RG.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/IPC.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/iN.ATAC-seq.narrowPeak',
         '/shen/shenlabstore3/ijones1/Song2020/ATAC/eN.ATAC-seq.narrowPeak')
PLAC_list <- c(         '/shen/shenlabstore3/ijones1/Song2020/PLAC/RG.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/IPC.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/iN.MAPS.peaks.mod.txt',
         '/shen/shenlabstore3/ijones1/Song2020/PLAC/eN.MAPS.peaks.mod.txt')

Cell_list <- c('RG_song','IPC_Song','iN_Song','eN_Song')

for (i in (1:length(PLAC_list))) {
    res <- c(enrichement_test('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Vista/file_8_negative.hg38.bed',cCRE_list[i],PLAC_list[i], head=F), Cell_list[i])
    
    df <- rbind(df,res)

}

colnames(df) <- c('cCRE_Overlap','Control_Overlap', 'cCRE_Non_Overlap','Control_Non_Overlap',
                      'OR','Low','High','p.val','Bed','cCRE','PLAC', 'Celltype')
df$p.val <-  as.numeric(df$p.val)
df$sig <- ''
df$sig[df$p.val <= 0.05] <- '*'
df$sig[df$p.val <= 0.01] <- '**'
df$sig[df$p.val <= 0.001] <- '***'



#### Plot results
df <- read.csv('Vista Enrichment - Song.csv')
head(df)

options(repr.plot.width=2.5, repr.plot.height=1.25)
ggplot(data=df, aes(y=Conditon, x=log2(OR), xmin=log2(Low), xmax=log2(High),color=Conditon )) +
  geom_point(size=2,position=position_dodge(width = .75)) + 
  # geom_linerange(size=1,position=position_dodge(width = 1), show.legend = FALSE) +
  geom_errorbarh(position=position_dodge(width = .75),height=.5) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  ggtitle('Vista Enhancer Enrichment') +
  ylab('')+
   xlim(-3.5, 3.5) +
  xlab('log2(Odds Ratio)') +
  theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 6),
       axis.title = element_text(size = 6),
       legend.position = "none") 
ggsave('Vista.Forest.Song.NH.6.3.24.pdf', width=2.5, height = 1.25)