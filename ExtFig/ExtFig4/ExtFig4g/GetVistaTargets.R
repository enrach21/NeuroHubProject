# Author: Ian Jones
# Email: Ian.Jones3@ucsf.edu
# Date: 11.21.23
# Goal: Given input of neuronal vista enhancers, overlap all regions with ATAC-seq peaks and UMRs/LMRs and link them to their target genes with PLAC-seq


#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)


#### Read in Vista Elements ####
Neu <- read.table('file_7_vistarPosNeu.hg38.bed')
dim(Neu)
head(Neu)

#### Read in Data ####

# Read in unique ATAC-seq peaks
ReadATAC <- function(dir, file){
    df<-read.table(paste0(dir,file))
    colnames(df) <- c('chr','start','end')
    return(df)
}

# Read in unique methylation peaks
ReadMRs <- function(dir, file){
    df<-read.table(paste0(dir,file), header=T)
    return(df)
}


# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

OC <- 'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

MG <- 'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.bedpe'

# Read in PLAC
OC_DF <- read.table(paste0(DIR,OC), fill=T, header=T)

MG_DF <- read.table(paste0(DIR,MG), fill=T, header=T)

vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)


# ATAC-seq Overlapping

# Location of called peaks
ATAC_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/'

# Read in Data
MG_ATAC <- ReadATAC(ATAC_dir,'MG.overlap.optimal_peak.unique.bed')
nrow(MG_ATAC)
OPC_ATAC <- ReadATAC(ATAC_dir,'OPC.overlap.optimal_peak.unique.bed')
nrow(OPC_ATAC)
oRG_ATAC <- ReadATAC(ATAC_dir,'oRG.4reps.overlap.optimal_peak.unique.bed')
nrow(oRG_ATAC)
vRG_ATAC <- ReadATAC(ATAC_dir,'vRG.4reps.overlap.optimal_peak.unique.bed')
nrow(vRG_ATAC)

# Location of called peaks
MRs_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Stephanie_Nov_2023/UMR_LMR/'

MG_LMRs <- ReadMRs(MRs_dir,'MG_UMRsLMRs.bed')
nrow(MG_LMRs)
OPC_LMRs <- ReadMRs(MRs_dir,'OPC_UMRsLMRs.bed')
nrow(OPC_LMRs)
oRG_LMRs <- ReadMRs(MRs_dir,'oRG_UMRsLMRs.bed')
nrow(oRG_LMRs)
vRG_LMRs <- ReadMRs(MRs_dir,'vRG_UMRsLMRs.bed')
nrow(vRG_LMRs)

# Read in 2kb bin
bin <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/mc_Levels/hg38.2kb.bed'

Tss <- '/shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed'

RNA <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/EdgeR/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv'



#### Function to link Vista regions to their target genes
Summarize_Vista <- function(Neu,ATAC,MR,PLAC,ratio, bin, TSS,RNA,RNA_Cell, output, output2) {
  
    # Read in 2kb bin
    BinnedGenomeDF <- read.table(bin)
    colnames(BinnedGenomeDF) <- c('chr','start','end')
    head(BinnedGenomeDF)
    
    # Read in TSS  
    TssDF <- read.table(Tss, header=FALSE)
    colnames(TssDF) <- c('chr','start','end', 'strand','ENST','ENSG','gene')
    
    # Overlap 5kb with TSS
    BinnedGenomeTssDF <- unique(bt.intersect(a=BinnedGenomeDF, b=TssDF, wa=T, wb=T)[,c(1:3,10)])
    colnames(BinnedGenomeTssDF) <- c('chr','start','end','Local_gene')
    # dim(BinnedGenomeTssDF)
    # head(BinnedGenomeTssDF)

    # Collapse TSS within same bin
    BinnedGenomeTssCollapsedDF <- BinnedGenomeTssDF %>% group_by(chr, start, end) %>% 
      mutate(Local_gene = as.character(Local_gene),
         Local_gene = ifelse(length(Local_gene)>=2, 
                          paste(Local_gene, collapse='|'), Local_gene) 
            ) %>%
      unique
    # dim(BinnedGenomeTssCollapsedDF)
    # head(BinnedGenomeTssCollapsedDF)

    # Merge with original bin
    BinnedGenomeAnnoDF <- full_join(BinnedGenomeDF,BinnedGenomeTssCollapsedDF, by=c('chr','start','end'))
    print('Dimensions of annotated binned genome...')
    dim(BinnedGenomeAnnoDF)

    
    # Overlap with Methyalted regions and accessible regions
    Neu_temp1 <- bt.intersect(Neu, ATAC,wa=T, f=ratio)
    Neu_temp2 <- bt.intersect(Neu_temp1, MR,wa=T, f=ratio)

    
    colnames(Neu_temp2) <- c('chr_vista','start_vista','end_vista','Name','Vista')


    # Get the target anchors
    df <- bt.intersect(Neu_temp2, bin, wa=T,wb=T, f=ratio)
    df$TargetAnchor <- paste0(df$V6,':',df$V7,'-',df$V8)
    df <- df[,c(1:4,9)]
    colnames(df) <- c('chr_vista','start_vista','end_vista','Name','TargetAnchor')
    head(df)
    dim(df)

    # Left join with PLAC-seq distally
    df2 <- unique(bt.intersect(a=df[,1:4], b=PLAC,wa=T, wb=T, f=ratio))
    colnames(df2) <- c(colnames(df)[1:4],colnames(PLAC)) 


    # Left join with PLAC-seq overlapping with promoters
    df3 <- unique(left_join(df, PLAC, by = 'TargetAnchor',relationship =  "many-to-many"))
    df3 <- df3[!is.na(df3$count),]

    df4 <- unique(full_join(Neu_temp2, df2))
    df4 <- unique(full_join(df4, df3))
    
    # Add which promoter the postive Vista overlaps with
    df5 <- left_join(df4, BinnedGenomeAnnoDF)
    
    write.csv(df5, output, row.names=F)
    
    # Return the target genes and local genes for each vista elements

    # Get the Target Genes of the Vista elements
    Target <- (df5[,c('chr_vista','start_vista','end_vista','Name','TargetGene')])

    # Separate the Target genes by '|'
    Target1 <- separate_rows(Target,TargetGene,sep="\\|")

    # Remove any Target gene that are NA
    Target2 <-   Target1[ !is.na(Target1$TargetGene),]

    # Remove any Target genes that are ''
    Target3 <-   Target2[ Target2$TargetGene != '',]

    # Filter for genes that are expressed
    rna <- read.csv(RNA, row.names = 1)
    rna <- unique(rna[,c(RNA_Cell,'hgnc_symbol')])
    colnames(rna) <- c(RNA_Cell,'TargetGene')
    head(rna)
    dim(rna)

    ### Filter for expressed genes
    rna <- rna[rna[RNA_Cell] >= 1 & !is.na(rna[RNA_Cell]),]
    dim(rna)

    Target4 <- inner_join(Target3,rna)

    Target4 <- Target4[,-6]

    # Recollapsed Target Genes
    Target5 <- Target4 %>% group_by(chr_vista, start_vista, end_vista, Name) %>% 
      mutate(TargetGene = as.character(TargetGene),
            TargetGene = ifelse(length(TargetGene)>=2, 
                    paste(unique(TargetGene), collapse='|'), TargetGene), 
            ) %>%
        unique


    # Get the Local Genes of the Vista elements and repeat the steps as above for local genes
    Local <- (df5[,c('chr_vista','start_vista','end_vista','Name','Local_gene')])

    Local1 <- separate_rows(Local,Local_gene,sep="\\|")

    Local2 <-  Local1[ !is.na(Local1$Local_gene),]

    Local3 <-  Local2[ Local2$Local_gene != '',]

    colnames(rna) <- c(RNA_Cell,'Local_gene')
    Local4 <- inner_join(Local3,rna)
    Local4 <- Local4[,-6]

    Local5 <- Local4 %>% group_by(chr_vista, start_vista, end_vista, Name) %>% 
      mutate(Local_gene = as.character(Local_gene),
         Local_gene = ifelse(length(Local_gene)>=2, 
                          paste(unique(Local_gene), collapse='|'), Local_gene) 
            ) %>%
        unique

    df6 <- full_join(Target5,Local5)

     write.csv(df6, output2, row.names=F)

    
    
}

#### Get reigons ####

ATAC_list <- list(MG_ATAC,OPC_ATAC,oRG_ATAC,vRG_ATAC)
MR_list <- list(MG_LMRs,OPC_LMRs,oRG_LMRs,vRG_LMRs)
PLAC_list <- list(MG_DF,OC_DF,oRG_DF,vRG_DF)
RNA_Cell_list <- list('MG','Oligo','oRG','vRG')
cell_list <- list('MG','OPC','oRG','vRG')
ratio <- 0.01


res <- data.frame()
for (x in 1:length(cell_list)) {
        ATAC <- ATAC_list[[x]]
        MR <- MR_list[[x]]
        PLAC <- PLAC_list[[x]]
        RNA_Cell <- RNA_Cell_list[[x]]
        cell <- cell_list[[x]]
    
    
        Summarize_Vista(Neu,ATAC,MR,PLAC,ratio, bin, Tss,RNA,RNA_Cell, output=paste0(cell,'_Vista_Interactions.csv'), output2=paste0(cell,'_Vista_Genes.csv'))
    }