#### Meta ####
# Author: Ian Jones
# Date: 2.2.23
# Email: Ian.Jones3@ucsf.edu

#### input ####
# df1 <- output from Count_TSS function
# df2 <- output from Count_TSS functions
# rpkm <- exprssion for iNGN2 Timecourse
# contrast <- which timepoint comparison default is week2 vs Day0

#### output ####


Plot_Corr <- function(df1, df2,
                      rpkm='/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/AVG_TMM_RPKM_NoLog_exp_geneID_4.13.23.txt', contrast=1){
    
    # Read in expression data
    RPKM <- read.csv(rpkm)
    
    # filter out genes not expressed at any timepoint
    RPKM_filt <- RPKM[
        (RPKM$MG > 1 |
        RPKM$Oligo > 1 |
        RPKM$oRG > 1 |
        RPKM$vRG > 1)
    ,]
    
    # Get difference between week2 and Day0
    if (contrast==1){
    # Get difference between MG
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$MG+1) - log2(RPKM_filt$Oligo+1)
    RPKM_filt <- RPKM_filt[
        (RPKM_filt$MG > 10 |
        RPKM_filt$Oligo > 10)
    ,]
    }
    # Get difference between week2 and Day3
    if (contrast==2){
    # Get difference between Week2 and Day3
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$MG+1) - log2(RPKM_filt$oRG+1)
    RPKM_filt <- RPKM_filt[
        (RPKM_filt$MG > 10 |
        RPKM_filt$oRG > 10 )
    ,]
    }
    # Get difference between Day3 and Day0
    if (contrast==3){
    # Get difference between Day3 and Day0
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$MG+1) - log2(RPKM_filt$vRG+1)
    RPKM_filt <- RPKM_filt[
        (RPKM_filt$MG > 10 |
        RPKM_filt$vRG > 10)
    ,]
    }
    if (contrast==4){
    # Get difference between Day3 and Day0
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$Oligo+1) - log2(RPKM_filt$oRG+1)
    RPKM_filt <- RPKM_filt[
        (
        RPKM_filt$Oligo > 10 |
        RPKM_filt$oRG > 10 )
    ,]
    }
     if (contrast==5){
    # Get difference between Day3 and Day0
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$Oligo+1) - log2(RPKM_filt$vRG+1)
    RPKM_filt <- RPKM_filt[
        (
        RPKM_filt$Oligo > 10 |
        RPKM_filt$vRG > 10)
    ,]
     }
    if (contrast==6){
    # Get difference between Day3 and Day0
    RPKM_filt$Log2_RPKM_Dif <- log2(RPKM_filt$oRG+1) - log2(RPKM_filt$vRG+1)
    RPKM_filt <- RPKM_filt[
        (
        RPKM_filt$oRG > 10 |
        RPKM_filt$vRG > 10)
    ,]
    }
    
    # Subset genenames and differencew
    RPKM_filt <- RPKM_filt[,c('hgnc_symbol','Log2_RPKM_Dif')]
    # print(head(RPKM_filt))
    
    # Get interaction changes of df2 over df1
    Int_df <- inner_join(df1, df2)
    colnames(Int_df) <- c('hgnc_symbol', 'Var1','Var2')
    Int_df$interaction_dif <- Int_df$Var2 -Int_df$Var1

    Int_df <- Int_df[,c('hgnc_symbol','interaction_dif')]
    # print(head(Int_df))
    
    df <- inner_join(RPKM_filt, Int_df)
    # print(head(df))
    
    # Remove blank gene elements
    df <- df[df$hgnc_symbol != '',]
    
    # Remove Inf and NA elements from RNA
    df <- df[!is.na(df$Log2_RPKM_Dif),]
    df <- df[df$Log2_RPKM_Dif != Inf,]
    # print(head(df))
    
    print(cor.test(df$interaction_dif, df$Log2_RPKM_Dif))
    
    print(cor.test(df$interaction_dif, df$Log2_RPKM_Dif)$p.value)
    
    p <- ggplot(df, aes(x=interaction_dif,y=Log2_RPKM_Dif)) +
    # geom_point()
      geom_pointdensity(size = .2) +
      scale_color_viridis() +
      theme_bw() +
       theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       # axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 8),
       legend.position = "none",
       axis.title=element_text(size=8)) +
    xlab(expression(Delta*" Interaction")) +
    ylab(expression(Delta*" log2(RPKM)")) + geom_smooth(method=lm, se=FALSE)

    return(list(df,p))
}