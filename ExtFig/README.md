# Characterizing the Non-Coding Genome of Glial Cells in the Developing Human Cortex
 
### ExtFig1: Reproducibility of multi-omic datasets
    ExtFig1a: Hierarhcial cluster of transcriptome created with DeSEQ2
    ExtFig1b: Spearman correlation of ATAC-seq reads within unified peakset
    ExtFig1c: CpG clustering (To be Done)
    ExtFig1d: Pearson correlation of PLAC-seq libraries created with HPrep
    ExtFig1e: Civersort prediction of cell distribution imputed from Nowakaski et al. 2016
    
### ExtFig2: Defining cCREs with ATAC-seq and WGBS.
    ExtFig2a: Barplot of how cCREs are defined
    ExtFig2bc: ATAC-seq signal (b) and Cpg methylation percentage (c) of LMAR, AR, and LMR cCREs    
    ExtFig2d: Upset plot of H3K4me3 anchors
    ExtFig2e: Barplot of LMAR overlapping H3K4me3 anchors
    
### ExtFig3: Features of the 3D epigenome and the role of cCRELMAR on chromatin interactions.
    ExtFig3a: Barplot of number of interactions per an anchor bin 
    ExtFig3b: Cumulative Distribution of interaction distance for all 4 cell types
    ExtFig3c: Overview of interaction clusters
    ExtFig3d: Violin plot of observed over expected counts for LMAR interaction vs other interactions
    ExtFig3e: Barplot showing the percentage of clusters containing a LMAR
    ExtFig3f: Violin plot of observed over pval for LMAR clusters vs other clusters
    
### ExtFig4: Vista enhancers involved in the 3D epigenome of the developing brain. 
    ExtFig4a: Upset plot of significant interactions
    ExtFig4b: Heatmap of correlation between RPKM and Interaction difference for all distal interacitons and distal interactions with LMAR in the distal bin.
    ExtFig4c: Overview of initial vista results
    ExtFig4d: Neural and negative vista elements for global analysis
    ExtFig4ef: Fisher exact test for association between interacting cCREs and vista enhancer function.
    ExtFig4g: Linking Neu Vista to their potential target genes.
    
### ExtFig5: Identifying differential chromatin accessibility, methylation, and motif enrichment 
    ExtFig5a: Venn diagram of overlap with vRG and oRG DAR with psuedo-bulk DARs from Ziffra et. al. Nature 2021.
    ExtFig5b: Heatmap of DEGs linked with DMRs
    ExtFig5c: HOMER motif enrichment results within vRG and oRG DARs
    ExtFig5d: Footprinting difference between OPC and MG
    ExtFig5e: Correltation of motif prediction difference and expression

### ExtFig6: Footprinting analysis
    ExtFig6ab: Footrpinting binding plots from TOBIAS for LHX2 and ASCL1
    ExtFig6c: Cyctoscape for ASCL1
    ExtFig6d: Violin plot for ASCL1 target genes
    
### ExtFig7: Overview of GKM model used to evaluate variants
    ExtFig7a: LDSC for distal PLAC-seq bins
    ExtFig7b: Area under the receiver operating characteristic curve for GKM-SVM model cell type accessiblity predicitons
        - Selecting_Train_and_test
            -Script to obtain train and testing regions
        - Model_Test
            - Script to evalute model when trained on all chromosome except chr 2 and evaluated on chr 2
        - Final_Model_Training
            - Code to perform final training of the model
    ExtFig7c: Overview of AD SNP prioritization
        - Analyze_SNPs.200bp.R
            - Script to obtain 200bp surrounding 
        - Model_Eval
            - Contains script for predicitons for all AD variants across the four celltypes
    ExtFig7d: Overview of SCZ prioritizaiton
        - Analyze_SNPs.200bp.R
            - Script to obtain 200bp surrounding 
        - Model_Eval
            - Contains script for predicitons for all AD variants across the four celltypes
    ExtFig7e: Comparision of GKM delta, explain and ISM prediction for AD variants using the MG model
    ExtFig7f: Comparision of GKM delta, explain and ISM prediction for SCZ variants using the vRG model

### ExtFig8: Prioritized variants for AD and SCZ
    ExtFig8a: AD variants predicted to lead to less accessibilty
    ExtFig8b: AD variants predicted to lead to greater accessibilty
    ExtFig8c: GKMexplain score within MG for rs636317
    ExtFig8e: SCZ variants predicted to lead to less accessibilty
    ExtFig8f: SCZ variants predicted to lead to greater accessibilty
    ExtFig8g: Browser of signal at rs4449074 locus for heterozygous and homozygous individuals

### ExtFig9: