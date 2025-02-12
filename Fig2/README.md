# Fig2: cCREsLMAR  are associated with transcriptional regulation and enhancer activity.
Scripts involved in the creation of figure two.

### Fig2a: 
- Count total number of XOR and AND PLAC-seq interactions
    - PLAC-seq processing can be found [here](../DataProcessing/PLAC-seq)
- Files:
    - PlotNumInt.r: Rscript that reads in sig interactions and plots bar plot
### Fig2b:
- Obtain number of XOR interactions with distal cCRE^LMAR
    - LMARs defined in [ExtFig2a](../ExtFig/ExtFig2/ExtFig2a)
- Files:
    - InteractionWithLMAR.r: Rscript that overlaps distal bins within XOR interations with cCREs
### Fig2c:
- Comparing contact frequency of XOR interaction with and without distal cCRE^LMAR
    - Expanded analysis to all cell types found in [ExtFig3d](../ExtFig/ExtFig3/ExtFig3d)
- Files:
    - InteractionStrengthWithLMAR.r: Rscript comparing interaction strength between cCRE^LMAR and others
### Fig2d:
- Heatmap of unique significant XOR interaction overlapping to a single cell type with distal cCRE^LMAR
    - Unique interactions shown in [ExtFig4a](../ExtFig/ExtFig4/ExtFig4a)
- Files:
    - Get.Counts.11.21.23.sh: Obtains the ATAC-seq average signal within cCRE^LMAR
    - ATAC_seq_methylation.sh: Obtains the average methylation level within cCRE^LMAR 
    - HeatMap.R: Plots all signals
### Fig2ef:
- Obtains correlation between interactions and expression changes
    - All combindation shown in [ExtFig4b](../ExtFig/ExtFig4/ExtFig4b)
- Files:
    - CorPlacRNA.R: Rscript to obtain correlation
        - The function Count_TSS_ATAC.V3.R gets count differences of intereactions with distal cCRE per TSS
        - The function Count_TSS.R gets count differences of all intereactions per TSS
        - The function Plot_Corr.NH.R gets expression difference and correlation
### Fig2g:
- Plots the RG vista results and how they overlap with epigenomic data
- Files:
    - RG.Vista.R: Plot Summary of ATAC, WGBS, and PLAC at tested regions
### Fig2h:
- Fisher test comparing different classification of cCRE
    - Results stored in .csv file.
- Files:
    - Enrichment.R: Plot forest plots of fisher exact test results
