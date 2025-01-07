# Fig2: cCREsLMAR  are associated with transcriptional regulation and enhancer activity.
Scripts involved in the creation of figure two.

### Fig2a: 
- Count total number of XOR and AND PLAC-seq interactions
    - PLAC-seq processing can be found [here](../DataProcessing/PLAC-seq)
### Fig2b:
- Obtain number of XOR interactions with distal cCRE^LMAR
### Fig2c:
- Comparing contact frequency of XOR interaction with and without distal cCRE^LMAR
### Fig2d:
- Heatmap of unique significant XOR interaction overlapping to a single cell type with distal cCRE^LMAR
    - Get.Counts.11.21.23.sh obtains the ATAC-seq average signal within cCRE^LMAR
    - ATAC_seq_methylation.sh obtains the average methylation level within cCRE^LMAR 
### Fig2ef:
- Obtains correlation between interactions and expression changes
    - CorPlacRNA.R is Rscript to obtain correlation
        - The function Count_TSS_ATAC.V3.R gets count differences of intereactions with distal cCRE per TSS
        - The function Count_TSS.R gets count differences of all intereactions per TSS
        - The function Plot_Corr.NH.R gets expression difference and correlation
### Fig2g:
- Plots the RG vista results and how they overlap with epigenomic data
### Fig2h:
- Fisher test comparing different classification of cCRE
    - Results stored in .csv file.
