# Fig1: Collecting cell types and annotating cCREs within the developing cortex 
Scripts involved in the creation of figure one.

### Fig1c: 
- Obtaining TMM-normalized RPKM and expression of Marker Genes
    - RNA-seq processing can be found [here](../DataProcessing/RNA-seq)
- Files:
    - Fig1C.sh: Bash script to create figre
    - Fig1C.R: Rscript that plots TMM-RPKM of marker genes    
### Fig1d:
- Venn of cCRE LMARs
    - ATAC-seq processing can be found [here](../DataProcessing/ATAC-seq)
    - WGBS Processing can be found [here](../DataProcessing/WGBS) 
    - LMARs defined in [ExtFig2a](../ExtFig/ExtFig2/ExtFig2a)
- Files:
    - Intervene.Both.sh: Bash script using [Intervene](https://intervene.readthedocs.io/en/latest/modules.html) pacakge 
### Fig1e:
- [HOMER](http://homer.ucsd.edu/homer/motif/) enrichment of cCRE LMARs for each cell type
- Files:
    - vRG_Homer.sh: enrich within cCREs^LMAR found in exclusivley vRG
    - oRG_Homer.sh: enrich within cCREs^LMAR f found in exclusivley oRG
    - OPC_Homer.sh: enrich within cCREs^LMAR f found in exclusivley OPC
    - MG_Homer.sh: enrich within cCREs^LMAR f found in exclusivley MG
    - PlotHomer.R: Rscript to plot p.values and expression of TF found to be enriched
