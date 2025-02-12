# Fig3:  Identifying transcription factors driving epigenomic changes between vRG and oRG.
Scripts involved in the creation of figure three.

### Fig3a: 
- Volcano plot of DEG
    - DARs used in this analysis are called [here](../DataProcessing/ATAC-seq/DAR)
- Files:
    - LinkDARtoDEGs.R: Get the genes that are interacting or overlap with DARs
    - Volcano.R: Plot DEG with number of DARs in contacts
### Fig3b:
- Tobias Motif binding differences
- Files:
    - diffTFBS_oRG_vRG_bindetect_results.txt: Tobias results without expression filtering
    - diffTFBS_oRG_vRG_bindetect_results: Tobias results with motif required to have greater than 10 rpkm
    - vRG_oRG_Volcano_Tobias.R: Rscript of plotting binding differences
### Fig3c:
- Heat map of significant TFs matched with the change in expression
- Files:
    - Exp.Binding.Cor.R: Rscript to output heatmap.   
### Fig3d:
- LHX2 network created with cytoscpape highlighting DEGs regulated by LHX2 binding sites
- Files:
    - CytoScape.LHX2.sh: Overlap LHX2 bound motifs with PLAC-seq and TSS
    - CytoScape.LHX2.R: Rscript to combine gene targets with DEGs
    - LHX2.cytoscape.csv: Output of DEGs regulated by LHX2 specifically in oRG
### Fig3e:
- Expression of all targets of the 126 LHX2 motifs predcited to be more bound in oRG
- Files:
    - LHX2_Target_Expression.R: Rscript comparing all gene targets of the 126 LHX2 motifs