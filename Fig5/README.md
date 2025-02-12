# Fig5: oRG cCREs are enriched for HARs
Scripts involved in the creation of figure five.

### Fig5a: 
- Random sampling of cCREs to determine HAR enrichment. Cell type sepecific cCREs obtain from [Fig1d](../Fig1/Fig1d)
- Files:
    - HAR_Sampling.R: Rscript to perfrom random sampling
### Fig5b:
- Investigated the association of HARs with vRG or oRG [DARs](../DataProcessing/ATAC-seq/DARs)
- Files:
    - Vista Enrichment - HAR_DAR.csv: Output of script
    - OR_DAR_HAR.R: Rscript to perform fish exact test on DARs overlap with HARs
### Fig5c:
- Use GKM to prioritze HAR variants predicted to distrupt accessiblity and subsequently obtained GO-terms of genes interacting with HARs
- Files:
    - GKM_Predictions: Code to get fasta file with variants and predict accessiblity changes
        1. Get bed file of accessible HARs 
        2. get_MAF.sh: obtain MAFs in all accessible HARs
        3. get_fasta.sh: Runs MSA view to obtain fasta file
        4. Shuffle.sh: Get Shuffled fasta files for controls
    - Gene_Ontology.R: Rscript to perform GeneOntology
    - GKM_Sig_Target_genes.txt: List of target Genes of perturbed HARs
    - oRG.GKM.HAR.Targets.go.csv: GeneOntology Output
### Fig5de:
- Browser showing predicted accessibility difference of humand and chimp sequence for HARsv2_2157
- Files:
    - Visualize_HAR_2157.py: Python script to visualize HARsv2_2157
### Fig5f:
- Predict disruption of ZEB1 motif by HARsv2_2157
- Files:
### Fig5g:
- Luciferase results
- Files:
    