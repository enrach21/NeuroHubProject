# Fig4: Prioritizing neuropsychiatric disorder variants within the developing brain.
Scripts involved in the creation of figure four.

### Fig4a: 
- LDSC analysis of cCRE^LNAR
- Files:
    - Get_LDSC_Bed.R: Rscript to get bed files for input of LDSC
    - LDSC_figure.csv: LDSC results
    - Plot_LDSC.R: Rscript to plot LDSC analysis
### Fig4c:
- For each cell type trained a model and evalute variants...
    - Train:
        - Train_kmer.sh: Obtained score for 11-mers for delta-SVM
        - Training_Full.sh: Training used for final model
        - Training.sh: Training used to evaluate chr 2
        - More info found in [ExtFig7b](../ExtFig/ExtFig7/ExtFig7b)
    - AD:
        - DeltavSVM_shuf.sh: Delta SVM scores for shuf controls
        - DeltaSVM.sh: Delta SVM scores for variants
        - fasta.txt: fasta files input into scripts
        - Get_Sig_SNP.py: Script to summaries significant variants
        - output.txt: Name of output files
        - Score.ISM.sh: Script to obtain ISM scores
        - Score.sh: Script to obtain GKMexplain scores
    - SCZ:
        - DeltavSVM_shuf.sh: Delta SVM scores for shuf controls
        - DeltaSVM.sh: Delta SVM scores for variants
        - fasta.txt: fasta files input into scripts
        - Get_Sig_SNP.py: Script to summaries significant variants
        - output.txt: Name of output files
        - Score.ISM.sh: Script to obtain ISM scores
        - Score.sh: Script to obtain GKMexplain scores
- Files:
    - Plot_rs4449074.py: Script to visulaize rs4449074 GKMexplain results
        