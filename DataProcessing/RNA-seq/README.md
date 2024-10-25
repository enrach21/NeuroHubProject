# Example of files on how RNA-seq was processed

## Steps
1. mapping
    - Fastp was used to trim and remove low quality reads
    - STAR mapped each reads
    - RSEM to quantified the gene expression