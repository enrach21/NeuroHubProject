# Example of files on how WGBS was processed

## Steps
1. mapping
2. Create Bigwig
    - Methyl.9.19.23_Cell.sh: File to filter for only CpG reads with greater than 10 reads
    - Maker_Browser.sh: Converts CpG with to bed format for washu browser
    - Make_Bigwig.sh: Converts CpG to bigwig files using in ExtFig2c