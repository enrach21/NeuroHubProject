source ~/.bashrc
conda activate seqkit

cat  oRG.90000.peaks.fa | seqkit grep -s -i -p n

#### Had to manual clean up the bed file ####

cat oRG.90000.peaks.fa | seqkit grep -s -i -v -p n > oRG.90000.peaks.filt.fa

# Each fasta is 18 lines
head oRG.90000.peaks.filt.fa -n 1350000 > oRG.75k.peaks.filt.fa


grep -A 17 chr2_ oRG.75k.peaks.filt.fa | sed '/^--/d' > oRG.chr2.75k.peaks.filt.fa 

grep -A 17 -e chr1 -e chr3_ -e chr4_ -e chr5_ -e chr6_ -e chr7_ -e chr8_ -e chr9_  -e chr20_ -e chr21_ -e chr22_ -e chrX_ -e chrY_ oRG.75k.peaks.filt.fa | sed '/^--/d'   > oRG.train.75k.peaks.filt.fa 



