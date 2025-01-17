source ~/.bashrc
conda activate seqkit

cat MG.neg.peaks.fa | seqkit grep -s -i -p n

#### Had to manual clean up the bed file ####

cat MG.neg.peaks.fa | seqkit grep -s -i -v -p n > MG.neg.peaks.filt.fa



# Each fasta is 18 lines
head MG.neg.peaks.filt.fa -n 1350000 > MG.neg.peaks.filt.75k.fa


grep -A 17 chr2_ MG.neg.peaks.filt.75k.fa | sed '/^--/d' > MG.chr2.75k.neg.filt.fa 

grep -A 17 -e chr1 -e chr3_ -e chr4_ -e chr5_ -e chr6_ -e chr7_ -e chr8_ -e chr9_  -e chr20_ -e chr21_ -e chr22_ -e chrX_ -e chrY_ MG.neg.peaks.filt.75k.fa | sed '/^--/d'   > MG.train.75k.neg.filt.fa 