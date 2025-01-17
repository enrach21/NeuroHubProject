source ~/.bashrc
conda activate seqkit

cat  OPC.90000.peaks.fa | seqkit grep -s -i -p n

#### Had to manual clean up the bed file ####

cat OPC.90000.peaks.fa | seqkit grep -s -i -v -p n > OPC.90000.peaks.filt.fa

# Each fasta is 18 lines
head OPC.90000.peaks.filt.fa -n 1350000 > OPC.75k.peaks.filt.fa


grep -A 17 chr2_ OPC.75k.peaks.filt.fa | sed '/^--/d' > OPC.chr2.75k.peaks.filt.fa 

grep -A 17 -e chr1 -e chr3_ -e chr4_ -e chr5_ -e chr6_ -e chr7_ -e chr8_ -e chr9_  -e chr20_ -e chr21_ -e chr22_ -e chrX_ -e chrY_ OPC.75k.peaks.filt.fa | sed '/^--/d'   > OPC.train.75k.peaks.filt.fa 
