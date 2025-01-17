# Given MACs file, sort by the MACs signal and then select he top 50k peaks. Afterwards 


#### Need to pick 50001 since 1 read had ambgiuous N ##### 

zcat MG.overlap.optimal_peak.narrowPeak.gz | sort -k5 | head -n 90000 | awk -v FS='\t' -v OFS='\t' -v OFMT='%f' '{print $1, ($2+$10-500), ($2+$10+500), ("peak_"NR) }' > MG.90000.peaks.bed
