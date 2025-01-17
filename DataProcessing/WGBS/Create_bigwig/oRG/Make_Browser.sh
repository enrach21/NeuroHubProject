#/bin/bash
#$ -l h_rt=4:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=200G
#$ -j y
#$ -r y

awk '$9=$2' oRG_CG.report.10.cov.bed | awk '{print $1,$2,$9,$6,$8,$3,($4+$5)}' | awk '{if ( $4 ~ /+/ ) {print $1,$2-1,$3,$4,$5,$6,$7} else {print $1,$2-1,$3,$4,$5,$6,$7}}' | sort -k1,1 -k2,2n | awk -v OFS='\t' '{ $1=$1; print}' > oRG_CG.report.10.cov.browser.bed

# bgzip vRG_CG.report.10.cov.browser.bed
# tabix -p bed vRG_CG.report.10.cov.browser.bed.gz