#/bin/bash
#$ -l h_rt=4:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=200G
#$ -j y
#$ -r y

# Author: Ian Jones
# Date: 7.12.23
# Determine methylation levels within 2kb



CHH=0
CHG=1
CG=1

# if [ $CHH ]; then 
#     zcat MG.CHH_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk '{ if (($4) != 0) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  MG_CHH.report.10.cov.bed

  #   Rscript Summary.R MG_CHH.report.10.cov.bed  MG_CHH.report.10.cov.csv
 #    fi
    
# if [ $CHG ]; then 
#    zcat MG.CHG_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk '{ if (($4) != 0) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  MG_CHG.report.10.cov.bed

 #    Rscript Summary.R MG_CHG.report.10.cov.bed  MG_CHG.report.10.cov.csv
 #    fi

if [ $CG ]; then 
    zcat MG.CG_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  MG_CG.report.10.cov.bed

    Rscript Summary.R MG_CG.report.10.cov.bed  MG_CG.report.10.cov.csv
    fi

