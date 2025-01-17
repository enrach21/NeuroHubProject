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

if [ $CHH ]; then 
    zcat vRG.CHH_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk '{ if (($4) != 0) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  vRG_CHH.report.10.cov.bed

    Rscript Summary.R vRG_CHH.report.10.cov.bed  vRG_CHH.report.10.cov.csv
    fi
    
if [ $CHG ]; then 
    zcat vRG.CHG_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk '{ if (($4) != 0) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  vRG_CHG.report.10.cov.bed

    Rscript Summary.R vRG_CHG.report.10.cov.bed  vRG_CHG.report.10.cov.csv
    fi

if [ $CG ]; then 
    zcat vRG.CG_report.txt.gz | awk '{ if (($4 + $5) >= 10) { print } }' | awk 'OFS="\t" {$1="chr"$1; print}' | awk '{$8 = $4 / ($4 + $5)}1'  >  vRG_CG.report.10.cov.bed

    Rscript Summary.R vRG_CG.report.10.cov.bed  vRG_CG.report.10.cov.csv
    fi

