#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-581

# Activate deeptools
conda activate MSA


index='index.txt'
sample=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $index)
echo ${sample}

pos='pos.txt'
num=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $pos)
echo ${num}

chrom='chr.txt'
chr=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $chrom)
echo ${chr}


msa_view --in-format MAF ${sample}.maf | awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' | grep -w 'hg38\|panTro4' -A 1 > ${sample}.msa

# Add the HAR to the fasta name
awk -v env_var=$sample '/^>/{print $0 "_" env_var; next}{print}' < ${sample}.msa | tr -d '-' > ${sample}.V2.msa 


conda activate snp_sites

snp-sites -v -o ${sample}.vcf ${sample}.msa

# rm ${sample}.msa 

#### Testing ### 
seq=$(head ${sample}.msa -n 2 | tail -n 1)
end=$(wc -l < $sample.vcf)

for ((i=5;i<=end;i++)); do
    echo $i
    # val is the current length in the vcf file
    val=$(cut -f 2 ${sample}.vcf | head -n $i | tail -n 1)
    echo $val
    # New length to be put in the vcf files
    val2=$(echo ${seq:0:$val} | tr -d '-' | wc -c)
    echo $val2
    
    awk -vOFS='\t' -v new=$val2 -v line=$i 'NR==line {$2=new}1' < $sample.vcf > $sample.temp.vcf && mv $sample.temp.vcf $sample.vcf
done

# Make vcf file start from chromosome position rather than the one location
awk -v val2=$chr -v val=$num -v env_var="$sample" '{if (NR>4) {$1 = val2}  if (NR>4) {$2 = $2 + val - 1} if (NR>4) {$3 = env_var "_"  (NR-4) }}1' OFS='\t' < ${sample}.vcf > ${sample}.V2.vcf

rm $sample.temp.vcf
rm $sample.vcf
rm ${sample}.msa


#### NEED TO RUN THIS MANUALLY afterwards
# for f in *.V2.msa; do head -n 2 "$f" >> human.fa ; done


#### NEED TO RUN THIS MANUALLY afterwards
# for f in *.V2.msa; do tail -n 2 "$f" >> chimp.fa ; done



#### NEED TO RUN THIS MANUALLY afterwards
# for f in *.V2.vcf; do tail -n +5 "$f" >> HAR.VAR.vcf; done