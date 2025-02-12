#/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 2-17
#$ -m ea
#$ -M Ian.Jones3@ucsf.edu


# activate enviroment:
source ~/.bashrc
conda activate base
date
hostname



# Define folders and load dependencies.
# Dir where fastq files are located
# LANE=/shen/shenlabstore/sequencing_data/YSJJ07/YSJJ07
# Sample name
INDEX=/shen/shenlabstore3/ijones1/NeuroHub/RNA-seq/Scripts/index.txt
SAMPLE=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX)

INDEX2=/shen/shenlabstore3/ijones1/NeuroHub/RNA-seq/Scripts/index2.txt
LANE=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX2)
# Type of sequencing
format=pe
# Dir where files will be output
RESULTS=/wynton/group/shen/NeuroHub.7.23.21/RNA-seq/STAR_RSEM_9_7_23
# Length fastq files will be trimmed to
trim=100



# activate enviroment:
source ~/.bashrc
conda activate rna_analysis

# Make directory for this sample
mkdir $RESULTS/$SAMPLE

# Make dir for fastp output
mkdir $RESULTS/$SAMPLE/fastp

# Read in input files
# Trim to 100bp
# Trim Adapters
/shen/shenlabstore3/ijones1/dependencies/fastp/fastp -i $LANE/${SAMPLE}*R1*fastq.gz -I $LANE/${SAMPLE}*R2*fastq.gz  -o $RESULTS/$SAMPLE/fastp/${SAMPLE}.${trim}bp.fastp.R1.fq.gz -O $RESULTS/$SAMPLE/fastp/${SAMPLE}.${trim}bp.fastp.R2.fq.gz -b ${trim} -B ${trim} -p -h "$RESULTS/$SAMPLE/fastp/${SAMPLE}.fastp.html" -j "$RESULTS/$SAMPLE/fastp/${SAMPLE}.fastp.json"

# Print out user input.
echo $format
which STAR
which rsem-calculate-expression


READ1=$RESULTS/$SAMPLE/fastp/${SAMPLE}*fastp*R1*f*.gz 
READ2=$RESULTS/$SAMPLE/fastp/${SAMPLE}*fastp.R2*f*.gz 
echo $READ1
echo $READ2

# Run analysis using STAR/RSEM.
conda activate RNAencode
cd $RESULTS/$SAMPLE
if [ $format == se ]; then 
	
	echo performing single-end analysis
	
	STAR --genomeDir /shen/shenlabstore3/ijones1/dependencies/ATAC_genomes/RSEM/ --readFilesCommand zcat \
        --readFilesIn $READ1 \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
        --outFilterMultimapNmax 20 --quantMode TranscriptomeSAM \
        --runThreadN 16 --outFileNamePrefix ${SAMPLE}. \
	--outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --genomeLoad LoadAndKeep

	rsem-calculate-expression --bam --no-bam-output -p 8 \
        --forward-prob 0 --estimate-rspd --calc-ci \
        ${SAMPLE}.Aligned.toTranscriptome.out.bam \
        shen/shenlabstore3/ijones1/dependencies/ATAC_genomes/RSEM/GRCH38.p13_gencode.v38 $SAMPLE >& ${SAMPLE}.rsem.log	

fi
if [ $format == pe ]; then
 
	echo performing paired-end analysis
	
 	STAR --genomeDir /shen/shenlabstore3/ijones1/dependencies/ATAC_genomes/RSEM/update --readFilesCommand zcat \
	--readFilesIn $READ1 $READ2 \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
	--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM \
	--runThreadN 16 --outFileNamePrefix ${SAMPLE}. \
	--outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --genomeLoad LoadAndKeep

	rsem-calculate-expression --bam --no-bam-output -p 8 \
	--paired-end --forward-prob 0 --estimate-rspd --calc-ci \
	${SAMPLE}.Aligned.toTranscriptome.out.bam \
	/shen/shenlabstore3/ijones1/dependencies/ATAC_genomes/RSEM/update/GRCH38.p13_gencode.v38.main $SAMPLE >& ${SAMPLE}.rsem.log

 fi


# Generate bedGraph files from BAM using STAR.
# Activate deeptools
conda activate deeptools
# make a new dir
mkdir visualization

# Read in Mapped Bed files
bam=$RESULTS/$SAMPLE/${SAMPLE}*.sortedByCoord.out.bam 
# Sort bam file
samtools sort ${bam} -o $RESULTS/$SAMPLE/${SAMPLE}.temp.sorted.bam
# Index the bam files
samtools index $RESULTS/$SAMPLE/${SAMPLE}.temp.sorted.bam

# Make bigwig
bamCoverage -b $RESULTS/$SAMPLE/${SAMPLE}.temp.sorted.bam --filterRNAstrand forward --normalizeUsing CPM -o $RESULTS/$SAMPLE/visualization/${SAMPLE}_hg38_cpm_forward.bw
bamCoverage -b $RESULTS/$SAMPLE/${SAMPLE}.temp.sorted.bam --filterRNAstrand reverse --normalizeUsing CPM -o $RESULTS/$SAMPLE/visualization/${SAMPLE}_hg38_cpm_reverse.bw
