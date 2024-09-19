#!/bin/bash
#$ -l h_rt=72:0:0
#$ -l mem_free=16G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M Ian.Jones3@ucsf.edu

##################################################################
###SET THESE VARIABLES FOR EACH SAMEPLE
##################################################################
# please provide the path to the following tools:
source ~/.bashrc
conda activate MAPS_update

python_path=/shen/shenlabstore3/ijones1/dependencies/conda/envs/encode.test/envs/MAPS_update/bin/python #should have pysam, pybedtools installed. bedtools, samtools should be in the path
Rscript_path=/shen/shenlabstore3/ijones1/dependencies/conda/envs/encode.test/envs/MAPS_update/bin/Rscript

sample_name="IJ374_YSJJ07_100bp.fastp.norm"
outdir="/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/HPREP/2023_10mil_usable"

feather=0 #If 1 run entire pipeline, if 0 start from short/long bed/bedpe files

input_bam_file="/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/Down_sample_feather_HPrep_122221_bam/oligo/IJ351_YSXY14_100bp.fastp.9.5mil.paired.srtn.rmdup.bam"
genomic_feat_filepath="/shen/shenlabstore3/ijones1/dependencies/MAPS_update/MAPS-master/Arima_Genomics/genomic_features/hg38_F_GC_M_arima_5Kb_el.txt"
macs2_filepath="/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/Merged/Anchor.Merged.9.7.23.final.bed"

bin_size=5000
binning_range=1000000
filter_file="None"
mapq=30
chr_count=22
length_cutoff=1000
threads=1


##################################################################
###DO NOT CHANGE
##################################################################
fdr=2
model="pospoisson"
sex_chroms_to_process="NA"
optical_duplicate_distance=0


##################################################################
###SET THESE VARIABLES ONLY IF feather = 0
##################################################################
long_bedpe_dir="/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Norm/feather_output/IJ374_YSJJ07_100bp.fastp.norm_20230908_140532/"
short_bed_dir="/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Norm/feather_output/IJ374_YSJJ07_100bp.fastp.norm_20230908_140532/"
##################################################################


DATE=`date '+%Y%m%d_%H%M%S'`
feather_output=$outdir"/feather_output/"$sample_name"_"$DATE
resolution=$(bc <<< "$bin_size/1000")
per_chr='True' # set this to zero if you don't want per chromosome output bed and bedpe files
feather_logfile=$feather_output"/"$sample_name".feather.log"
# Had to change this
cwd="/shen/shenlabstore3/ijones1/dependencies/HPRep/bin"


HPRep_output=$outdir"/HPRep_output/"$sample_name"_"$DATE"/"

if [ $feather -eq 1 ]; then
	mkdir -p $feather_output
	$python_path $cwd/feather/feather_pipe preprocess -o $feather_output -p $sample_name -i $input_bam_file -q $mapq -l $length_cutoff -t $threads -c $per_chr -a $macs2_filepath -d $optical_duplicate_distance
	qc_filename=$feather_output/$sample_name".feather.qc"
	temp_qc_file=$feather_output/tempfiles/$sample_name".feather.qc.modified"
	sed -r 's/  +/\t/g' $qc_filename > $temp_qc_file
	sed -r -i 's/ /\_/g' $temp_qc_file
	cut -f 1-2 $temp_qc_file > $temp_qc_file".cut"
	sed -i 's/\_$//g' $temp_qc_file".cut"
	paste -d"\t" $cwd/feather/qc_template.txt $temp_qc_file".cut" > $temp_qc_file".cut.tmp"
	awk '{print $1,$3,$2}' $temp_qc_file".cut.tmp" >> $qc_filename".tsv"
	sed -i 's/ /\t/g' $qc_filename".tsv"
	cp "$(readlink -f $0)" $feather_output"/execution_script_copy"
	chmod 777 $feather_output
	mkdir -p $HPRep_output
	echo "$sample_name $HPRep_output $macs2_filepath $genomic_feat_filepath $feather_output"/" $feather_output"/" $bin_size $chr_count $maps_output $sex_chroms_to_process"
	$python_path $cwd/HPRep/make_and_xor_runfile.py $sample_name $HPRep_output $macs2_filepath $genomic_feat_filepath $feather_output"/" $feather_output"/" $bin_size $chr_count $HPRep_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	$python_path $cwd/HPRep/split_and_xor.py $HPRep_output"HPRep_"$sample_name".split"
        $Rscript_path $cwd/HPRep/find_anchor_bins.R $macs2_filepath $bin_size $chr_count ${HPRep_output}$sample_name"."$resolution"k.anchors"
        $Rscript_path $cwd/HPRep/normalize_bins.R $HPRep_output reg_raw $sample_name"."$resolution"k" $bin_size $chr_count "None" "pospoisson" $HPRep_output
fi

	
if [ $feather -eq 0 ]; then
	mkdir -p $HPRep_output
	echo "$sample_name $HPRep_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process"
	$python_path $cwd/HPRep/make_and_xor_runfile.py $sample_name $HPRep_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $HPRep_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	$python_path $cwd/HPRep/split_and_xor.py $HPRep_output"HPRep_"$sample_name".split"
	$Rscript_path $cwd/HPRep/find_anchor_bins.R $macs2_filepath $bin_size $chr_count ${HPRep_output}$sample_name"."$resolution"k.anchors"
	$Rscript_path $cwd/HPRep/normalize_bins.R $HPRep_output reg_raw $sample_name"."$resolution"k" $bin_size $chr_count "None" "pospoisson" $HPRep_output 
fi


