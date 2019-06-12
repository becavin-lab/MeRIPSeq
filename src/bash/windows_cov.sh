#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh windows_cov.sh folder annotation_file data peak_technique folder_script
#
# PARAMETERS:
# folder   folder of your data
# annotation_file   Name of the annotation file
# genome_file    Name of the genome used
# peak_technique   name of the technique to run
# folder_script   folder where to find the python scripts
#
###################################################################################

# Setup variables
folder=$1
annotation_file=$2
data=$3
peak_technique=$4
folder_script=$5

echo $peak_technique
###################################################################################

if [[ " ${peak_technique} " == *" POI "* ]]; then
	module add bedtools
	source activate py27
  	sliding_annot=${folder}/Genome/$annotation_file.exon.slidingwindows.gtf
	gene_annot=${folder}/Genome/$annotation_file.annotation.exon.gtf
	wig_data=$folder/PeakDetection/RawCov/$data.wig
	echo "Calculate median window coverage for POI peak detection technique : " $data
   	output=${folder}/PeakDetection/RawCov/HTSeq_Median_${data}_windows.txt
   	python $folder_script/src/python/htseq_window.py -a $sliding_annot -i $wig_data -o $output
	echo "Calculate median transcript coverage for POI peak detection technique : " $data
	output=${folder}/PeakDetection/RawCov/HTSeq_Median_${data}_genes.txt
 	python $folder_script/src/python/htseq_transcript.py -a $gene_annot -i $wig_data -o $output
fi

if [[ (" ${peak_technique} " == *" Fisher "* ) || (" ${peak_technique} " == *" RPMF "*)]]; then
	echo "Count number of reads in window for Fisher or RPMF peak detection technique : " $data
	sliding_annot=${folder}/Genome/$annotation_file.gene.slidingwindows.gtf
	output=${folder}/PeakDetection/RawCov/HTSeq_Count_${data}_windows.txt
	htseq-count -i ID -t window -s no -m union --nonunique all -f 'bam' "${folder}/Mapping/${data}.bam" "$sliding_annot" > "$output"
fi