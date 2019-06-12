#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh peak_detection.sh folder data ip_data input_data annotation_file peak_technique
#
# PARAMETERS:
# folder   folder of your data
# data    name of the dataset
# ip_data    name of the corresponding IP (MeRIPSeq) dataset
# input_data    name of the corresponding Input (RNASeq) dataset
# peak_technique   name of peak detection technique to use
#
###################################################################################

#  Setup true or false each step of the workflow to run it or not
MACS2_STEP=true
FISHER_STEP=true
RPMF_STEP=true
POI_STEP=true

# Setup variables
folder=$1
biocond=$2
ip_data=$3
input_data=$4
peak_technique=$5
script_folder=$6

WINDOW_CUTOFF=10  # 10 reads minimum in windows


###################################################################################

if $MACS2_STEP; then
	if [[ " $peak_technique " == *" MACS2 "* ]]; then
		outfolder=$folder/PeakDetection/temp/
		module add macs/2.1.0
		echo "Run MACS2 peak detection for" $biocond
		macs2 callpeak -t $folder/Mapping/${ip_data}.bam -c $folder/Mapping/${input_data}.bam -n $biocond -f 'BAM' -g 282000000 --nomodel --outdir $outfolder
	fi
fi

if $POI_STEP; then
	if [[ " $peak_technique " == *" POI "* ]]; then
	  	python $script_folder/src/python/peak_calling_poi.py -p "$folder" -b $biocond -i $ip_data -c $input_data -w $WINDOW_CUTOFF
	fi
fi

if $RPMF_STEP; then
	if [[ " $peak_technique " == *" RPMF "* ]]; then
		python $script_folder/src/python/peak_calling_rpmf.py -p "$folder" -b $biocond -i $ip_data -c $input_data -w $WINDOW_CUTOFF
	fi
fi

if $FISHER_STEP; then
	if [[ " $peak_technique " == *" Fisher "* ]]; then
		python $script_folder/src/python/peak_calling_fisher.py -p "$folder" -b $biocond -i $ip_data -c $input_data -w $WINDOW_CUTOFF
	fi
fi
