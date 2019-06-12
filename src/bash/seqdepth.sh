#!/bin/sh
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh seqdepth.sh folder dataset
#
# PARAMETERS:
# folder   folder of your data
# scratch   scratch folder
# dataset  name of the biocondition
#
###################################################################################
set -e

#  Setup true or false each step of the workflow to run it or not
RAW_STEP=false
FASTQ_STEP=false
BAM_STEP=true
WIG_STEP=false

# Setup variables
folder=$1
scratch=$2
data=$3
rawdata=$4
#module add samtools
###################################################################################
if $RAW_STEP; then
	sizeGZ=$(ls -al -h $folder/RNASeq_raw/all_gz/${rawdata}.fastq.gz)
	sizeFQ=$(ls -al -h $scratch/RNASeq_raw/${data}.fastq)
fi

if $FASTQ_STEP; then
	# Print the number of raw reads in the Fastq files
	LINE_COUNT=$(wc -l $scratch/RNASeq_raw/${data}.fastq)             #Counts the number of lines in fastq file
	PART=(${LINE_COUNT//\t/ })                          # Removes the file name from the output and only returns number of lines
	rawRead=$((PART[0]/4))
	echo $data\\t$READ_COUNT >> $folder/Seqdepth/NbRawReads.txt

	# Print the number of reads in the Fastq files
	LINE_COUNT=$(wc -l $scratch/FastQ/${data}.fastq)             #Counts the number of lines in fastq file
	PART=(${LINE_COUNT//\t/ })                          # Removes the file name from the output and only returns number of lines
	trimRead=$((PART[0]/4))
	echo $data\\t$READ_COUNT >> $folder/Seqdepth/NbReadsTrimmed.txt
fi

if $BAM_STEP; then
	# calculate sequencing depth
	if [ -e $folder/Mapping/${data}.bam ]; then
	    numberReads=$(samtools view -c $folder/Mapping/${data}.bam)
	    readsLength=$(samtools view $folder/Mapping/${data}.bam |  awk '{sum+=length($10)} END { print sum/NR}')
	    sizeTranscriptome=282000000
	    sequencingDepth=$(echo 'scale=5; '$numberReads' * '$readsLength' / '$sizeTranscriptome| bc)
	    echo $data\\t$numberReads >> $folder/Seqdepth/STAR_nbReads.txt
	fi
fi

if $WIG_STEP; then
	sizewig=$(ls -al -h /pasteur/projets/policy01/m6aAkker/PeakDetection/RawCov/${data}.wig)
fi

echo $data $rawdata $sizeGZ $sizeFQ $rawRead $trimRead $numberReads $sequencingDepth $sizewig >> $folder/Seqdepth/Sequencing_Summary.txt


