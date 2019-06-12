#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh run.sh folder data annotation_file genome_file
# 
# PARAMETERS:
# folder   folder of your data
# data    name of the dataset
# annotation_file   Name of the annotation file
# genome_file    Name of the genome used
###################################################################################
set -e

#  Setup true or false each step of the workflow to run it or not
GZIP_STEP=true
TRIM_STEP=true
FASTQC_STEP=true
STAR_STEP=true
BAM_STEP=true
MOVE_STEP=true
WIG_STEP=true

# Setup variables
folder=$1
data=$2
annotation_file=$3
genome_file=$4

###################################################################################

# Test if fastq files exist
# Your sequencing data should be in : $folder/RNASeq_raw/
FILE=$folder'/RNASeq_raw/'$data'.fastq.gz'
if ! [ -f $FILE ]; then
	echo "ERROR - File $FILE does not exist"
   	exit 1 
fi

# GZip - Extract fastq
if $GZIP_STEP; then
	echo "Extract archive " $data'.fastq.gz'
	gzip -k -d $FILE
fi

# AlienTrimmer - Reads Trimming
if $TRIM_STEP; then
	adapter_fasta_file=$folder/Genome/alienTrimmerPF8contaminants.fasta
	AlienTrimmer_Soft="java -jar /opt/AlienTrimmer_0.4.0/src/AlienTrimmer.jar"
	echo "Trim reads in " $data'.fastq' " with " $adapter_fasta_file
	$AlienTrimmer_Soft -i $folder/RNASeq_raw/$data.fastq -c $adapter_fasta_file -o $folder/FastQ/$data.fastq
fi

# FastQC - Control Reads quality
if $FASTQC_STEP; then
	echo "Run FastQC for " $data'.fastq'
	fastqc -o $folder/FastQC/ $folder/FastQ/$data.fastq
fi

if $STAR_STEP; then
	echo "Map reads from:" $data'.fastq'
	# prepare variables
	genome=GRCm38.primary_assembly.genome.fa
	genomeFolder=$folder/Genome
	annotation=$genomeFolder/$annotation_file.annotation.gtf
	mappingFolder=$folder/Mapping/STAR_results/${data}
	if [ ! -d $mappingFolder ]; then
		mkdir $mappingFolder
	fi
	readFile=$folder/FastQ/${data}.fastq
	chromSize=$genomeFolder/chrNameLength.txt
	wigfile=$mappingFolder/STAR_results/${data}/${data}Signal.Unique.str1.out.wig
	bwfile=$mappingFolder/${data}.bw

	# run STAR mapper
	STAR --genomeDir ${genomeFolder} --outFileNamePrefix ${mappingFolder}/${dataName} \
	--sjdbGTFfile ${annotation} --sjdbOverhang 100 --readFilesIn ${readFile} --outSAMtype BAM SortedByCoordinate --outWigType wiggle \
	--outWigStrand Unstranded --outWigNorm RPM
	# # create bigwig
	echo "CreateBigWig file" $bwfile
	#/pasteur/homes/cbecavin/opt/wigtools/
	module load  wigToBigWig
	wigToBigWig ${wigfile} ${chromSize} ${bwfile}
fi

if $BAM_STEP; then
	# prepare variables
	mappingFolder=$folder/Mapping
	bamFileUnfiltered=${mappingFolder}/STAR_results/${dataName}/${dataName}Aligned.sortedByCoord.out.bam
	bamFileFiltered=${mappingFolder}/${data}
	echo "Filter and index BAM file" $bamFileFiltered
	# filter bam files by removing data with score = 0
	samtools view -b -q 1 ${bamFileUnfiltered} > ${bamFileFiltered}.bam
	#samtools sort ${bamFileFiltered}.bam ${bamFileFiltered}
	samtools index ${bamFileFiltered}.bam
fi

if $WIG_STEP; then
	echo "Calculate raw coverage wig file for POI peak detection technique : " $data
	bamFile=${mappingFolder}/${data}.bam
	wig_data=$folder/PeakDetection/RawCov/$data.wig
  	bedtools genomecov -d -split -ibam $bamFile | awk '$3!=0 {print $0}' > $wig_data
fi