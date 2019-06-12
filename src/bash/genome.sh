#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh genome.sh annotation_file genome_file
#
# PARAMETERS:
# folder   folder of your data
# annotation_file   Name of the annotation file
# genome_file    Name of the genome used
#
###################################################################################

#  Setup true or false each step of the workflow to run it or not
DOWNLOAD_ANNOTATION=true
DOWNLOAD_GENOME=true
STAR_GENOME=true
PARSE_ANNOTATION=true
BOWTIE_RRNA=true
CREATE_WINDOWS=true
SALMON_INDEX=true
DOWNLOAD_METDB=true

# Setup variables
genome_folder=$1/Genome
script_folder=$2
annotation_file=$3
genome_file=$4.fa
entrez_file=$3.metadata.EntrezGene
python --version

###################################################################################

if $DOWNLOAD_ANNOTATION; then
	echo "Download Genome annotation from Gencode website"
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M13/$annotation_file.annotation.gtf.gz
	gzip -d $annotation_file.annotation.gtf.gz
	mv $annotation_file.annotation.gtf $genome_folder/$annotation_file.annotation.gtf

	echo "Download Entrez information from Gencode website"
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M13/$entrez_file.gz
	gzip -d $entrez_file.gz
	mv $entrez_file $genome_folder/$entrez_file
fi

if $DOWNLOAD_GENOME; then
	echo "Download Genome sequence from Gencode website"
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M13/$genome_file.gz
	gzip -d $genome_file.gz
	mv $genome_file $genome_folder/$genome_file
fi

if $STAR_GENOME; then
	echo "Create STAR Genome for Mapping"
	STAR --runMode genomeGenerate --genomeDir $genome_folder --genomeFastaFiles $genome_folder/$genome_file
fi

if $PARSE_ANNOTATION; then
	echo "Create Annotation files"
	python $script_folder/src/python/parse_annotation.py -p $genome_folder/ -a $annotation_file
fi

if $BOWTIE_RRNA; then
	echo "Create rRNA files"
	python $script_folder/src/python/create_rRNA_genome.py -p $genome_folder/ -a $annotation_file -o bed
	sort -k1,1 -k2,2n $genome_folder/$annotation_file'.annotation.rRNA.bed' > $genome_folder/$annotation_file'.annotation.rRNA_sort.bed'
	bedtools merge -i $genome_folder/$annotation_file'.annotation.rRNA_sort.bed' -s > $genome_folder/$annotation_file'.annotation.rRNA_merge.bed'
	python $script_folder/src/python/create_rRNA_genome.py -p $genome_folder/ -a $annotation_file -g $genome_file -o fasta
	rm $genome_folder/$annotation_file'.annotation.rRNA.bed'
	rm $genome_folder/$annotation_file'.annotation.rRNA_sort.bed'
	rm $genome_folder/$annotation_file'.annotation.rRNA_merge.bed'
	echo "Build Bowtie database for rRNA"
	bowtie2-build $genome_folder/$annotation_file.rRNA.fasta $genome_folder/$annotation_file_rRNA
fi

if $SALMON_INDEX; then
	echo "Create index for Salmon search"
	# build index
	transcripts_fasta=$genomeFolder/$annotation_file.transcripts.fa
	fasta=$genome_folder/$annotation_file.annotation.transcript.gfa
	salmon index -t $fasta -i transcripts_index -k 31
fi

if $CREATE_WINDOWS; then
	echo "Create 100bp windows files"
	python $script_folder/src/python/create_windows.py -p $genome_folder/ -a $annotation_file
fi

if $DOWNLOAD_METDB; then
	echo "Download MeRIP-Seq data from MetDB-V2"
	wget http://180.208.58.19/metdb_v2/download/overall/pc_ep_mouse_mm10.gz
	gzip -d pc_ep_mouse_mm10.gz
	mv pc_ep_mouse_mm10 $genome_folder/pc_ep_mouse_mm10

	echo "Download TREW data from MetDB-V2"
	wget http://180.208.58.19/metdb_v2/download/overall/trew_mouse_mm10.gz
	gzip -d trew_mouse_mm10
	mv trew_mouse_mm10 $genome_folder/trew_mouse_mm10

	echo "Download CLIP-Seq data from MetDB-V2"
	wget http://180.208.58.19/metdb_v2/download/overall/sb_m6a_mouse_mm10.gz
	gzip -d sb_m6a_mouse_mm10.gz
	mv sb_m6a_mouse_mm10 $genome_folder/sb_m6a_mouse_mm10

	echo "Prepare .bed files for MeTDB-V2 files"
	python $script_folder/src/python/reference_peak.py -p "$genome_folder"
fi


