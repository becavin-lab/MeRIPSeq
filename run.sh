#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
# GitLab: https://gitlab.pasteur.fr/cbecavin/MeRIPSeq/
#
# Executable (2 versions are available one for local computing, another for slurm cluster computing)
# For local analysis
# sh run.sh project_name
# For example: sh run.sh Liver
# For Slurm cluster
# nohup ./RunSlurmExample.sh > runServer.log 2>&1 &
# RunSlurmExample.sh should be executable
# chmod 777 RunSlurmExample.shFor Slurm cluster
#
# PARAMETERS:
# project-name   -  Name of the project
# The list of file should be in your_folder/Exp_Design/project_name_exp_design.txt
# Your sequencing data should be in : your_folder/RNASeq_raw/
# our experiment design file should be a tab-separated table with at least 5 columns:
# -----------------------------------------------------------------------
# DataName	IP_dataset	Input_dataset	BioCond	Seq
# Am_ZT3_1	Am_ZT3_IP_1	Am_ZT3_Input_1	Am	1
# Am_ZT13_1	Am_ZT13_IP_1	Am_ZT13_Input_1	Am	2
# CONV_ZT3_1	CONV_ZT3_IP_1	CONV_ZT3_Input_1	CONV	1
# CONV_ZT13_4	CONV_ZT13_IP_4	CONV_ZT13_Input_4	CONV	2
# -----------------------------------------------------------------------
# (See example/Liver_exp_design.txt) 
#
###################################################################################
# PIPELINE Description:
# Setup true or false each step of the workflow to run it or not
# INSTALL first all dependencies: See INSTALL.txt
# 1 - setup.sh - Prepare Analysis by creating all necessary folders
SETUP_STEP=true
# 2 - genome.sh - Prepare genome sequence, annotation, and windows for peak detection
GENOME_STEP=true
# Run on every dataset:
# WARNING - Steps 3 to 7 are run separately for each datasets.
# WARNING - They might take a ot of time to perform, so they can be parallelized on a cluster.
# 3 - Trimming, mapping, BAM file filtering and quality control, And WIG calculation for peak detection
MAPPING_STEP=true
# 4 - Calculate SeqDepth for peak detection techniques
SEQDEPTH_STEP=false
# 5 - HTSeq - Count number of reads per genes
HTSEQ_STEP=false
# 6 - coverage_window.py - Count number of reads per window
WINDOW_STEP=false
# 7 - peak_detection.py - Run Fisher, POI,and RPMF peak detection techniques
# 8 - MACS2 - Run MACS2 peak detection
PEAK_DETECTION_STEP=false
# Run only one time per project_name 
# 9 - Search for median coverage for max coverage search
MEDIAN_COV_STEP=false
# 10 - MultiQC for quality control
MULTIQC_STEP=false
# 11 - finalize peak detection by filtering out "bad" peaks, searchning for max coverage and peak center
PEAK_FINALIZE_STEP=false
# 12 - (Optionnal) Annotate all peaks from all different techniques
PEAK_ANNOTATE_STEP=false
# 13 - Search max coverage of every peaks and center peak position on it. 
SEARCH_MAX_STEP=false
# 14 - Regroup peaks from all the detection techniques, annotate them, find overlap position with genes and referent MeRIP-Seq, CLIPSeq and TREW
REGROUP_PEAK_STEP=false
# Run on every dataset:
# 15 - HTSeq - Count number of reads in each methylation sites detected
HTSEQ_PEAK_STEP=false
# Run only one time per project_name 
# 16 - Perform differential methylation sites analysis in R
DIFF_PEAK_STEP=false
# 17 - Run differential splicing events
SPLICING_STEP=false
# 18 - Run Motif presence algorithm
MOTIF_STEP=false
# 19 - Run GuitarPlot for every list of peak
GUITAR_STEP=false
# 20 - (Optional) Compare two peak lists to Ref peaks
VENN_STEP=false
#
###################################################################################
set -e

# Setup variables
# folder of all datasets
folder=/Users/cbecavin/Documents/m6aAkker/m6aAnalysis
# project_name given in argument
project_name=$1
# exp_design_file with datasets and model paremters description
exp_design_file=$folder/Exp_Design/$project_name'_exp_design.txt'
# annotation file for your genome
annotation_file=gencode.vM13
# genome name
genome_file=GRCm38.primary_assembly.genome
# Choose peak detection technique to run
# defaut : peak_detections=("MACS2" "RPMF" "POI" "Fisher")
peak_detections=("MACS2" "RPMF" "POI" "Fisher")
# Peak size to keep after max coverage detection
peak_size=150
# min number of peak detection technique by which a peak have to be detected to be included in peak diff expression analysis 
min_number_technique=3

# different bed_name to use for peak detection
bed_name_raw=Raw
bed_name_max=MaxValues
bed_name_regroup=${bed_name_max}_${min_number_technique}
final_bed_name=MaxMaxValues_${min_number_technique}

###################################################################################
# Run the workflow
echo ""
echo "####################################  MERIP-SEQ WORKFLOW  ###############################################"
echo "Run MeRIP-Seq analysis"
echo "Author: Christophe Bécavin, Anne Biton, Hub Bioinformatic and Biostatistic, Institut Pasteur"
# set the directory for your data
echo "Perform analysis in " $folder
echo "Exprimental design file : " $folder'/Exp_Design/'$project_name'_exp_design.txt'
echo "Your sequencing data should be in :" $folder'/RNASeq_raw/'
echo "Using genome annotation: " $annotation_file
echo "And genome sequence: " $genome_file
echo "Run peak detection using" ${#peak_detections[@]} "techniques:" ${peak_detections[@]}

echo "##############################  READ EXPERIMENTAL DESIGN  #########################################"
source activate py27  # Activate python version 2.7
DATANAMES=()
IP_DATA=()
INPUT_DATA=()
if [ ! -f $exp_design_file ]; then
	echo "ERROR - Cannot find Experimental design file: " $exp_design_file
	exit 1
else
	echo "Read Experimental design file: " $exp_design_file
	index=0
	while read p; do
		dataname=$(echo "$p" | awk '{print $1}')
		ip_data=$(echo "$p" | awk '{print $2}')
		input_data=$(echo "$p" | awk '{print $3}')
		if [[ "$dataname" != "BioCond" ]]; then
			echo $dataname;
			DATANAMES[index]=$dataname
			IP_DATA[index]=$ip_data
			INPUT_DATA[index]=$input_data
			let "index += 1"
		fi
	done < $exp_design_file
fi

echo "Found" ${#DATANAMES[@]} "biological conditions"
echo ${#IP_DATA[@]} "MeRIP-SEQ (IP) datasets and" ${#INPUT_DATA[@]} "RNASeq (Input) datasets"
echo "##############################  RUN BIOINFORMATIC PIPELINE  #########################################"
echo "##############################  SETUP STEP   ######################################"  ##
# STEP : setup.sh - Prepare Analysis
if $SETUP_STEP; then
	sh src/bash/setup.sh $folder
fi

echo "##############################  GENOME STEP   ######################################"  ##
# STEP : genome.sh - Prepare genome sequence, annotation, and windows for peak detection
if $GENOME_STEP; then
	sh src/bash/genome.sh $folder $annotation_file $genome_file
fi

echo "##############################  MAPPING STEP  ######################################"  ##
# STEP : Trimming, mapping, Quality Control
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $MAPPING_STEP; then
	for index in ${!DATANAMES[@]}; do
		ip_data=${IP_DATA[index]}
		input_data=${INPUT_DATA[index]}
		echo "Trimming, mapping, Quality Control of " $ip_data'.fastq.gz'
		sh src/bash/mapping.sh $folder $ip_data $annotation_file $genome_file
		echo "Trimming, mapping, Quality Control of " $input_data'.fastq.gz'
		sh src/bash/mapping.sh $folder $input_data $annotation_file $genome_file
	done
fi

echo "##############################  SEQ DEPTH STEP  ######################################"  ##
# STEP : SeqDepth - Control number of reads and calculate sequencing depth of each dataset (Important for Fisher test peak detection)
if $SEQDEPTH_STEP; then
	for index in ${!DATANAMES[@]}; do
		ip_data=${IP_DATA[index]}
		echo "Run Sequencing depth calculation for" $ip_data
		sh src/bash/seqdepth.sh $folder $ip_data
		input_data=${INPUT_DATA[index]}
		echo "Run Sequencing depth calculation for" $input_data
		sh src/bash/seqdepth.sh $folder $input_data		
	done
fi

echo "##############################  HTSeq STEP  ######################################"  ##
# STEP : HTSeq - Count number of reads per genes
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $HTSEQ_STEP; then
	for index in ${!DATANAMES[@]}; do
		input_data=${INPUT_DATA[index]}
		genome_annot=${folder}/Genome/$annotation_file.annotation.exon.gtf
		output=${folder}/Expression/HTSeq_${input_data}.txt
		echo "Run HTSeq count of genes for" $input_data
		htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' "${folder}/Mapping/${input_data}.bam" "$genome_annot" > "$output"
	done
fi

echo "##############################  COV WINDOWS  STEP  ######################################"
# STEP : Windows vocerage - Count number of reads per window, median expression of window and transcript (for POI detection)
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $WINDOW_STEP; then
	for index in ${!DATANAMES[@]}; do
		ip_data=${IP_DATA[index]}
		input_data=${INPUT_DATA[index]}
		if [[ " ${peak_detections[*]} " == *" POI "* ]]; then
		  	peak_technique=POI
			echo "Run windows and genes coverage of" $ip_data "for" $peak_technique "peak detection technique"
			sh src/bash/windows_cov.sh $folder $annotation_file $ip_data $peak_technique
			echo "Run windows and genes coverage of" $input_data "for" $peak_technique "peak detection technique"
			sh src/bash/windows_cov.sh $folder $annotation_file $input_data $peak_technique
		fi

		if [[ " ${peak_detections[*]} " == *" Fisher "* || " ${peak_detections[*]} " == *" RPMF "* ]]; then
			peak_technique=RPMF
			echo "Run windows coverage of" $ip_data "for" $peak_technique "peak detection technique"
			sh src/bash/windows_cov.sh $folder $annotation_file $ip_data $peak_technique
			echo "Run windows coverage of" $input_data "for" $peak_technique "peak detection technique"
			sh src/bash/windows_cov.sh $folder $annotation_file $input_data $peak_technique
		fi
	done
fi

echo "##############################  PEAK DETECTION STEP  ######################################"
# STEP : PeakDetection - Run the different peak detection techniques
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $PEAK_DETECTION_STEP; then
	for index in ${!DATANAMES[@]}; do
		dataname=${DATANAMES[index]}
		ip_data=${IP_DATA[index]}
		input_data=${INPUT_DATA[index]}
		for index_technique in ${!peak_detections[@]}; do
			peak_technique=${peak_detections[index_technique]}
			echo "Run" $peak_technique "peak detection for" $dataname
			sh src/bash/peak_detection.sh $folder $dataname $ip_data $input_data $peak_technique
		done
	done
fi

echo "##############################  MULTIQC STEP  ######################################"  ##
# STEP : MultiQC 
if $MULTIQC_STEP; then
	outdir=$folder/MultiQC/
	multiqc -i MERIPSeq -o $outdir .
fi

echo "##############################  MEDIAN COVERAGE STEP  ######################################"  ##
# STEP : MultiQC 
if $MULTIQC_STEP; then
	python src/python/create_mean_wig.py -p $folder -e $project_name -b all
fi


echo "##############################  FINALIZE PEAK DETECTION STEP ######################################"
# STEP : PeakDetection 2 - Finalize peak detection by regrouping all the datasets together
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $PEAK_FINALIZE_STEP; then
	for index_technique in ${!peak_detections[@]}; do
		peak_technique=${peak_detections[index_technique]}
		echo "Finalize" $peak_technique "peak detection for" $project_name
		python src/python/peak_finalize.py -p "$folder" -e $project_name -t $peak_technique -a $annotation_file
	done
fi

echo "##############################  PEAK ANNOTATION STEP ######################################"
# STEP : PeakAnnotation - Annotater les differents résultats de detection de pics
# WARNING - Optional step
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $PEAK_ANNOTATE_STEP; then
	for index_technique in ${!peak_detections[@]}; do
		peak_technique=${peak_detections[index_technique]}
		echo "Annotate peaks detected by" $peak_technique "for" $project_name
		python src/python/peak_annotation.py -p "$folder" -e $project_name -t $peak_technique -a $annotation_file -g $genome_file -b $bed_name_raw
	done
fi

echo "##############################  SEARCH MAX COVERAGE STEP ######################################"
# STEP : Search Max Coverage of every peak and extract peak_size peaks centered on it for 
if $SEARCH_MAX_STEP; then
	for index_technique in ${!peak_detections[@]}; do
		peak_technique=${peak_detections[index_technique]}
		echo "Search Max Coverage of every peak and extract "$peak_size" peaks centered on it for "$peak_technique " in " $project_name
		python src/python/find_peak_max.py -p $folder -e $project_name -t $peak_technique -s $peak_size -g $genome_file -b $bed_name_raw
	done
fi

echo "##############################  REGROUP PEAK DETECTION STEP  ######################################"
# PeakRegroup - Regroup the result from different peak detection techniques
if $REGROUP_PEAK_STEP; then
	echo "Regroup peaks from different peak detection technique:" ${peak_detections[@]}
	technique_bed=''
	for index_technique in ${!peak_detections[@]}; do
		peak_technique=${peak_detections[index_technique]}
		technique_bed="$technique_bed $folder/PeakDetection/$peak_technique/${project_name}_${peak_technique}_${bed_name_max}.bed"
	done
	raw_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name_regroup}_regroup_raw.bed
	sort_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name_regroup}_regroup_sort.bed
	merge_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name_regroup}_regroup.bed
	final_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name_regroup}.bed
	cat $technique_bed > "$raw_bed"
	sort -k1,1 -k2,2n $raw_bed > $sort_bed
	bedtools merge -i "$sort_bed" -c 4,5 -o distinct > $merge_bed
	awk '{OFS = \"\\t\"; if(split(\$5,t,\",\")>=${min_number_technique}){print \$1,\$2,\$3,\$4,\$5,\"+\"}}' $merge_bed > $final_bed
	rm $sort_bed
	rm $raw_bed
	rm $merge_bed

	echo "Search Max Coverage of every peak and extract "$peak_size" peaks centered on it, save to "${project_name}_${bed_name_regroup}
	python src/python/find_peak_max.py -p $folder -e $project_name -t All -s $peak_size -g $genome_file -b $bed_name_regroup
	
	echo "Annotate peaks detected by at least "$min_number_technique" for "${project_name}_${bed_name_regroup_max}
	python src/python/peak_annotation.py -p "$folder" -e $project_name -t '' -a $annotation_file -g $genome_file -b $bed_name_regroup_max
fi

echo "##############################  HTSEQ PEAK STEP  ######################################"
# HTSeq - Run the different peak detection techniques
# WARNING - This step is very long and should be parallelized for each datanameition and run on a cluster
if $HTSEQ_PEAK_STEP; then
	# suffix of the peak file
	final_bed_name=MaxMaxValues_${min_number_technique}
	for index in ${!DATANAMES[@]}; do
		ip_data=${IP_DATA[index]}
		input_data=${INPUT_DATA[index]}
		project_name=Liver
		genome_annot=${folder}/PeakDetection/Peaks/${project_name}_${bed_name_regroup_max}.gtf
		output=${folder}/PeakDetection/RawCov/HTSeq_${project_name}_${ip_data}_${bed_name_regroup_max}_peaks.txt
		echo "Run HTSeq count for peaks for" $ip_data
		htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' "${folder}/Mapping/${ip_data}.bam" "$genome_annot" > "$output"
		output=${folder}/PeakDetection/RawCov/HTSeq_${project_name}_${input_data}_${bed_name_regroup_max}_peaks.txt
		echo "Run HTSeq count for peaks for" $input_data
		htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' "${folder}/Mapping/${input_data}.bam" "$genome_annot" > "$output"
	done
fi

echo "##############################  DIFFERENTIAL PEAK  STEP  ######################################"
# STEP : Differential peak analysis
if $DIFF_PEAK_STEP; then
	source activate py36 # because some packages R do not work with python 2.7, remove if you have no problem)
	echo "Run peak differential expression analysis for" $project_name
	folder_diff=$folder/PeakDiffExpression/
	Rscript src/R/diff_methyl.R ${folder_diff} ${project_name} ${annotation_file} ${bed_name_regroup_max}
fi


echo "##############################  DIFF SPLICING STEP  ######################################"
# STEP : Differential splicing analysis - study difference of exon
if $DIFF_PEAK_STEP; then
	echo "Run Differential splicing analysis for" $project_name
	dir_salmon=${folder}/Salmon
	index=${folder}/Genome/transcripts_index
	for index in ${!DATANAMES[@]}; do
		echo "Run Salmon analysis for "$input_data
		salmon quant -i $index -l A -r ${folder}/FastQ/${input_data}.fastq  --validateMappings -o ${dirout}/$dataName
	done
	echo "Finalize diff splicing analysis for" $project_name
	gtf_salmon=${folder}/Genome/$annotation_file.annotation.gtf
	Rscript src/R/analysis_salmon_dexseq.R ${project_name}  ${dir_salmon} 7 ${gtf_salmon}
fi

echo "##############################  MOTIF PRESENCE ANALYSIS STEP  ######################################"
# STEP : Motif presence analysis
if $MOTIF_STEP; then
	echo "Search for Motif presence"
	motif_file_description=$folder/example/${project_name}_Motif.txt
	db_mouse=$folder/Genome/meme_db
	meme_minw=4
	meme_maxw=7
	meme_nmotifs=3
	meme_mod=anr
	peak_name=${project_name}_${bed_name_regroup_max}
	if [ ! -f $motif_file_description ]; then
		echo "ERROR - Cannot find Motif file list: " $motif_file_description
		exit 1
	else
		echo "Read Motif list file: " $motif_file_description
		index=0
		while read p; do
			type_data=$(echo "$p" | awk '{print $2}')		
			file_name=${peak_name}_${type_data}
			output=$folder/PeakDiffExpression/${peak_name}/Motif/results/${file_name}
			fasta_file=$folder/PeakDiffExpression/${peak_name}/Motif/fasta/${file_name}.fasta
			meme-chip -oc ${output} -time 300 -order 1 -meme-mod $meme_mod -meme-minw $meme_minw -meme-maxw $meme_maxw -meme-nmotifs $meme_nmotifs -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ${fasta_file}
			echo "Finish Motif search for" ${file_name}
		done
	fi

	echo "Regroup motifs"
	python src/python/motif/parse_motif_search.py -p $folder -e $project_name -b $bed_name_regroup_max

	echo "Create Fimo script"
	fimo_file=${folder}'PeakDiffExpression/Motif/Fimo_'${peak_name}'.sh'
	python src/python/motif/create_fimo_sh -p $folder -e $project_name -b $bed_name_regroup_max -f $fimo_file
	echo "Run Fimo script"
	sh $fimo_file

	echo "Calculate motif table"
	peak_name=${project_name}_${bed_name_regroup_max}
	python src/python/motif/create_motif_table.py -p $folder -e $peak_name -m 'Motif'_$peak_name

	echo "Create motif presence heatmap"
	Rscript src/R/motif_analysis.R ${project_name} ${bed_name_regroup_max}

	echo "Create motif figure with all motifs classified"
	python src/python/motif/regroup_figures_motif.py -p $folder -e $project_name -b $bed_name_regroup_max

fi

echo "##############################  GUITARPLOT ANALYSIS STEP  ######################################"
# STEP : GuitarPlot analysis (long processing time)
if $GUITAR_STEP; then
	echo "Process GuitarPlot"
	guitarplot_file=$folder/example/${project_name}_GuitarPlot.txt
	if [ ! -f $guitarplot_file ]; then
		echo "ERROR - Cannot find GuitarPlot file list: " $guitarplot_file
		exit 1
	else
		echo "Read GuitarPlot list file: " $guitarplot_file
		index=0
		while read p; do
			exp_design_name=$(echo "$line" | awk '{print $1}')
			bed_name=$(echo "$line" | awk '{print $2}')
			type_plot=$(echo "$line" | awk '{print $3}')

			file_name=${exp_design_name}_${bed_name}_${type_plot}
			Rscript src/R/run_guitar_plot.R $exp_design_name $bed_name $type_plot
		done
	fi
fi

echo "##############################  PEAK VENN STEP  ######################################"
# STEP : Compare listof peaks together
if $VENN_STEP; then
	Rscript src/R/venns.R
fi
