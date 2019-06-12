#!/bin/bash
# Title: MeRIP-Seq pipeline analysis
# Author: Christophe Bécavin, Hub Bioinformatic and Biostatistic
# 28 rue du Docteur Roux, Institut Pasteur, Paris
#
# Executable:
# sh setup.sh folder
#
# PARAMETERS:
# folder   folder of your data
#
###################################################################################

# Setup variables
folder=$1

###################################################################################
echo "Setup folders for the analysis"

temp_folder=$folder'/Exp_Design/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/RNASeq_raw/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Genome/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/temp/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Seqdepth/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/FastQ/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/RSeqQC/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/FastQC/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/MultiQC/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Mapping/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Salmon/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/MappingRNA/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Mapping/STAR_results/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Expression/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/RawCov/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/temp/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/Peaks/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/POI/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/RPMF/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/Fisher/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDetection/MACS2/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/All_Figures/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/All_Figures/HeatMap/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/All_Figures/PCA/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/All_Figures/Motif/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/PeakDiffExpression/All_Figures/GuitarPlot/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi

temp_folder=$folder'/Motif/'
if [ ! -d $temp_folder ]; then
	mkdir $temp_folder
fi
