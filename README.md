# MeRIP-Seq pipeline analysis

Author: Christophe Bécavin, Anne Biton, Hub Bioinformatic and Biostatistic

28 rue du Docteur Roux, Institut Pasteur, Paris

GitLab: https://gitlab.pasteur.fr/cbecavin/MeRIPSeq/

Executable (2 versions are available one for local computing, another for slurm cluster computing)
### For local analysis

`sh run.sh project_name`

For example: sh run.sh Liver

### For Slurm cluster

`nohup ./RunSlurmExample.sh > runServer.log 2>&1 &`

RunSlurmExample.sh should be executable `chmod 777 RunSlurmExample.sh`

## PARAMETERS:

project-name   -  Name of the project

The list of file should be in your_folder/Exp_Design/project_name_exp_design.txt
Your sequencing data should be in : your_folder/RNASeq_raw/
our experiment design file should be a tab-separated table with at least 5 columns:

|DataName|	IP_dataset|	Input_dataset|	BioCond|	Seq|
| ------ | ------ | ------ | ------ | ------ |
|Am_ZT3_1|	Am_ZT3_IP_1|	Am_ZT3_Input_1|	Am|	1|
|Am_ZT13_1|	Am_ZT13_IP_1|	Am_ZT13_Input_1|	Am|	2|
|CONV_ZT3_1|	CONV_ZT3_IP_1|	CONV_ZT3_Input_1|	CONV|	1|
|CONV_ZT13_4|	CONV_ZT13_IP_4|	CONV_ZT13_Input_4|	CONV|	2|

(See example/Liver_exp_design.txt)

## PIPELINE Description:

Setup true or false each step of the workflow to run it or not

INSTALL first all dependencies: See INSTALL.txt

1 - setup.sh - Prepare Analysis by creating all necessary folders

2 - genome.sh - Prepare genome sequence, annotation, and windows for peak detection

### Run on every dataset:

WARNING - Steps 3 to 7 are run separately for each datasets.

WARNING - They might take a ot of time to perform, so they can be parallelized on a cluster.

3 - Trimming, mapping, BAM file filtering and quality control, And WIG calculation for peak detection

4 - Calculate SeqDepth for peak detection techniques

5 - HTSeq - Count number of reads per genes

6 - coverage_window.py - Count number of reads per window

7 - peak_detection.py - Run Fisher, POI,and RPMF peak detection techniques

8 - MACS2 - Run MACS2 peak detection

### Run only one time per project_name

9 - Search for median coverage for max coverage search

10 - MultiQC for quality control

11 - finalize peak detection by filtering out "bad" peaks, searchning for max coverage and peak center

12 - (Optionnal) Annotate all peaks from all different techniques

13 - Search max coverage of every peaks and center peak position on it.

14 - Regroup peaks from all the detection techniques, annotate them, find overlap position with genes and referent MeRIP-Seq, CLIPSeq and TREW

### Run on every dataset:

15 - HTSeq - Count number of reads in each methylation sites detected

### Run only one time per project_name

16 - Perform differential methylation sites analysis in R

17 - Run differential splicing events

18 - Run Motif presence algorithm

19 - Run GuitarPlot for every list of peak

20 - (Optional) Compare two peak lists to Ref peaks