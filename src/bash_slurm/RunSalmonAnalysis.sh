#!/bin/bash                                                                                                                                                                                                               
#SBATCH --job-name salmonDEXSeq                                                                                                                                                                                       
#SBATCH --output temp/salmonDEXSeq_liver_tx5inMinCond.out                                                                                                                                                                         
#SBATCH --error temp/salmonDEXSeq_liver_tx5inMinCond.err                                                                                                                                                                                 
#SBATCH --array=1-1                                                                                                                                                                                                    
#SBATCH --mail-type FAIL                                                                                                                                                                                                
#SBATCH --qos hubbioit                                                                                                                                                                                           
#SBATCH --mem=30GB                                                                                                                                                                                                     
#SBATCH -p hubbioit                                                                                                                          
#SBATCH --cpus-per-task 7

exp_design_name=$1

folder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq
gtf=$folder/Genome/gencode.vM13.annotation.gtf


module load R/3.5.2
#Rscript  /pasteur/projets/policy01/m6aAkker/MeRIPSeq/scratch/transcript_usage/analysis_salmon.R Cecum_all  /pasteur/projets/policy01/m6aAkker/salmon/ 7 $gtf
Rscript  $folder/MeRIPSeq/src/R/analysis_salmon_dexseq.R $exp_design_name $folder//Salmon/ 7 $gtf
#Rscript  /pasteur/projets/policy01/m6aAkker/MeRIPSeq/scratch/transcript_usage/analysis_salmon_dexseq.R Liver_all  /pasteur/projets/policy01/m6aAkker/salmon/ 7 $gtf
#Rscript --no-restore  /pasteur/projets/policy01/m6aAkker/MeRIPSeq/scratch/transcript_usage/analysis_salmon.R Liver_all_T13 /pasteur/projets/policy01/m6aAkker/salmon/ 7 $gtf
