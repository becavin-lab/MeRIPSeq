#!/bin/bash
#SBATCH --job-name Regroup
#SBATCH --output temp/Regroup.out
#SBATCH --error temp/Regroup.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=5GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
project_name=$1
bed_name=$2
min_nb_peaktechnique=$3

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${project_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${project_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
annotation_file=gencode.vM13
genome_file=GRCm38.primary_assembly.genome

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add bedtools/2.25.0
source activate py27
echo Regroup peaks from different peak detection technique: ${peak_detections[@]}
technique_bed=''
peak_detections=(MACS2 RPMF POI Fisher)
for index_technique in \${!peak_detections[@]}; do
 	peak_technique=\${peak_detections[index_technique]}
 	technique_bed=\"\$technique_bed $folder/PeakDetection/\$peak_technique/${project_name}_\${peak_technique}_${bed_name}.bed\"
done
raw_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name}_regroup_raw.bed
sort_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name}_regroup_sort.bed
merge_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name}_regroup_merge.bed
final_bed=$folder/PeakDetection/Peaks/${project_name}_${bed_name}_${min_nb_peaktechnique}.bed
cat \$technique_bed > \$raw_bed
sort -k1,1 -k2,2n \$raw_bed > \$sort_bed
bedtools merge -i \$sort_bed -c 4,5 -o distinct > \$merge_bed
awk '{OFS = \"\\t\"; if(split(\$5,t,\",\")>=${min_nb_peaktechnique}){print \$1,\$2,\$3,\$4,\$5,\"+\"}}' \$merge_bed > \$final_bed
rm \$sort_bed
rm \$raw_bed
rm \$merge_bed
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
