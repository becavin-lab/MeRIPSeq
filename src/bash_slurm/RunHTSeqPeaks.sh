#!/bin/bash
#SBATCH --job-name HTPeak
#SBATCH --output temp/HTPeak.out
#SBATCH --error temp/HTPeak.err
#SBATCH --mail-type FAIL
#SBATCH --qos=hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
project_name=$2
gtf_name=$3

# prepare sbatch varialbes
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${gtf_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${gtf_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
#genome_annot=${folder}/PeakDetection/Peaks/${project_name}_MaxValue.gtf
genome_annot=${folder}/PeakDetection/Peaks/${project_name}_${gtf_name}.gtf
output=${folder}/Expression/HTSeq_${project_name}_${dataName}_${gtf_name}_peaks.txt
		
# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add HTSeq/0.9.1
# count reads per peak
htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${dataName}.bam $genome_annot > $output
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
