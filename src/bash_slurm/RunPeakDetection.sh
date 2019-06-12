#!/bin/bash
#SBATCH --job-name PDetect
#SBATCH --output temp/PDetect.out
#SBATCH --error temp/PDetect.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
peak_technique=$2
dataName=$(echo "$line" | awk '{print $1}')
ip_data=$(echo "$line" | awk '{print $2}')
input_data=$(echo "$line" | awk '{print $3}')

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
annotation_file=gencode.vM13
genome_file=GRCm38.primary_assembly.genome

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
#source activate py27
sh $scriptFolder/src/bash/peak_detection.sh $folder $dataName $ip_data $input_data $peak_technique $scriptFolder
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
