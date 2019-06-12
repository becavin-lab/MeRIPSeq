#!/bin/bash
#SBATCH --job-name Fimo
#SBATCH --output temp/Fimo.out
#SBATCH --error temp/Fimo.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=20GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker/

script=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# Create script SH to run in qsub
scriptSH="""#!/bin/bash -l
module add meme/4.10.1
$script
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
