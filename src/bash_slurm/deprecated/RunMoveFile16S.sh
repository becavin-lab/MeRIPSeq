#!/bin/bash
#SBATCH --job-name move
#SBATCH --array=1-8
#SBATCH --output temp/move.out
#SBATCH --error temp/move.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=5GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker/

line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
OLD=$(echo "$line" | awk '{print $1}')
dataName=$(echo "$line" | awk '{print $1}')

NEW=$(echo "$line" | awk '{print $2}')

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log
# Create script SH to run in qsub
scriptSH="""#!/bin/bash -l

cp /pasteur/projets/policy01/BioIT/amine/jabs/nov17/${OLD} $folder/16S/${NEW}

"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
