#!/bin/bash
#SBATCH --job-name Rscript
#SBATCH --array=1-168
#SBATCH --output temp/Rscript.out
#SBATCH --error temp/Rscript.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 3
#SBATCH --mem=5GB
#SBATCH -p hubbioit

#folder=/pasteur/projets/policy01/m6aAkker
folder=/pasteur/scratch/cbecavin/sample
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare variables
mkdir /pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
script_r=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/scripts/$dataName.R
# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
source activate rscript
Rscript $script_r
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
