#!/bin/bash
#SBATCH --job-name WinCovP
#SBATCH --output temp/WinCovP.out
#SBATCH --error temp/WinCovP.err
#SBATCH --mail-type FAIL
#SBATCH --qos fast
#SBATCH --cpus-per-task 1
#SBATCH --mem=30GB
#SBATCH -p common

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare sbatch variables
mkdir $folder/temp/${SLURM_JOB_NAME}
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
annotation_file=gencode.vM13
genome_file=GRCm38.primary_assembly.genome

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
# count median reads per transcript
module add bedtools
source activate py27
#calculate cov for POI
sh $scriptFolder/src/bash/windows_cov.sh $folder $annotation_file ${dataName} POI $scriptFolder
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
