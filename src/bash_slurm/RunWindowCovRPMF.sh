#!/bin/bash
#SBATCH --job-name WinCovR
#SBATCH --array=1-1
#SBATCH --output temp/WinCovR.out
#SBATCH --error temp/WinCovR.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=50GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
#dataName=Conv_B_Input
#echo $dataName
#echo $folder

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
# calculate cov for RPMF and Fisher
module add HTSeq/0.9.1
sh $scriptFolder/src/bash/windows_cov.sh $folder $annotation_file ${dataName} RPMF $scriptFolder
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
