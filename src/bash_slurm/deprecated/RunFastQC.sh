#!/bin/bash
#SBATCH --job-name FastQC
#SBATCH --array=1-186
#SBATCH --output temp/FastQC.out
#SBATCH --error temp/FastQC.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=5GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/cbecavin/m6aAkker

dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
#dataName=Conv_B_Input
#echo $dataName
#echo $folder

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
fastqFile=$scratch/FastQ/${dataName}.fastq
fastqcFolder=$folder/FastQC/


# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add fastqc
fastqc -o $fastqcFolder $fastqFile
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
