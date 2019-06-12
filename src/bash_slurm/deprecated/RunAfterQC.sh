#!/bin/bash
#SBATCH --job-name AfterQC
#SBATCH --array=1-49
#SBATCH --output temp/AfterQC.out
#SBATCH --error temp/AfterQC.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/cbecavin/sample

line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
dataName=$(echo "$line" | awk '{print $2}')

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
fastqFile=$scratch/RNASeq_raw/${dataName}.fastq
afterFold=$scratch/AfterQC/
fastqcFolder=$folder/FastQC/


# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
source activate py27
python /pasteur/homes/cbecavin/opt/AfterQC/after.py -1 $fastqFile -g $afterFold/good -b $afterFold/bad -r $afterFold/QC
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
