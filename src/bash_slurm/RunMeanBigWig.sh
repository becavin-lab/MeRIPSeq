#!/bin/bash
#SBATCH --job-name MeanBW
#SBATCH --output temp/MeanBW.out
#SBATCH --error temp/MeanBW.err
#SBATCH --array=1-1
#SBATCH --mail-type FAIL
#SBATCH --qos fast #hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
project_name=$1
biocond=$2

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${project_name}_${biocond}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${project_name}_${biocond}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

chromSize=$folder/Genome/chrNameLength.txt

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
source activate wiggle
python $scriptFolder/src/python/create_mean_wig.py -p $folder -e $project_name -b $biocond
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}