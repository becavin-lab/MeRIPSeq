#!/bin/bash
#SBATCH --job-name SeqDep
#SBATCH --output temp/SeqDep.out
#SBATCH --error temp/SeqDep.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem=5GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/abiton/m6aAkker/
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
dataName=$(echo "$line" | awk '{print $2}')
rawData=$(echo "$line" | awk '{print $1}')

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable


# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
sh $scriptFolder/src/bash/seqdepth.sh $folder $scratch $dataName $rawData
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
