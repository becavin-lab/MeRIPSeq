#!/bin/bash
#SBATCH --job-name move
#SBATCH --array=1-1
#SBATCH --output temp/move.out
#SBATCH --error temp/move.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=5GB
#SBATCH -p hubbioit

scratch=/pasteur/scratch/cbecavin/m6aAkker/
folder=/pasteur/projets/policy01/m6aAkker/
final=/pasteur/scratch/cbecavin/sample/

#line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
#OLD=$(echo "$line" | awk '{print $1}')
#dataName=$(echo "$line" | awk '{print $1}')

#NEW=$(echo "$line" | awk '{print $2}')

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# Create script SH to run in qsub
scriptSH="""#!/bin/bash -l
module add meme/4.10.1
sh $folder/motif.sh
#cp $scratch/RNASeq_raw/${OLD}.fastq $final/RNASeq_raw/${NEW}.fastq
#cp $scratch/FastQ/${OLD}.fastq $final/FastQ/${NEW}.fastq
#cp $folder/Mapping/${OLD}.bam $final/Mapping/${NEW}.bam
#cp $folder/Mapping/${OLD}.bam.bai $final/Mapping/${NEW}.bam.bai
#cp $folder/Mapping/STAR_results/${OLD}/${OLD}Log.final.out $final/Mapping/${NEW}Log.final.out
#cp $folder/Mapping/STAR_results/${OLD}/${OLD}Log.out $final/Mapping/${NEW}Log.out
#cp $folder/Mapping/STAR_results/${OLD}/${OLD}Log.progress.out $final/Mapping/${NEW}Log.progress.out
#cp $folder/MappingRNA/${OLD}.log $final/MappingRNA/${NEW}.log

#cp $folder/Expression/HTSeq_${OLD}.txt $final/Expression/HTSeq_${NEW}.txt
#cp $folder/Expression/HTSeq_${OLD}_peaks.txt $final/Expression/HTSeq_${NEW}_peaks.txt

"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
