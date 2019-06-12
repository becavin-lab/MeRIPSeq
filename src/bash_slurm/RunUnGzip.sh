#!/bin/bash
#SBATCH --job-name ungzip
#SBATCH --array=1-48
#SBATCH --output temp/ungzip.out
#SBATCH --error temp/ungzip.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=5GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker/
scratch=/pasteur/scratch/users/cbecavin/
echo ${SLURM_ARRAY_TASK_ID}
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
rawName=$(echo "$line" | awk '{print $1}')
dataName=$(echo "$line" | awk '{print $2}')

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

echo $rawName
echo $dataName
echo $folder
# variable for Bowtie = mappin on rRNA genome
genomeFolder=$folder/Genome/
fastqcFolder=$folder/FastQC/
readFile=$scratch/FastQ/${dataName}.fastq
genomeRRNA=gencode.vM13.rRNA
mappingRNAFolder=$folder/MappingRNA/${dataName}

# Create script SH to run in qsub
scriptSH="""#!/bin/bash -l
module add AlienTrimmer/0.4.0                                                                                                                                                  
cp $folder/RNASeq_raw/all_gz/$rawName.fastq.gz $scratch/RNASeq_raw/all_gz/$rawName.fastq.gz
unpigz -p 3 -k $scratch/RNASeq_raw/all_gz/$rawName.fastq.gz
mv $scratch/RNASeq_raw/all_gz/$rawName.fastq $scratch/RNASeq_raw/$dataName.fastq
AlienTrimmer -i $scratch/RNASeq_raw/${dataName}.fastq -c $folder/Genome/alienTrimmerPF8contaminants.fasta -o $scratch/FastQ/${dataName}.fastq                                  
rm $scratch/RNASeq_raw/${dataName}.fastq 

# map reads on rRNA 
module add bowtie2
# # map datasets
bowtie2 -p 3 --local -x $genomeFolder/$genomeRRNA -U $readFile -S $mappingRNAFolder.sam 2>$mappingRNAFolder.log

## Run FastQC
module add fastqc
fastqc -o $fastqcFolder $readFile


"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
