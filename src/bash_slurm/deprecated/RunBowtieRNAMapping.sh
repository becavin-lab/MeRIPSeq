#!/bin/bash
#SBATCH --job-name BMapping
#SBATCH --array=1-8
#SBATCH --output temp/Mapping_fastq.out
#SBATCH --error temp/Mapping_fastq.err
#SBATCH --mail-type FAIL
#SBATCH --qos=hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=20GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/cbecavin/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
#dataName=Makegenome
#dataName=Conv_B_Input
#echo $dataName
#echo $folder

# prepare sbatch varialbes
mkdir $folder/temp/${SLURM_JOB_NAME}
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
genomeRRNA=gencode.vM13.rRNA
genomeFolder=$folder/Genome/
mappingRNAFolder=$folder/MappingRNA/${dataName}
readFile=$scratch/FastQ/${dataName}.fastq

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add bowtie2
# create genome
#bowtie2-build $genomeFolder/$genome.fasta $genomeFolder/$genome
# map datasets
bowtie2 -p 3 --local -x $genomeFolder/$genome -U $readFile -S $mappingFolder.sam 2>$mappingFolder.log
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
