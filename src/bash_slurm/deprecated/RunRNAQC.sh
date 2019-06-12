#!/bin/bash
#SBATCH --job-name RNAQC
#SBATCH --array=1-49
#SBATCH --output temp/RNAQC.out
#SBATCH --error temp/RNAQC.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem=30GB
#SBATCH -p hubbioit

#folder=/pasteur/projets/policy01/m6aAkker
folder=/pasteur/scratch/cbecavin/sample
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
dataName=$(echo "$line" | awk '{print $2}')

# prepare variables
mkdir /pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
bamfile=$folder/Mapping/${dataName}
genome=/pasteur/projets/policy01/m6aAkker/Genome
# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add picard-tools/1.94
module add samtools
AddOrReplaceReadGroups I=${bamfile}.bam O=${bamfile}_GATK.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
samtools index ${bamfile}_GATK.bam

if [ ! -d -o ${folder}/RNAQC/${dataName} ]; then
	mkdir ${folder}/RNAQC/${dataName}
fi

module add bwa
module add RNA-SeQC
RNA-SeQC -t ${genome}/gencode.vM13.transcript.gtf -singleEnd -r ${genome}/GRCm38.primary_assembly.genome.fa -o ${folder}/RNAQC/${dataName} -s \"${dataName}|${bamfile}_GATK.bam|NA\"

"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
