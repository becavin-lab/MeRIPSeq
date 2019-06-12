#!/bin/bash
#SBATCH --job-name RSeqQC
#SBATCH --array=1-120
#SBATCH --output temp/RSeqQC.out
#SBATCH --error temp/RSeqQC.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem=40GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
#folder=/pasteur/scratch/cbecavin/sample
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare variables
mkdir /pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
genomeAnnot=/pasteur/projets/policy01/m6aAkker/Genome/mm10_GENCODE_VM11_basic.bed
bamfile=$folder/Mapping/${dataName}.bam
outputDist=$folder/RSeqQC/${dataName}_distrib.txt
outputTIN=$folder/RSeqQC/${dataName}_tin
outputStat=$folder/RSeqQC/${dataName}_stat.txt
outputDupl=$folder/RSeqQC/${dataName}_dupl
outputQual=$folder/RSeqQC/${dataName}_qual
outputGene=$folder/RSeqQC/${dataName}_gene


# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add R
module add RSeQC

read_distribution.py  -i $bamfile -r $genomeAnnot > $outputDist
#bam_stat.py  -i $bamfile > $outputStat
#read_duplication.py -i $bamfile -o $outputDupl
#read_quality.py -i $bamfile -o $outputQual
#tin.py -i $bamfile -r $genomeAnnot
#geneBody_coverage.py -r $genomeAnnot -i $bamfile  -o $outputGene
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${outputStat} bash ${shFile}
