#!/bin/bash
#SBATCH --job-name HTSeq
#SBATCH --array=1-88
#SBATCH --output temp/HTSeq.out
#SBATCH --error temp/HTSeq.err
#SBATCH --mail-type FAIL
#SBATCH --qos=hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare sbatch varialbes
mkdir $folder/temp/${SLURM_JOB_NAME}
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
input_data=${dataName}
#genome_annot=${folder}/Genome/gencode.vM13.annotation.transcript.gtf
#genome_annot=${folder}/Genome/gencode.vM13.annotation.exon.gtf
genome_annot=${folder}/Genome/gencode.vM13.annotation.gtf
#output=${folder}/Expression/HTSeq_Transcript_${input_data}.txt
#output=${folder}/Expression/HTSeq_Exon_${input_data}.txt
output=${folder}/Expression/HTSeq_${input_data}.txt

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add HTSeq/0.9.1
# count reads per genes
#htseq-count -t exon -i exon_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${input_data}.bam $genome_annot > $output
#htseq-count -t exon -i transcript_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${input_data}.bam $genome_annot > $output
htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${input_data}.bam $genome_annot > $output
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
