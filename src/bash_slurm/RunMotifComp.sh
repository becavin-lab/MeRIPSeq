#!/bin/bash
#SBATCH --job-name MotifComp
#SBATCH --output temp/MotifComp.out
#SBATCH --error temp/MotifComp.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
peak_name=$(echo "$line" | awk '{print $1}')
type_data=$(echo "$line" | awk '{print $2}')

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variables
#fasta_file=$folder/PeakDiffExpression/$peak_name/Motif/fasta/${file_name}.fasta
fasta_file=$folder/PeakDiffExpression/$peak_name/Motif/fasta/${type_data}.fasta
if [ ! -d $folder/PeakDiffExpression/$peak_name/Motif/results ]; then
	mkdir $folder/PeakDiffExpression/$peak_name/Motif/results
fi
output=$folder/PeakDiffExpression/$peak_name/Motif/results/${type_data}
db_mouse=$folder/Genome/meme_db
meme_minw=4
meme_maxw=7
meme_nmotifs=3
meme_mod=anr

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add meme/4.10.1
meme-chip -oc $output -time 300 -order 1 \
-meme-mod $meme_mod -meme-minw $meme_minw -meme-maxw $meme_maxw -meme-nmotifs $meme_nmotifs \
-dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 $fasta_file
echo Finish Motif search
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
