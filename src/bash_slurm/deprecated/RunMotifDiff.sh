#!/bin/bash
#SBATCH --job-name MotifDiff
#SBATCH --output temp/MotifDiff.out
#SBATCH --error temp/MotifDiff.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
peak_name=$1


# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${peak_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${peak_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variables
if [ ! -d $folder/PeakDiffExpression/$peak_name/Motif/results ]; then
	mkdir $folder/PeakDiffExpression/$peak_name/Motif/results
fi
PATH_MOTIF=/pasteur/projets/policy01/m6aAkker/PeakDiffExpression/$peak_name/Motif
PATH_RESULT=${PATH_MOTIF}/results/${peak_name}'_NEG'
PATH_NEG_fasta=$PATH_MOTIF/fasta/${peak_name}'_NoMeth.fasta'
PATH_fasta=$PATH_MOTIF/fasta/${peak_name}'_Meth.fasta'

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add meme/4.10.1
meme-chip -oc $PATH_RESULT -time 300 -neg $PATH_NEG_fasta -index-name meme-chip.html -order 1 -meme-mod anr -meme-minw 4 -meme-maxw 7 -meme-nmotifs 3 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 $PATH_fasta
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
