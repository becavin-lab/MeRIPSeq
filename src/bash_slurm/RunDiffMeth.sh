#!/bin/bash
#SBATCH --job-name DiffMeth
#SBATCH --output temp/DiffMeth.out
#SBATCH --error temp/DiffMeth.err
#SBATCH --qos=fast #hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 3
#SBATCH --mem=20GB
#SBATCH -p dedicated #hubbioit

folder=/pasteur/projets/policy01/m6aAkker
exp_design_name=$1
bed_name=$2

script_r=$folder/MeRIPSeq/src/diff_methyl.R
file_name=${exp_design_name}_${bed_name}

# prepare variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module load R/3.4.1
#source activate rscript
Rscript $script_r $exp_design_name $bed_name
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
