#!/bin/bash
#SBATCH --job-name Guitar
#SBATCH --output temp/Guitar.out
#SBATCH --error temp/Guitar.err
#SBATCH --qos=hubbioit
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 3
#SBATCH --mem=20GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
line=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
exp_design_name=$(echo "$line" | awk '{print $1}')
bed_name=$(echo "$line" | awk '{print $2}')
type_plot=$(echo "$line" | awk '{print $3}')

file_name=${exp_design_name}_${bed_name}_${type_plot}

# prepare variables
mkdir /pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}
shFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=/pasteur/projets/policy01/m6aAkker/temp/${SLURM_JOB_NAME}/${file_name}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# variable
script_r=$folder/MeRIPSeq/src/run_guitar_plot.R
# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
source activate rscript
Rscript $script_r $exp_design_name $bed_name $type_plot
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
