#!/bin/bash
#SBATCH --job-name Final
#SBATCH --output temp/Final.out
#SBATCH --error temp/Final.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 1
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
project_name=$1
peak_technique=$2

# prepare sbatch variables
mkdir $folder/temp/${SLURM_JOB_NAME}
shFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
annotation_file=gencode.vM13
genome_file=GRCm38.primary_assembly.genome

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add bedtools/2.25.0
module load Python/2.7.11
#source activate py27
echo Finalize $peak_technique peak detection for $project_name
python $scriptFolder/src/python/peak_finalize.py -p $folder -e $project_name -t $peak_technique -a $annotation_file -b Raw
echo End finalize for $peak_technique
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
