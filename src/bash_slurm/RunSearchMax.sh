#!/bin/bash
#SBATCH --job-name Max
#SBATCH --output temp/Max.out
#SBATCH --error temp/Max.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=10GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
scriptFolder=/pasteur/projets/policy01/m6aAkker/MeRIPSeq/
project_name=$1
peak_technique=$2
bed_name=$3

# prepare sbatch variables
mkdir $folder/temp/${SLURM_JOB_NAME}
shFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${peak_technique}_${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
annotation_file=gencode.vM13
genome_file=GRCm38.primary_assembly.genome
size=150

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
#source activate py27
echo SearchMax $peak_technique peak detection for $project_name
python $scriptFolder/src/python/find_peak_max.py -p $folder -e $project_name -t $peak_technique -s $size -g $genome_file -b $bed_name
echo End SearchMax for $peak_technique
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
