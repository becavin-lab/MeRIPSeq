#!/bin/bash
#SBATCH --job-name MultiQC
#SBATCH --output temp/MultiQC.out
#SBATCH --error temp/MultiQC.err
#SBATCH --mail-type FAIL
#SBATCH --qos hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=30GB
#SBATCH -p hubbioit

folder=/pasteur/projets/policy01/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)
#dataName=Conv_B_Input
#echo $dataName
#echo $folder

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

outdir=$folder/MultiQC/

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add MultiQC
#Cecum all
multiqc -x *ZT13* -x *ZT3* -x *CONV_I* -x *vanco_I* -x *abx_I* -x *Ec* -x *ex_GF* -f -i Yo -n Yo -o $outdir .
#multiqc -x *ZT13* -f -i Cecum -n Cecum -o $outdir .

#Liver all
#multiqc -x *GF_I* -x *CONV_I* -x *Am_I* -x *_ZT3_* -x *vanco_I* -x *abx_I* -f -i Liver -n Liver -o $outdir .
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
