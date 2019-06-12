#!/bin/bash
#SBATCH --job-name salmon
#SBATCH --output /pasteur/projets/policy01/m6aAkker/temp/salmon.out
#SBATCH --error /pasteur/projets/policy01/m6aAkker/temp/salmon.err
#SBATCH --mail-type FAIL
#SBATCH --qos fast
#SBATCH --cpus-per-task 6
#SBATCH --mem=10GB
#SBATCH -p common
#SBATCH --array=1-141%47

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/users/abiton/abiton/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)


runTrimming=False

echo $dataName
dataName2=`echo $dataName | sed 's/B1770.2-//'  | sed 's/B1770-//' | sed 's/_R1_001.*//' | sed 's/_R1_001.*//' `

echo $dataName
echo $dataName2

# prepare sbatch varialbes
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
        mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# general variables
genome=GRCm38.primary_assembly.genome.fa
genomeFolder=$folder/Genome/
readFile=$scratch/FastQ/${dataName}.fastq


transcripts_fasta=$genomeFolder/gencode.vM13.transcripts.fa
dirout=$folder/salmon 
index=$genomeFolder/transcripts_index
#salmon=/pasteur/homes/abiton/tools/Salmon-latest_linux_x86_64/bin/salmon
fasta=$genome_folder/gencode.vM13.annotation.transcript.gfa


mkdir $dirout
# build index
#salmon index -t $fasta -i transcripts_index -k 31

echo """""" > ${shFile}

if $runTrimming; then
    # if needed, run Alien trimmer as it was done for the mapping, for the files not yet available in scratch
    scriptSH1="""
    #!/bin/sh -l                                                                                                                                                             
    module add AlienTrimmer/0.4.0                                                                                                                                                  
    gunzip $scratch/RNASeq_raw/${dataName}.fastq.gz                                                                                                                                
    AlienTrimmer -i $scratch/RNASeq_raw/${dataName}.fastq -c $folder/Genome/alienTrimmerPF8contaminants.fasta -o $scratch/FastQ/${dataName}.fastq                                  
    #rm $scratch/RNASeq_raw/${dataName}.fastq
    """
    echo  """$scriptSH1""" >> ${shFile}
fi

scriptSH="""
echo 'Processing sample ${dataName}'
if [ ! -f $dirout/${dataName}/quant.sf ]; then
#salmon quant -i $salmonFolder -l A \
#         -r ${fn} \
#         -p 6 -o $dirout/${samp}_quant
salmon quant -i $index -l A -r $scratch/FastQ/${dataName}.fastq  --validateMappings -p ${SLURM_CPUS_PER_TASK} -o $dirout/$dataName
fi
"""

# save script and run with qsub                                                                                                                                                                                                                       
echo """$scriptSH""" >> ${shFile}

#echo $scriptSH                                                                                                                                                                                                                                   
# Run script                                                                                                                                                                                                                                            
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
