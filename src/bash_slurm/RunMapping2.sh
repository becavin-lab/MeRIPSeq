#!/bin/bash
#SBATCH --job-name Map
#SBATCH --output temp/Map.out
#SBATCH --error temp/Map.err
#SBATCH --mail-type FAIL
#SBATCH --qos normal #fast
#SBATCH --cpus-per-task 3
#SBATCH --mem=30GB
#SBATCH -p common #dedicated #common 
#SBATCH --array=1-2%2

folder=/pasteur/projets/policy01/m6aAkker
scratch=/pasteur/scratch/abiton/m6aAkker
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

dataName2=`echo $dataName | sed 's/B1770.2-//'  | sed 's/B1770-//' | sed 's/_R1_001.*//' | sed 's/_R1_001.*//' `

echo $dataName
echo $dataName2
#dataName=Makegenome
#dataName=Conv_B_Input
#echo $dataName
#echo $folder

# prepare sbatch varialbes
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# general variables
genome=GRCm38.primary_assembly.genome.fa
genomeFolder=$folder/Genome/
fastqcFolder=$folder/FastQC/
annotation=$folder/Genome/gencode.vM13.annotation.gtf
readFile=$scratch/FastQ/${dataName}.fastq
#readFile=$folder/RNASeq_raw/all_gz/${dataName}.fastq.gz

# variables for STAR
mappingFolder=$folder/Mapping/STAR_results/${dataName}
chromSize=$folder/Genome/chrNameLength.txt
wigfile=$folder/Mapping/STAR_results/${dataName}/${dataName}Signal.Unique.str1.out
bamFileUnfiltered=$folder/Mapping/STAR_results/${dataName}/${dataName}Aligned.sortedByCoord.out.bam
bamFileFiltered=$folder/Mapping/${dataName}
genome_annot=${folder}/Genome/gencode.vM13.annotation.exon.gtf
wig_data=$folder/PeakDetection/RawCov/${dataName}.wig

# variable for Bowtie = mappin on rRNA genome
genomeRRNA=gencode.vM13.rRNA
mappingRNAFolder=$folder/MappingRNA/${dataName}

# variables for HTSeq
outputHTSeq=${folder}/Expression/HTSeq_${dataName}.txt

# variables for RSeqQC
genomeRSeQC=/pasteur/projets/policy01/m6aAkker/Genome/mm10_GENCODE_VM11_basic.bed
bamfile=$folder/Mapping/${dataName}.bam
outputDist=$folder/RSeqQC/${dataName}_distrib.txt
outputTIN=$folder/RSeqQC/${dataName}_tin
outputStat=$folder/RSeqQC/${dataName}_stat.txt
outputDupl=$folder/RSeqQC/${dataName}_dupl
outputQual=$folder/RSeqQC/${dataName}_qual
outputGene=$folder/RSeqQC/${dataName}_gene


# Create script SH to run in qsub
scriptSH="""
#!/bin/sh -l
module add AlienTrimmer/0.4.0
#gunzip $scratch/RNASeq_raw/${dataName}.fastq.gz
#AlienTrimmer -i $scratch/RNASeq_raw/${dataName}.fastq -c $folder/Genome/alienTrimmerPF8contaminants.fasta -o $scratch/FastQ/${dataName}.fastq
#rm $scratch/RNASeq_raw/${dataName}.fastq 

module load STAR/2.5.0a
# # Create genome
# #STAR --runThreadN ${SLURM_CPUS_PER_TASK} --runMode genomeGenerate --genomeDir $genomeFolder --genomeFastaFiles $genomeFolder/$genome

# #Map dataset
 mkdir ${mappingFolder}

# # # Run basic STAR with normalized WIG
# STAR --runThreadN ${SLURM_CPUS_PER_TASK} --genomeDir ${genomeFolder} --outFileNamePrefix ${mappingFolder}/${dataName} \
#   --sjdbGTFfile ${annotation} \
#    --sjdbOverhang 100 --readFilesIn ${readFile}\
#    --outSAMtype BAM SortedByCoordinate --outWigType wiggle \
#    --outWigStrand Unstranded --outWigNorm RPM
# STAR --runThreadN ${SLURM_CPUS_PER_TASK} --genomeDir ${genomeFolder}  --outFileNamePrefix ${mappingFolder}/${dataName2} \
#    --sjdbGTFfile ${annotation} \
#    --readFilesCommand zcat\
#    --sjdbOverhang 100 --readFilesIn ${readFile} \
#    --outSAMtype BAM SortedByCoordinate --outWigType wiggle \
#    --outWigStrand Unstranded --outWigNorm RPM


# # # create bigwig
#module load wigToBigWig
#wigToBigWig ${wigfile}.wig ${chromSize} ${wigfile}.bw
##/pasteur/homes/cbecavin/opt/wigtools/wigToBigWig ${wigfile}.wig ${chromSize} ${wigfile}.bw
#mv ${wigfile}.bw $folder/Mapping/${dataName}.bw

# # filter bam files by removing data with score = 0
 if [ ! -f ${bamFileFiltered}.bam.bai ]; then
module add samtools
samtools view -b -q 1 ${bamFileUnfiltered} > ${bamFileFiltered}.bam
# #samtools sort ${bamFileFiltered}.bam ${bamFileFiltered}
samtools index ${bamFileFiltered}.bam
fi

# Calculate raw coverage wig file for POI peak detection technique : 
module add bedtools
bedtools genomecov -d -split -ibam $folder/Mapping/${dataName}.bam | awk '\$3!=0 {print \$0}' > ${wig_data}

# map reads on rRNA 
module add bowtie2
# create genome
#bowtie2-build $genomeFolder/$genome.fasta $genomeFolder/$genome
# map datasets
#bowtie2 -p 3 --local -x $genomeFolder/$genomeRRNA -U $readFile -S $mappingRNAFolder.sam 2>$mappingRNAFolder.log

# Run FastQC
#module add fastqc
#fastqc -o $fastqcFolder $readFile

# Count reads per genes
module add HTSeq/0.9.1
# # count reads per genes
 htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${dataName}.bam $genome_annot > $outputHTSeq

# Run quality control on bam file
#module add R
#module add RSeQC
#read_distribution.py  -i $bamfile -r $genomeRSeQC > $outputDist
#bam_stat.py  -i $bamfile > $outputStat
#read_duplication.py -i $bamfile -o $outputDupl
#read_quality.py -i $bamfile -o $outputQual
#tin.py -i $bamfile -r $genomeRSeQC
#geneBody_coverage.py -r $genomeRSeQC -i $bamfile  -o $outputGene

echo \"Everything is finished\"
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
