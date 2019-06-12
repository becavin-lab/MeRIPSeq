#!/bin/bash
#sh DetectCrashedJobs.sh > TarsProblem.log

logline="overall alignment rate"
uniq="Uniquely mapped reads %"
multiple="% of reads mapped to multiple loci"
toomany="% of reads mapped to too many loci"
short="% of reads unmapped: too short"

#gc="%GC	"

rRNAfolder=/pasteur/projets/policy01/m6aAkker/MappingRNA
starfolder=/pasteur/projets/policy01/m6aAkker/Mapping/STAR_results
fastqcfolder=/pasteur/projets/policy01/m6aAkker/FastQC


# for dataName in $(ls $rRNAfolder/*.log)
# do	
# 	search=$(grep "${logline}" $dataName)
# 	length=${#search}
# 	repl="% overall alignment rate"
# 	rep=""
# 	echo $dataName $search
	
# done

for dataName in $(ls $starfolder/*/*Log.final.out)
do	
	search1=$(grep "${uniq}" $dataName)
	search2=$(grep "${multiple}" $dataName)
	search3=$(grep "${toomany}" $dataName)
	search4=$(grep "${short}" $dataName)
	# length=${#search}
	# repl="% overall alignment rate"
	# rep=""
	echo $dataName $search1 $search2 $search3 $search4
	
done

# for dataName in $(ls $fastqcfolder/*_fastqc/fastqc_data.txt)
# do	
# 	search1=$(grep "${gc}" $dataName)
	
# 	# length=${#search}
# 	# repl="% overall alignment rate"
# 	# rep=""
# 	echo $dataName $search1
	
# done

