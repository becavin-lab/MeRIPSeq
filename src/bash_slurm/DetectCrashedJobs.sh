#!/bin/bash
#sh DetectCrashedJobs.sh > TarsProblem.log


#type=Map
#logline="Everything is finished"

#type=WinCovP
#logline="Median transcript expression calculated"

#type=WinCovR
#logline="SAM alignments  processed."

#type=HTSeq
#logline="SAM alignments  processed."

#type=PDetect  # for MACS2
#logline="Done!"

#type=PDetect  # for Fisher
#logline="P-value correction calculated"

#type=PDetect  # for RPMF
#logline="RPMF Score calculated"

#type=PDetect  # for POI
#logline="POI Score calculated"

type=Motif
logline="Finish Motif search"

#type=MotifComp
#logline="Finish Motif search"

#type=Fimo  # for RPMF
#logline="CANCELLED AT "

#type=HTPeak
#type=HTSeqRun
#logline="SAM alignments  processed."
#logline="GFF lines processed."

tempfolder=/pasteur/projets/policy01/m6aAkker/temp/$type

for dataName in $(ls $tempfolder/*.log)
do	
	search=$(grep "${logline}" $dataName)
	length=${#search}
	if [ $length != 0 ] 
	then
		echo $search $dataName
		shfile=${dataName/.log/.sh}
		#echo $shfile
		rm ${dataName}
		rm ${shfile}
	fi
done


