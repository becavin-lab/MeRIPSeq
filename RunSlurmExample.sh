#!/bin/bash
# run the script using
# Run.sh should be executable
# chmod 777 Run.sh
# nohup ./Run.sh > Run.log 2>&1 &

folder=/pasteur/projets/policy01/m6aAkker
folderScript=$folder/MeRIPSeq/src/bash_slurm
#project=CecAm
#project=Cecum
#project=Liver
project=LiverZTall
listfile=$folder/ExpDesign/${project}_List.txt
nbarray=$(wc -l $listfile | awk '{print $1}')
YOUR_EMAIL=christophe.becavin@pasteur.fr

echo $nbarray

# run Mapping
sbatch --array=${array} $folderScript/RunMapping.sh $folder/example/${project}_List.txt

# seq depth for creating Seqdepth/STAR_nbReads.txt
sbatch --array=${array} $folderScript/RunSeqDepth.sh $folder/example/RenameData.txt

# count number of reads per windows, and calculate median coverage for POI
sbatch --array=${array} $folderScript/RunWindowCovRPMF.sh $folder/example/${project}_List.txt
sbatch --array=${array} $folderScript/RunWindowCovPOI.sh $folder/example/${project}_List.txt

# create median wig for Searching max courage
sbatch $folderScript/RunMeanBigWig.sh $project all

# Run Peak detection for every dataset
sbatch --array=${array} $folderScript/RunPeakDetection.sh $folder/example/${project}_exp_design.txt MACS2
sbatch --array=${array} $folderScript/RunPeakDetection.sh $folder/example/${project}_exp_design.txt POI
sbatch --array=${array} $folderScript/RunPeakDetection.sh $folder/example/${project}_exp_design.txt RPMF
sbatch --array=${array} $folderScript/RunPeakDetection.sh $folder/example/${project}_exp_design.txt Fisher

# Finalize peaks detection : filter by peak occurence and regroup with bedtools
sbatch $folderScript/RunFinalize.sh $project POI
sbatch $folderScript/RunFinalize.sh $project MACS2
sbatch $folderScript/RunFinalize.sh $project RPMF
sbatch $folderScript/RunFinalize.sh $project Fisher
# for every technique ->  
# final_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + ‘_Raw.bed'

# Annotate peaks, search for overlapping genes, create gif and bed file, overlap with ref peaks
sbatch $folderScript/RunAnnotation.sh $project Fisher Raw
sbatch $folderScript/RunAnnotation.sh $project POI Raw
sbatch $folderScript/RunAnnotation.sh $project RPMF Raw
sbatch $folderScript/RunAnnotation.sh $project MACS2 Raw
# final_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + ‘_Raw.bed’

# Search for Max coverage = summit of the peaks
bed_name=Raw
sbatch $folderScript/RunSearchMax.sh $project MACS2 $bed_name
sbatch $folderScript/RunSearchMax.sh $project RPMF $bed_name
sbatch $folderScript/RunSearchMax.sh $project POI $bed_name
sbatch $folderScript/RunSearchMax.sh $project Fisher $bed_name
 # final_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + ‘_MaxValues.bed’

# Regroup peaks from different techniques
min_number_technique=3
bed_name=Raw
sbatch $folderScript/RunRegroup.sh $project $bed_name $min_number_technique
sbatch $folderScript/RunAnnotation.sh $project All ${bed_name}_$min_number_technique
bed_name=MaxValues
sbatch $folderScript/RunRegroup.sh $project $bed_name $min_number_technique
sbatch $folderScript/RunAnnotation.sh $project All ${bed_name}_$min_number_technique

# Search for max coverage center of he peak
bed_name_max=${bed_name}_${min_number_technique}
sbatch $folderScript/RunSearchMax.sh $project All $bed_name_max

#Run HTSeq for peaks
bed_name=MaxMaxValues
bed_name_peak=${bed_name}_$min_number_technique
sbatch --array=${array} $folderScript/RunHTSeqPeaks.sh $folder/example/${project}_List.txt $project $bed_name_peak

# Run differential analysis
bed_name=MaxMaxValues_3
sbatch $folderScript/RunDiffMeth.sh $project $bed_name

# Run Salmon Analysis
sbatch --array=${array} $folderScript/RunSalmon.sh $folder/example/${project}_List.txt
sbatch $folderScript/RunSalmonAnalysis.sh ${project}

# Run Motif search and Fimo
sbatch --array=1-108 RunMotif.sh $folder/example/${project}_Motif.txt
echo "Regroup motifs"
python $folder/MeRIPSeq/src/python/motif/parse_motif_search.py -p $folder -e $project_name -b $bed_name_regroup_max

echo "Create Fimo script"
fimo_file=${folder}'PeakDiffExpression/Motif/Fimo_'${peak_name}'.sh'
python $folder/MeRIPSeq/src/python/motif/create_fimo_sh -p $folder -e $project_name -b $bed_name_regroup_max -f $fimo_file

echo "Run Fimo script"
sbatch --array=1-3 $folderScript/RunFimo.sh $fimo_file

echo "Calculate motif table"
peak_name=${project_name}_${bed_name_regroup_max}
python $folder/MeRIPSeq/src/python/motif/create_motif_table.py -p $folder -e $peak_name -m 'Motif'_$peak_name

echo "Create motif presence heatmap"
Rscript $folder/MeRIPSeq/src/R/motif_analysis.R ${project_name} ${bed_name_regroup_max}

echo "Create motif figure with all motifs classified"
python $folder/MeRIPSeq/src/python/motif/regroup_figures_motif.py -p $folder -e $project_name -b $bed_name_regroup_max

# run Guitarplot and motif search
sbatch --array=1-117 $folderScript/RunGuitarPlot.sh $project/example/${project}_GuitarPlot.txt

# Run Venn diagram for Liver vs Cecum comparison
Rscript $folder/MeRIPSeq/src/R/venns.R


echo "Your MeRIPSeq pipeline ended "$(date "+%A %B %d at %H:%M:%S") | mailx -s "MeRIPSeq pipeline" $YOUR_EMAIL