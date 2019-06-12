rm(list=ls())                 # remove all the objects from the R session
library(Guitar)
general_path <- "/pasteur/projets/policy01/m6aAkker/"
#general_path <- "/Users/cbecavin/Documents/m6aAkker/"
rscript_folder <- paste(general_path,"MeRIPSeq/src/m6aAnalysis/",sep="")
pkgload::load_all(path = rscript_folder, reset = TRUE)

args <- commandArgs(trailingOnly = TRUE)
exp_design_name <- args[1]
bed_name <- args[2]
type_plot <- args[3]

type_plot = gsub("\r","",type_plot)

#exp_design_name <- "LivOld"
#bed_name <- "Raw"
#type_plot <- "CONV_ZT13-Am_ZT13"

peak_name <- peak_name <- paste(exp_design_name,"_",bed_name,sep="")
annotation_name <- "gencode.vM13"

path_guitar <- paste(general_path,"PeakDiffExpression/",peak_name,"/GUITARPlot/",sep="")
setwd(path_guitar)
print("Calculate GuitarPlots in")
figure_folder <- paste("figures/",sep="")
if (!dir.exists(figure_folder)){
    dir.create(figure_folder)
}
print(getwd())
print(path_guitar)
print(peak_name)
print(type_plot)

# load targe_epi
target_epi_filename = paste(general_path,"PeakDiffExpression/",peak_name,"/",exp_design_name,"_target.txt",sep = "")
target_epi_fold <- read.table(file=target_epi_filename, header = TRUE)
print("load contrast")
load_contrast()

#all_biocond = paste0(all_biocond, ".T3")
all_biocond = gsub('Liver_', '', all_biocond)
# all_biocond = gsub('_T3', '.T3', all_biocond)
# all_biocond = gsub('_T13', '.T13', all_biocond)
# all_biocond <- all_biocond[grep("T13",all_biocond)]

# create design file
if(type_plot == "All_only"){
  list_names <- list("")
  list_files <- c(paste("../../PeakDetection/Peaks/",peak_name,".bed",sep = ""))
}else if(type_plot == "All"){
  list_names <- list("","_Meth","_NoMeth","_MethGene","_MethNoGene")
  list_files <- c(paste("../../PeakDetection/Peaks/",peak_name,".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[2]],".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[3]],".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[4]],".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[5]],".bed",sep = ""))
}else if(type_plot == "Diff_Meth"){
  list_names <- paste(all_biocond,sep="")
  list_files <- paste(path_guitar,"bed/",list_names,".bed",sep="")
}else if(type_plot == "Diff_MethNoGene"){
  list_names <- paste(all_biocond,"_MethNoGene",sep="")
  list_files <- paste(path_guitar,"bed/",list_names,".bed",sep="")
}else if(type_plot == "Diff_MethvsGene"){
  list_names <- paste(all_biocond,"_MethvsGene",sep="")
  list_files <- paste(path_guitar,"bed/",list_names,".bed",sep="")
}else{
  list_names <- c(type_plot,paste(type_plot,"_MethvsGene",sep=""),paste(type_plot,"_MethNoGene",sep=""))
  list_files <- paste(path_guitar,"bed/",list_names,".bed",sep="")
}
print(list_names)
print(list_files)

# Init Guitar plot
print("Init GUITAR plot")
guitar_file <-paste(annotation_name,"_guitar.RData",sep="")
if(!file.exists(guitar_file)){
  gencodeVM13 <- GenomicFeatures::makeTxDbFromGFF(paste(general_path,"Genome/", annotation_name,".annotation.gtf",sep=""),format = "gtf")
  gc_mm10_txdb <- Guitar::makeGuitarCoordsFromTxDb(gencodeVM13)
  save(gc_mm10_txdb,file = guitar_file)
}else{
  load(guitar_file)
}

feature_mm10 <- list()
for(k in 1:length(list_files)){
  filename <- list_files[k]
  if(file.exists(filename)){
    print(paste("Add:",list_names[k],"from",filename))
    bed_peaks <- rtracklayer::import.bed(filename)
    name <- paste(list_names[k],"(",length(bed_peaks),")",sep="")
    feature_mm10[[name]] = bed_peaks
  }
}

# Display guitar plot
print("Display GUITAR plot")
png(filename = paste(figure_folder,"GUITAR_",peak_name,"_",type_plot,".png",sep=""), width = 2000, height = 1000,
    res = 150)
Guitar::GuitarPlot(gfeatures = feature_mm10,GuitarCoordsFromTxDb = gc_mm10_txdb,rescaleComponent=TRUE)
dev.off()
