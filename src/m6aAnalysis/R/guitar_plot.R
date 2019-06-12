#' Title
#'
#' @return
#' @export
#'
#' @examples
guitar_plot <- function(){

  print("Plot GUITAR plot for DiffSummary lists")

  # Init Guitar plot
  print("Init GUITAR plot")
  guitar_file <-paste(general_path,"/temp/",annotation_name,"_guitar.RData",sep="")
  if(!file.exists(guitar_file)){
    gencodeVM13 <- GenomicFeatures::makeTxDbFromGFF(paste("../../Genome/", annotation_name,".annotation.gtf",sep=""),format = "gtf")
    gc_mm10_txdb <- Guitar::makeGuitarCoordsFromTxDb(gencodeVM13)
    save(gc_mm10_txdb,file = guitar_file)
  }else{
    load(guitar_file)
  }

  print("Calculate GuitarPlots in")
  figure_folder <- paste(path_guitar,"figures/",sep="")
  if (!dir.exists(figure_folder)){
    dir.create(figure_folder)
  }


  # Set data for Guitar plot
  #list_name <- list("")
  list_names <- list("","_Meth","_MethGene","_MethNoGene")

  #final_bed <- read.delim(paste("../../PeakDetection/Peaks/",peak_name,"_Final.bed",sep = ""),header = FALSE)
  #final_bed <- final_bed[sample(rownames(final_bed),5000),]
  #write.table(final_bed,paste("../../PeakDetection/Peaks/",peak_name,"_Final_Sample.bed",sep = ""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  list_files <- c(paste(path_guitar,"bed/",peak_name,".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[2]],".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[3]],".bed",sep = ""),
                  paste(path_guitar,"bed/",peak_name,list_names[[4]],".bed",sep = ""))
  feature_mm10 <- list()


  for(k in 1:length(list_files)){
    filename <- list_files[[k]]
    if(file.exists(filename)){
      print(filename)
      bed_peaks <- rtracklayer::import.bed(filename)
      feature_mm10[[list_names[[k]]]] = bed_peaks
    }
  }
  # Display guitar plot
  print("Display GUITAR plot")
  png(filename = paste(figure_folder,"GUITAR_",peak_name,".png",sep=""), width = 2000, height = 1000,
      res = 150)
  Guitar::GuitarPlot(gfeatures = feature_mm10,GuitarCoordsFromTxDb = gc_mm10_txdb,rescaleComponent=TRUE)
  dev.off()

  # diff_folder <- paste("DiffMethyl",sep="")
  # for(i in 1:length(contr_Epi_names)){
  #   feature_mm10 <- list()
  #   final_names <- list()
  #   # #feature_mm10 <- list(bed_peaks,bed_meth,bed_meth_gene,bed_meth_nogene)
  #   #list_type <- list("","_Epi", "_MethvsGene","_Meth+Gene-","_Meth+Gene+"
  #   #                  ,"_Meth-Gene+","_Meth-Gene-","_MethNoGene")
  #   list_type <- list("","_MethvsGene","_MethNoGene")

  #   folderContrast <- paste(diff_folder,"/",contr_Epi_names[i],"/",sep="")
  #   for(k in 1:length(list_type)){
  #     type <- list_type[[k]]
  #     norm.filename = paste(folderContrast,contr_Epi_names[i],type,".bed",sep = "")
  #     if(file.exists(norm.filename)){
  #       bed_peaks <- rtracklayer::import.bed(norm.filename)
  #       feature_mm10[[paste(contr_Epi_names[i],type,sep = "")]] = bed_peaks
  #     }
  #   }
  #   png(filename = paste(figure_folder,"GUITAR_",contr_Epi_names[i],".png",sep=""), width = 2000, height = 1000,res = 150)
  #   Guitar::GuitarPlot(gfeatures = feature_mm10, GuitarCoordsFromTxDb = gc_mm10_txdb,rescaleComponent=TRUE)
  #   dev.off()
  # }
}
