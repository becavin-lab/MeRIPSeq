plot_description <- function(){
  # description plots
  varInt <- "BioCond"      # factor of interest
  majSequences <- descriptionPlots(counts=counts_epi, group=target_epi[,varInt],
                                   col=color_biocond[levels(target_epi[,varInt])])
  image_folder = paste(project_dir,"FiguresEpi/",sep="")
  if (!dir.exists(image_folder)){
    file.rename(from="figures/",to=image_folder)
  }else{
    unlink(image_folder, recursive=TRUE)
    file.rename(from="figures/",to=image_folder)
  }
  write.table(counts_epi[order(counts_epi[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, peak_name,"_Count_Epi.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")

  counts_epiIP <- counts_epi[, grep('_IP',colnames(counts_epi))]
  majSequences <- descriptionPlots(counts=counts_epiIP, group=target_epi[match(colnames(counts_epiIP),target_epi$labels),varInt],
                                   col=color_biocond[unique(as.character(target_epi[match(colnames(counts_epiIP),target_epi$labels),varInt]))])
  image_folder = paste(project_dir,"FiguresEpiIP/",sep="")
  if (!dir.exists(image_folder)){
    file.rename(from="figures/",to=image_folder)
  }else{
    unlink(image_folder, recursive=TRUE)
    file.rename(from="figures/",to=image_folder)
  }
  write.table(counts_epiIP[order(counts_epi[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, peak_name,"_Count_EpiIP.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")

  counts_epiInput <- counts_epi[, grep('_.nput',colnames(counts_epi))]
  majSequences <- descriptionPlots(counts=counts_epiInput, group=target_epi[match(colnames(counts_epiInput),target_epi$labels),varInt],
                                   col=color_biocond[unique(as.character(target_epi[match(colnames(counts_epiInput),target_epi$labels),varInt]))])
  image_folder = paste(project_dir,"FiguresEpiInput/",sep="")
  if (!dir.exists(image_folder)){
    file.rename(from="figures/",to=image_folder)
  }else{
    unlink(image_folder, recursive=TRUE)
    file.rename(from="figures/",to=image_folder)
  }
  write.table(counts_epiInput[order(counts_epi[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, peak_name,"_Count_EpiInput.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")


  majSequences <- descriptionPlots(counts=counts_RNA_I, group=target_RNA_I[,varInt],
                                   col=color_biocond[unique(as.character(target_RNA_I[match(colnames(counts_RNA_I),target_RNA_I$labels),varInt]))])
  image_folder = paste(project_dir,"FiguresRNAInput/",sep="")
  if (!dir.exists(image_folder)){
    file.rename(from="figures/",to=image_folder)
  }else{
    unlink(image_folder, recursive=TRUE)
    file.rename(from="figures/",to=image_folder)
  }
  write.table( counts_RNA_I[order(counts_RNA_I[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, peak_name,"_Count_RNA_Input.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")

  majSequences <- descriptionPlots(counts=counts_RNA_IP, group=target_RNA_IP[,varInt],
                                   col=color_biocond[unique(as.character(target_RNA_IP[match(colnames(counts_RNA_IP),target_RNA_IP$labels),varInt]))])
  image_folder = paste(project_dir,"FiguresEpiIP/",sep="")
  if (!dir.exists(image_folder)){
    file.rename(from="figures/",to=image_folder)
  }else{
    unlink(image_folder, recursive=TRUE)
    file.rename(from="figures/",to=image_folder)
  }
  write.table( counts_RNA_IP[order(counts_RNA_IP[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, peak_name,"_Count_RNA_IP.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  if (!dir.exists(figure_folder)){
    dir.create(figure_folder)
  }

  # majSequences <- descriptionPlots(counts=counts_Trans_I, group=target_RNA_I[,varInt], col=color_biocond)
  # image_folder = paste(project_dir,"FiguresTranscrI/",sep="")
  # if (!dir.exists(image_folder)){
  #   file.rename(from="figures/",to=image_folder)
  # }else{
  #   unlink(image_folder, recursive=TRUE)
  #   file.rename(from="figures/",to=image_folder)
  # }
  # write.table(counts_Trans_I[order(counts_Trans_I[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, project_name,"_Count_Trans_I.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  # if (!dir.exists(figure_folder)){
  #   dir.create(figure_folder)
  # }
  #
  # majSequences <- descriptionPlots(counts=counts_Exon_I, group=target_RNA_I[,varInt], col=color_biocond)
  # image_folder = paste(project_dir,"FiguresExonI/",sep="")
  # if (!dir.exists(image_folder)){
  #   file.rename(from="figures/",to=image_folder)
  # }else{
  #   unlink(image_folder, recursive=TRUE)
  #   file.rename(from="figures/",to=image_folder)
  # }
  # write.table(counts_Exon_I[order(counts_Exon_I[,1],decreasing = TRUE)[1:50],], file=paste(image_folder, project_name,"_Count_Exon_I.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  # if (!dir.exists(figure_folder)){
  #   dir.create(figure_folder)
  # }

}
