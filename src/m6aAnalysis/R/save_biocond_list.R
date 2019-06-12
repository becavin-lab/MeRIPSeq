save_biocond_list <- function(){
  print("Run save_biocond_list")
  # # # Play with gene lists
  if ('Time' %in% colnames(target_epi)) {
    target_epi$BioCond <- (paste0(as.character(target_epi$BioCond), '.', as.character(target_epi$Time)))
    target_epi$labels <- sub('-T', '.T',as.character(target_epi$labels))
  }
  # } else if (exp_design_name == 'Liver_all_T13') {
  #   target_epi$labels <- sub('-T', '.T',as.character(target_epi$labels))
  #   rownames(target_epi) <-    target_epi$labels
  # }

  bioconds <- unique(as.character(target_epi$BioCond))
  data_to_bioconds <- list()
  norm_median_epi <- norm_voom_epi
  for(i in 1:length(bioconds)){
    data_biocond <- as.character(target_epi$labels[target_epi$BioCond==bioconds[i] & target_epi$TypeData=="IP"])
    data_to_bioconds[[bioconds[i]]] = data_biocond
    norm_median_epi[,bioconds[i]] = apply(norm_median_epi[,data_biocond], 1, function(x) median(x))
  }

  norm_median_epi <- norm_median_epi[,bioconds]
  norm_median_epi$ID <- rownames(norm_median_epi)
  norm_median_epi <- merge(norm_median_epi, peaks, by.all = "ID")
  rownames(norm_median_epi) <-norm_median_epi$ID
  #hist(norm_median_epi[,1],breaks=100)
  #hist(norm_median_epi$CONV_ZT13,breaks=100)
  #hist(norm_median_epi$GF_ZT13,breaks=100)

  # norm_sd_epi <- norm_voom_epi
  # for(i in 1:length(bioconds)){
  #   data_biocond <- as.character(target_epi$labels[target_epi$BioCond==bioconds[i] & target_epi$TypeData=="IP"])
  #   data_to_bioconds[[bioconds[i]]] = data_biocond
  #   norm_sd_epi[,bioconds[i]] = apply(norm_sd_epi[,data_biocond], 1, function(x) IQR(x))
  # }
  # norm_sd_epi <- norm_sd_epi[,bioconds]
  # norm_sd_epi$ID <- rownames(norm_sd_epi)
  # norm_sd_epi <- merge(norm_sd_epi, peaks, by.all = "ID")
  # rownames(norm_sd_epi) <-norm_sd_epi$ID

  # save everything in bed file for guitarplots
  value_cutoff = 2
  path_guitar <- paste(general_path,"PeakDiffExpression/",peak_name,"/GUITARPlot/bed/",sep="")
  list_type <- list("","_MethvsGene","_MethNoGene")
  for(k in 1:length(list_type)){
    type <- list_type[[k]]
    for(i in 1:length(bioconds)){
      biocond = bioconds[i]
      list_peaks <- diff_peaks
      if(type == "_MethvsGene"){
        list_peaks <- diff_peaks_genes
      }else if(type == "_MethNoGene"){
        list_peaks <- diff_peaks[!diff_peaks_NOgenes]
      }
      norm_median_epi_filter <- norm_median_epi[list_peaks,]
      norm_median_epi_filter <- norm_median_epi_filter[norm_median_epi_filter[,biocond] > value_cutoff,]
      rownames(norm_median_epi_filter) = norm_median_epi_filter$ID
      bed_col <- c("chromo_peak","begin_peak","end_peak","ID","Motif","strand")
      if(nrow(norm_median_epi_filter) > 0){
        # save in bedfile for GUITAR plot
        bed_file <- unique(norm_median_epi_filter[,bed_col])
        norm.filename = paste(path_guitar,biocond,type,".bed",sep = "")
        write.table(bed_file, file=norm.filename, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
      }
    }
  }

  # Extract fasta sequences
  path_motif <- paste(general_path,"PeakDiffExpression/",peak_name,"/Motif/fasta/",sep="")
  list_type <- list("","_MethvsGene","_MethNoGene")
  for(k in 1:length(list_type)){
    type <- list_type[[k]]
    for(i in 1:length(bioconds)){
      biocond = bioconds[i]
      list_peaks <- diff_peaks
      if(type == "_MethvsGene"){
        list_peaks <- diff_peaks_genes
      }else if(type == "_MethNoGene"){
        list_peaks <- diff_peaks[!diff_peaks_NOgenes]
      }
      norm_median_epi_filter <- norm_median_epi[list_peaks,]
      matrix <- norm_median_epi_filter[norm_median_epi_filter[,biocond] > value_cutoff,]
      seq_fasta <- fasta[names(fasta) %in% matrix$ID]
      # save all peaks in fasta
      #print(paste("Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
      norm.filename = paste(folderFasta,biocond,type,".fasta",sep = "")
      seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)

    }
  }
}
