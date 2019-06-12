diff_summary <- function(){
  # Diff Summary
  print("Save summary contrasts")
  folderBed <- "GUITARPlot/bed/"

  diff_peaks_NOgenes <- diff_peaks[!(diff_peaks %in% diff_peaks_genes)]
  diff_peaks_NO_Meth <- peaks$ID[!(peaks$ID %in% diff_peaks)]

  list_matrix <- list(diff_peaks,diff_Epi,diff_peaks_genes,diff_peaks_NOgenes,diff_peaks_pos_neg,diff_peaks_pos_pos,diff_peaks_neg_neg,diff_peaks_neg_pos)
  list_name <- list("_Meth","_Meth_Epi","_MethGene","_MethNoGene","_Meth+Gene-","_Meth+Gene+","_Meth-Gene-","_Meth-Gene+")
  folderContrast <- paste("DiffSummary/",sep="")
  if (!dir.exists(folderContrast)){
    dir.create(folderContrast)
  }

  folderFasta <- paste("Motif/fasta/",sep="")
  if (!dir.exists(folderFasta)){
    dir.create(folderFasta)
  }
  assign("folderFasta",folderFasta, pos = globalenv())

  summary_table <- c("Type","Nb peaks","Nb genes")
  for(k in 1:length(list_name)){
      matrix <- fold_epi[list_matrix[[k]],]
      rownames(matrix) = matrix$ID
      type <- list_name[[k]]
      norm.filename = paste(folderContrast,peak_name,type,".xls",sep = "")
      write.table(matrix, file=norm.filename, quote=FALSE, sep = "\t", row.names = FALSE)
      norm.filename = paste(folderContrast,peak_name, type,"_gene.txt",sep = "")
      write.table(unique(matrix$gene_name), file=norm.filename, quote=FALSE, row.names = FALSE, col.names = FALSE)
      norm.filename = paste(folderContrast,peak_name,type,"_uniprot.txt",sep = "")
      write.table(unique(matrix$UniprotID), file=norm.filename, quote=FALSE, row.names = FALSE, col.names = FALSE)

      # save bed file and fasta
      bed_col <- c("chromo_peak","begin_peak","end_peak","ID","Motif","strand")
      bed_file <- unique(matrix[bed_col])
      if(length(bed_file$chromo_peak)!=0){
        norm.filename = paste(folderBed,peak_name,type,".bed",sep = "")
        write.table(bed_file, file=norm.filename, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

        }

      # save fasta
      # save in fasta for motif search
        seq_fasta <- fasta[names(fasta) %in% matrix$ID]
        # save all peaks in fasta
        #print(paste("Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
        norm.filename = paste(folderFasta,peak_name,type,".fasta",sep = "")
        seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)

        if(type=="_Meth"){
          type_overlaps = c("3UTR","5UTR","CDS","Intron")
          for(type_overlap in type_overlaps){
            matrix_type = matrix[grep(type_overlap,matrix$Type_overlap),]
            seq_fasta <- fasta[names(fasta) %in% matrix_type$ID]
            #print(paste(type_overlap,"Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
            norm.filename = paste(folderFasta,peak_name,type,"_",type_overlap,".fasta",sep = "")
            seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)
          }
        }

    # print summary
    summary <- c(type,length(matrix$ID),length(unique(matrix$gene_name)))
    summary_table <- rbind(summary_table,summary)
    #print(summary)
  }
  write.table(summary_table, file=paste(peak_name,"_Diff_Summary.txt",sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  # save peaks not diff meth
  matrix <- fold_epi[diff_peaks_NO_Meth,]
  type <- "_NoMeth"
  seq_fasta <- fasta[names(fasta) %in% matrix$ID]
  bed_col <- c("chromo_peak","begin_peak","end_peak","ID","Motif","strand")
  bed_file <- unique(matrix[bed_col])
  if(length(bed_file$chromo_peak)!=0){
    norm.filename = paste(folderBed,peak_name,type,".bed",sep = "")
    write.table(bed_file, file=norm.filename, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

  }

  # save all peaks in fasta
  print("Save all peaks in fasta")
  #print(paste("Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
  norm.filename = paste(folderFasta,peak_name,type,".fasta",sep = "")
  seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)

  # plot histogram of motif distrib
  list_matrix <- list(diff_peaks,diff_peaks_genes,diff_peaks_NOgenes)
  list_name <- list("_Meth","_MethGene","_MethNoGene")
  motif_distrib = c()
  motif_length = c()
  motif_type = c()
  for(k in 1:length(list_name)){
    matrix <- fold_epi[list_matrix[[k]],]
    type <- list_name[[k]]
    # Put motif presence in a vector
    if(length(matrix$Motif)!=0){
      for(w in 1:length(matrix$Motif)){
        motif_distrib <- rbind(motif_distrib, matrix$Motif[w])
        motif_length <- rbind(motif_length, matrix$length_peak[w])
        motif_type <- rbind(motif_type, paste(peak_name,type,sep = ""))
      }
    }
  }
  motif_df = data.frame(motif_distrib,motif_type,motif_length)
  png(filename = paste(figure_folder,"MOTIF_",peak_name,"_distrib.png",sep=""), width = 2000, height = 1000,res = 150)
  ggplot2::ggplot(motif_df, ggplot2::aes(x = motif_distrib, fill = motif_type)) + ggplot2::geom_density(alpha = 0.3) + ggplot2::facet_wrap(~ motif_type)
  dev.off()
  png(filename = paste(figure_folder,"MOTIF_",peak_name,"_dist_length.png",sep=""), width = 2000, height = 1000,res = 150)
  ggplot2::ggplot(motif_df, ggplot2::aes(x = motif_distrib/motif_length, fill = motif_type)) + ggplot2::geom_density(alpha = 0.3) + ggplot2::facet_wrap(~ motif_type)
  dev.off()
  png(filename = paste(figure_folder,"MOTIF_",peak_name,"_length.png",sep=""), width = 2000, height = 1000, res = 150)
  ggplot2::ggplot(motif_df, ggplot2::aes(x = motif_length, fill = motif_type)) + ggplot2::geom_density(alpha = 0.3) + ggplot2::facet_wrap(~ motif_type)
  dev.off()
  ggplot2::ggplot(motif_df, ggplot2::aes(x = motif_distrib, fill = motif_type)) + ggplot2::geom_histogram(binwidth = 1,alpha = 0.3)
  #ggplot(patients, aes(x = age)) + geom_density() + facet_wrap(~ study)

  #dev.off()
}


##### OBSOLETE
# # save fasta by cutting every peaks in overlapping region 2*motif_len (default = 50)
# seq_fasta <- fasta[names(fasta) %in% matrix$ID]
# motif_len <- 25 # half of motif_len desired
# new_fasta <- list()
# for(peak_id in names(seq_fasta)){
#   seq <- seq_fasta[[peak_id]]
#   center <- floor(length(seq)/2)
#   start_pos <- center - motif_len
#   end_pos <- center + motif_len
#   first_interval <- c(max(start_pos,1), min(end_pos,length(seq)))
#   # go on the left
#   interval <- first_interval
#   index_seq = 0
#   while((interval[2]-interval[1]) == 2*motif_len){
#     #print(interval)
#     index_seq = index_seq + 1
#     new_id <- paste(peak_id,"_l_",index_seq,sep="")
#     new_seq <- as.SeqFastadna(seq[start_pos:end_pos], name = new_id, Annot = new_id)
#     new_fasta[[new_id]] <- new_seq
#     start_pos <- max(start_pos - motif_len,1)
#     end_pos <- min(end_pos - motif_len,length(seq))
#     interval <- c(start_pos, end_pos)
#   }
#   # go on the right
#   start_pos <- center - motif_len
#   end_pos <- center + motif_len
#   first_interval <- c(max(start_pos,1), min(end_pos,length(seq)))
#   interval <- first_interval
#   index_seq = 0
#   while((interval[2]-interval[1]) == 2*motif_len){
#     #print(interval)
#     index_seq = index_seq + 1
#     new_id <- paste(peak_id,"_r_",index_seq,sep="")
#     new_seq <- as.SeqFastadna(seq[start_pos:end_pos], name = new_id, Annot = peak_id)
#     new_fasta[[new_id]] <- new_seq
#     start_pos <- max(start_pos + motif_len,1)
#     end_pos <- min(end_pos + motif_len,length(seq))
#     interval <- c(start_pos, end_pos)
#   }
# }
# print(paste("Extract:",length(new_fasta),"sub-sequences of length",motif_len,"from",nrow(matrix),"rows matrix"))
# norm.filename = paste(folderFasta,"Const_",contr_Epi_names[i],type,".fasta",sep = "")
# write.fasta(sequences=new_fasta,names=names(new_fasta),file.out=norm.filename)

