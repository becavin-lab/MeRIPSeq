
save_contrast <- function(){
  print("Save contrasts")
  diff_folder <- paste("DiffMethyl",sep="")
  if (!dir.exists(diff_folder)){
     dir.create(diff_folder)
  }
  folderBed <- "GUITARPlot/"
  if (!dir.exists(folderBed)){
    dir.create(folderBed)
  }
  folderBed <- "GUITARPlot/bed/"
  if (!dir.exists(folderBed)){
    dir.create(folderBed)
  }
  diff_peaks <- character()
  diff_genes <- character()
  diff_peaks_genes <- character()
  diff_peaks_pos_neg <- character()
  diff_peaks_pos_pos <- character()
  diff_peaks_neg_neg <- character()
  diff_peaks_neg_pos <- character()
  diff_Epi <- character()

  # Load fasta for peaks
  folderFasta <- paste("Motif/fasta/",sep="")
  if (!dir.exists(folderFasta)){
    dir.create(folderFasta)
  }

  motif_distrib = c()
  motif_length = c()
  motif_type = c()
  motif_name = c()
  summary_table <- c("Comparison","Type","Nb peaks","Nb genes")
  for(i in 1:length(contr_Epi_names)){
    folderContrast <- paste(diff_folder,"/",contr_Epi_names[i],"/",sep="")
    if (!dir.exists(folderContrast)){
      dir.create(folderContrast)
    }
    best <- limma::topTable(fit_bayes_epi, coef = contr_Epi_names[i], adjust="BH", p.value = 0.05, lfc = 1, number=Inf)
    best$ID = row.names(best)
    best_fold <- limma::topTable(fit_bayes_epi_fold, coef = contr_Epi_names[i], adjust="BH", p.value = 0.05, lfc = 1, number=Inf)
    best_fold$ID = row.names(best_fold)
    best_RNA_I <- limma::topTable(fit_bayes_RNA_I, coef = contr_Epi_names[i], adjust="BH", p.value = 0.05, lfc = 1, number=Inf)
    best_RNA_I$Gene_ID = row.names(best_RNA_I)

    diffmethyl = contr_Epi_names[i]
    diffgene = paste(contr_Epi_names[i],"_Input",sep="")
    fold_epi_pvalue <- fold_epi[best$ID,]
    fold_epi_pvalue_gene <- merge(fold_epi_pvalue,best_RNA_I,by.x = "Gene_ID",by.y = "Gene_ID",all.x=FALSE,all.y=FALSE)
    fold_epi_pvalue_NO_gene <- fold_epi_pvalue[!(fold_epi_pvalue$ID %in% fold_epi_pvalue_gene$ID),]

    fold_epi_pos_neg <- fold_epi_pvalue_gene[fold_epi_pvalue_gene[,diffmethyl] > 1 & fold_epi_pvalue_gene[,diffgene] < -1,]
    fold_epi_pos_pos <- fold_epi_pvalue_gene[fold_epi_pvalue_gene[,diffmethyl] > 1 & fold_epi_pvalue_gene[,diffgene] > 1,]
    fold_epi_neg_pos <- fold_epi_pvalue_gene[fold_epi_pvalue_gene[,diffmethyl] < -1 & fold_epi_pvalue_gene[,diffgene] > 1,]
    fold_epi_neg_neg <- fold_epi_pvalue_gene[fold_epi_pvalue_gene[,diffmethyl] < -1 & fold_epi_pvalue_gene[,diffgene] < -1,]

    diff_peaks <- c(diff_peaks, best$ID)
    diff_peaks_genes <- c(diff_peaks_genes, fold_epi_pvalue_gene$ID)
    diff_peaks_pos_neg <- c(diff_peaks_pos_neg, fold_epi_pos_neg$ID)
    diff_peaks_pos_pos <- c(diff_peaks_pos_pos, fold_epi_pos_pos$ID)
    diff_peaks_neg_neg <- c(diff_peaks_neg_neg, fold_epi_neg_neg$ID)
    diff_peaks_neg_pos <- c(diff_peaks_neg_pos, fold_epi_neg_pos$ID)
    diff_Epi <- c(diff_Epi, best_fold$ID)
    diff_genes <- c(diff_genes, best_RNA_I$Gene_ID)

    list_matrix <- list(fold_epi[best$ID,],fold_epi[best_fold$ID,],fold_RNA[best_RNA_I$Gene_ID,]
                        ,fold_epi_pvalue_gene,fold_epi_pos_neg,fold_epi_pos_pos,fold_epi_neg_pos,fold_epi_neg_neg,fold_epi_pvalue_NO_gene)

    list_type <- list("","_Epi"    ,"_Input","_MethvsGene","_Meth+Gene-","_Meth+Gene+"
                      ,"_Meth-Gene+","_Meth-Gene-","_MethNoGene")

    file_to_plot <- fold_epi
    file_to_plot$Color <- "Grey"
    file_to_plot[file_to_plot$ID %in% fold_epi_pvalue_gene$ID,"Color"] = "Red"
    file_to_plot[file_to_plot$ID %in% fold_epi_pvalue_NO_gene$ID,"Color"] = "Blue"
    png(filename=paste(figure_folder,"/",peak_name,"_ScatterPlot_",contr_Epi_names[i],".png",sep=""),width = 800, height = 800)
    plot(file_to_plot[,diffmethyl],file_to_plot[,diffgene],col = file_to_plot$Color,
         type = "p",main=diffmethyl, xlim=c(-8,8),ylim=c(-8,8),xlab = "logFC of Methyl sites in IP", ylab = "logFC of gene in Input", pch=19)
    dev.off()

    for(k in 1:length(list_type)){
      matrix <- list_matrix[[k]]
      rownames(matrix) = matrix$ID
      type <- list_type[[k]]
      norm.filename = paste(folderContrast,contr_Epi_names[i],type,".xls",sep = "")
      write.table(matrix, file=norm.filename, quote=FALSE, sep = "\t", row.names = FALSE)
      norm.filename = paste(folderContrast,contr_Epi_names[i],type,"_gene.txt",sep = "")
      write.table(unique(matrix$gene_name), file=norm.filename, quote=FALSE, row.names = FALSE, col.names = FALSE)
      norm.filename = paste(folderContrast,contr_Epi_names[i],type,"_uniprot.txt",sep = "")
      write.table(unique(matrix$UniprotID), file=norm.filename, quote=FALSE, row.names = FALSE, col.names = FALSE)
      bed_col <- c("chromo_peak","begin_peak","end_peak","ID","Motif","strand")
      if(type != "_Input" & nrow(matrix) > 0){
        #print("Save bed and fasta")
        # save in bedfile for GUITAR plot
        bed_file <- unique(matrix[bed_col])
        norm.filename = paste(folderBed,contr_Epi_names[i],type,".bed",sep = "")
        write.table(bed_file, file=norm.filename, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

          seq_fasta <- fasta[names(fasta) %in% matrix$ID]
          # save all peaks in fasta
          #print(paste("Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
          norm.filename = paste(folderFasta,contr_Epi_names[i],type,".fasta",sep = "")
          seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)

          # type_overlaps = c("3UTR","5UTR","CDS","Intron")
          # for(type_overlap in type_overlaps){
          #   matrix_type = matrix[grep(type_overlap,matrix$Type_overlap),]
          #   seq_fasta <- fasta[names(fasta) %in% matrix_type$ID]
          #   print(paste(type_overlap,"Extract:",length(seq_fasta),"from",nrow(matrix),"rows matrix"))
          #   norm.filename = paste(folderFasta,contr_Epi_names[i],type,"_",type_overlap,".fasta",sep = "")
          #   seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)
          # }


      }
      summary <- c(contr_Epi_names[i],type,length(rownames(matrix)),length(unique(matrix$gene_name)))
      #print(summary)
      summary_table <- rbind(summary_table,summary)
    }

    # draw motif distrib
    list_matrix <- list(fold_epi[best$ID,],fold_epi_pvalue_gene,fold_epi_pvalue_NO_gene)
    list_name <- list("Meth","MethGene","MethNoGene")
    for(k in 1:length(list_name)){
      matrix <- list_matrix[[k]]
      type <- list_name[[k]]
      # Put motif presence in a vector
      if(length(matrix$Motif)!=0){
        for(w in 1:length(matrix$Motif)){
          motif_distrib <- rbind(motif_distrib, matrix$Motif[w])
          motif_length <- rbind(motif_length, matrix$length_peak[w])
          motif_type <- rbind(motif_type, type)
          motif_name <- rbind(motif_name, contr_Epi_names[i])
        }
      }
    }
  }

  diff_peaks <- unique(diff_peaks)
  diff_genes <- unique(diff_genes)
  diff_peaks_genes <- unique(diff_peaks_genes)
  diff_peaks_no_genes <- diff_peaks %in% diff_peaks_genes
  diff_peaks_pos_neg <- unique(diff_peaks_pos_neg)
  diff_peaks_pos_pos <- unique(diff_peaks_pos_pos)
  diff_peaks_neg_neg <- unique(diff_peaks_neg_neg)
  diff_peaks_neg_pos <- unique(diff_peaks_neg_pos)
  diff_Epi <- unique(diff_Epi)


  write.table(summary_table, file=paste(peak_name,"_Diff.txt",sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  assign("diff_peaks",diff_peaks, pos = globalenv())
  assign("diff_genes",diff_genes, pos = globalenv())
  assign("diff_Epi",diff_Epi, pos = globalenv())
  assign("diff_peaks_genes",diff_peaks_genes, pos = globalenv())
  assign("diff_peaks_NOgenes",diff_peaks_no_genes, pos = globalenv())
  assign("diff_peaks_pos_neg",diff_peaks_pos_neg, pos = globalenv())
  assign("diff_peaks_pos_pos",diff_peaks_pos_pos, pos = globalenv())
  assign("diff_peaks_neg_neg",diff_peaks_neg_neg, pos = globalenv())
  assign("diff_peaks_neg_pos",diff_peaks_neg_pos, pos = globalenv())

}
