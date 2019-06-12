calculate_contrast_exon <- function(){

  # Get all comparisons
  # You can create your own contrast HERE
  # contrasts = c("Cond2-Cond1","Cond3-Cond2")
  all_biocond <- unique(as.character(target_epi_fold$BioCond))
  contrasts <- combn(x = all_biocond, m = 2)
  contrasts <- paste(contrasts[2,],"-",contrasts[1,],sep="")
  print(exp_design_name == "CecAm")
  if(exp_design_name == "CecAm"){
    contrasts <- c("CONV-GF", "CONV-abx", "CONV-ex_GF", "ex_GF-GF" , "ex_GF-abx" , "GF-abx", "CONV-vanco", "ex_GF-vanco", "vanco-abx", "vanco-GF", "CONV-Am", "ex_GF-Am", "GF-Am", "Am-abx", "vanco-Am")
  }
  contr_Epi_names <- contrasts
  contr_RNA_I_names <- paste(contrasts,"_Input",sep="")
  # replace first constrast by appropriate names
  if(exp_design_name == "CecAm"){
    contrasts <- c("CONV-GF", "CONV", "CONV-ex_GF", "ex_GF-GF" , "ex_GF" , "GF", "CONV-vanco", "ex_GF-vanco", "vanco", "vanco-GF", "CONV-Am", "ex_GF-Am", "GF-Am", "Am", "vanco-Am")
  }else{
    for(index in 2:length(all_biocond)){
      contrasts[index-1] = all_biocond[index]
    }
  }

  # Fit contrast
  print("Fit contrast")
  print(contrasts)
  contr_matrix_epi <- makeContrasts(contrasts = contrasts, levels=colnames(model_epi))
  vfit_epi <- contrasts.fit(fit_epi, contrasts=contr_matrix_epi)
  fit_bayes_epi <- eBayes(vfit_epi)
  epi <- data.frame(fit_bayes_epi$coefficients)
  colnames(epi) = contr_Epi_names
  epi$ID = rownames(epi)

  contr_matrix_epi_fold <- makeContrasts(contrasts = contrasts, levels=colnames(model_epi_fold))
  vfit_epi_fold <- contrasts.fit(fit_epi_fold, contrasts=contr_matrix_epi_fold)
  fit_bayes_epi_fold <- eBayes(vfit_epi_fold)
  fold_epi <- data.frame(fit_bayes_epi_fold$coefficients)
  colnames(fold_epi) = paste(contr_Epi_names,"_epi",sep="")
  fold_epi$ID = rownames(fold_epi)

  # create big table will all logfc and pvalue
  epi_pvalue <- data.frame(fit_bayes_epi$p.value)
  epi_adjpvalue <- data.frame(sapply(epi_pvalue, function(col) p.adjust(col, "BH")))
  colnames(epi_adjpvalue) = paste(contr_Epi_names,"_pval",sep="")
  rownames(epi_adjpvalue) <- rownames(epi_pvalue)
  epi_adjpvalue$ID = rownames(epi_adjpvalue)
  epi_merge <- merge(epi,epi_adjpvalue, by="ID")

  fold_epi_pvalue <- data.frame(fit_bayes_epi_fold$p.value)
  fold_epi_adjpvalue <- data.frame(sapply(fold_epi_pvalue, function(col) p.adjust(col, "BH")))
  rownames(fold_epi_adjpvalue) <- rownames(fold_epi_pvalue)
  colnames(fold_epi_adjpvalue) = paste(contr_Epi_names,"_epi_pval",sep="")
  fold_epi_adjpvalue$ID <- rownames(fold_epi_adjpvalue)

  fold_epi_merge <- merge(fold_epi,fold_epi_adjpvalue, by="ID")
  table_fold_pvalue <- merge(epi_merge,fold_epi_merge, by="ID")
  table_fold_pvalue <- merge(table_fold_pvalue, peaks_to_gene, by = "ID")

  # make contrast for Input
  # calc contrast
  contr_matrix_RNA_I <- makeContrasts(contrasts = contrasts, levels=colnames(model_RNA_I))
  vfit_RNA_I <- contrasts.fit(fit_RNA_I, contrasts=contr_matrix_RNA_I)
  fit_bayes_RNA_I <- eBayes(vfit_RNA_I)
  # get logFC and pvalue
  fold_RNA_I <- data.frame(fit_bayes_RNA_I$coefficients)
  colnames(fold_RNA_I) = contr_RNA_I_names
  fold_RNA_I$Gene_ID = rownames(fold_RNA_I)
  fold_RNA_I_pvalue <- data.frame(fit_bayes_RNA_I$p.value)
  fold_RNA_I_adjpvalue <- data.frame(sapply(fold_RNA_I_pvalue, function(col) p.adjust(col, "BH")))
  colnames(fold_RNA_I_adjpvalue) = paste(contr_RNA_I_names,"_pval",sep="")
  rownames(fold_RNA_I_adjpvalue) <- rownames(fold_RNA_I)
  fold_RNA_I_adjpvalue$Gene_ID = rownames(fold_RNA_I_adjpvalue)
  fold_RNA <- merge(fold_RNA_I,fold_RNA_I_adjpvalue, by="Gene_ID")

  # make contrast for Transcript
  # calc contrast
  contr_matrix_Trans_I <- makeContrasts(contrasts = contrasts, levels=colnames(model_RNA_I))
  vfit_Trans_I <- contrasts.fit(fit_Trans_I, contrasts=contr_matrix_Trans_I)
  fit_bayes_Trans_I <- eBayes(vfit_Trans_I)
  # get logFC and pvalue
  fold_Trans_I <- data.frame(fit_bayes_Trans_I$coefficients)
  colnames(fold_Trans_I) = contr_RNA_I_names
  fold_Trans_I$Transcript_ID = rownames(fold_Trans_I)
  fold_Trans_I_pvalue <- data.frame(fit_bayes_Trans_I$p.value)
  fold_Trans_I_adjpvalue <- data.frame(sapply(fold_Trans_I_pvalue, function(col) p.adjust(col, "BH")))
  colnames(fold_Trans_I_adjpvalue) = paste(contr_RNA_I_names,"_pval",sep="")
  rownames(fold_Trans_I_adjpvalue) <- rownames(fold_Trans_I)
  fold_Trans_I_adjpvalue$Transcript_ID = rownames(fold_Trans_I_adjpvalue)
  fold_Trans <- merge(fold_Trans_I,fold_Trans_I_adjpvalue, by="Transcript_ID")

  # make contrast for Transcript
  # calc contrast
  contr_matrix_Exon_I <- makeContrasts(contrasts = contrasts, levels=colnames(model_RNA_I))
  vfit_Exon_I <- contrasts.fit(fit_Exon_I, contrasts=contr_matrix_Exon_I)
  fit_bayes_Exon_I <- eBayes(vfit_Exon_I)
  # get logFC and pvalue
  fold_Exon_I <- data.frame(fit_bayes_Exon_I$coefficients)
  colnames(fold_Exon_I) = contr_RNA_I_names
  fold_Exon_I$Exon_ID = rownames(fold_Exon_I)
  fold_Exon_I_pvalue <- data.frame(fit_bayes_Exon_I$p.value)
  fold_Exon_I_adjpvalue <- data.frame(sapply(fold_Exon_I_pvalue, function(col) p.adjust(col, "BH")))
  colnames(fold_Exon_I_adjpvalue) = paste(contr_RNA_I_names,"_pval",sep="")
  rownames(fold_Exon_I_adjpvalue) <- rownames(fold_Exon_I)
  fold_Exon_I_adjpvalue$Exon_ID = rownames(fold_Exon_I_adjpvalue)
  fold_Exon <- merge(fold_Exon_I,fold_Exon_I_adjpvalue, by="Exon_ID")

  # merge all tables together
  print("Merge all tables")
  #fold_epi <- merge(epi,fold_epi, by="ID")
  fold_epi <- merge(table_fold_pvalue,fold_RNA, by="Gene_ID", all.x = TRUE)
  fold_epi <- merge(fold_epi, peaks, by="ID", all.x = TRUE)
  fold_epi$Gene_ID.y <-NULL
  fold_epi$Gene_ID <- fold_epi$Gene_ID.x
  fold_epi$Gene_ID.x <- NULL
  rownames(fold_epi) = fold_epi$ID
  fold_RNA <- merge(fold_RNA,genes,by.all="Gene_ID")
  rownames(fold_RNA) = fold_RNA$Gene_ID
  fold_Exon <- merge(fold_Exon,exons, by="Exon_ID", all = T)
  colnames(fold_Exon) <- gsub("_Input","_Exon",colnames(fold_Exon))
    fold_Trans <- merge(fold_Trans,transcripts, by="Transcript_ID", all = T)
  colnames(fold_Trans) <- gsub("_Input","_Trans",colnames(fold_Trans))

  # fold_all <- merge(fold_Exon,fold_Trans,by = "Transcript_ID", all = T)
  # fold_all <- merge(fold_epi, fold_all, by="Gene_ID", all =T)

  norm_voom_epi$ID <- NULL
  norm_voom_epi$ID <- rownames(norm_voom_epi)
  count_epi <- merge(norm_voom_epi,peaks,by.all="ID")
  norm_voom_RNA_I$ID <- NULL
  norm_voom_RNA_I$ID <- rownames(norm_voom_RNA_I)
  count_RNA <- merge(norm_voom_RNA_I,genes,by.x="ID",by.y="Gene_ID")
  write.table(fold_epi, file=paste(project_name,"_Epi_Fold.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(count_epi, file=paste(project_name,"_Epi_Count.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(fold_RNA, file=paste(project_name,"_RNA_Fold.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(count_RNA, file=paste(project_name,"_RNA_Count.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")

  assign("fold_epi",fold_epi, pos = globalenv())
  assign("fold_RNA",fold_RNA, pos = globalenv())
  assign("contr_Epi_names",contr_Epi_names, pos = globalenv())
  assign("fit_bayes_epi",fit_bayes_epi, pos = globalenv())
  assign("fit_bayes_epi_fold",fit_bayes_epi_fold, pos = globalenv())
  assign("fit_bayes_RNA_I",fit_bayes_RNA_I, pos = globalenv())

}
