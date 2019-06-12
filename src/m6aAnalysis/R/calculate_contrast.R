calculate_contrast <- function(){
  # Fit contrast
  print("Fit contrast")
  print(contrasts)
  contr_matrix_epi <- limma::makeContrasts(contrasts = contrasts, levels=colnames(model_epi))
  vfit_epi <- limma::contrasts.fit(fit_epi, contrasts=contr_matrix_epi)
  fit_bayes_epi <- limma::eBayes(vfit_epi, robust = TRUE)
  epi <- data.frame(fit_bayes_epi$coefficients)
  colnames(epi) = contr_Epi_names
  epi$ID = rownames(epi)

  contr_matrix_epi_fold <- limma::makeContrasts(contrasts = contrasts, levels=colnames(model_epi_fold))
  vfit_epi_fold <- limma::contrasts.fit(fit_epi_fold, contrasts=contr_matrix_epi_fold)
  fit_bayes_epi_fold <- limma::eBayes(vfit_epi_fold, robust = TRUE)
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
  contr_matrix_RNA_I <- limma::makeContrasts(contrasts = contrasts, levels=colnames(model_RNA_I))
  vfit_RNA_I <- limma::contrasts.fit(fit_RNA_I, contrasts=contr_matrix_RNA_I)
  fit_bayes_RNA_I <- limma::eBayes(vfit_RNA_I, robust = TRUE)
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

  # fold_all <- merge(fold_Exon,fold_Trans,by = "Transcript_ID", all = T)
  # fold_all <- merge(fold_epi, fold_all, by="Gene_ID", all =T)

  norm_voom_epi$ID <- NULL
  norm_voom_epi$ID <- rownames(norm_voom_epi)
  count_epi <- merge(norm_voom_epi,peaks,by.all="ID")
  norm_voom_RNA_I$ID <- NULL
  norm_voom_RNA_I$ID <- rownames(norm_voom_RNA_I)
  count_RNA <- merge(norm_voom_RNA_I,genes,by.x="ID",by.y="Gene_ID")
  write.table(fold_epi, file=paste(peak_name,"_Epi_Fold.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(count_epi, file=paste(peak_name,"_Epi_Count.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(fold_RNA, file=paste(peak_name,"_RNA_Fold.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(count_RNA, file=paste(peak_name,"_RNA_Count.xls",sep=""), quote=FALSE, row.names = FALSE, sep = "\t")

  # plot histograms of p-values
  pdf(paste0(figure_folder,peak_name,"_Epi_Fold_pvalues.pdf"), 10,10)
  par(mfrow = c(4,4))
  for(i in 1:ncol(fit_bayes_epi_fold$p.value)) hist(fit_bayes_epi_fold$p.value[,i], breaks = 70, main = colnames(fit_bayes_epi_fold$p.value)[i])
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_Epi_Count_pvalues.pdf"), 10,10)
  par(mfrow = c(4,4))
  for(i in 1:ncol(fit_bayes_epi$p.value)) hist(fit_bayes_epi$p.value[,i], breaks = 70, main = colnames(fit_bayes_epi$p.value)[i])
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_RNA_Input_pvalues.pdf"), 10,10)
  par(mfrow = c(4,4))
  for(i in 1:ncol(fit_bayes_RNA_I$p.value)) hist(fit_bayes_RNA_I$p.value[,i], breaks = 70, main = colnames(fit_bayes_RNA_I$p.value)[i])
  dev.off()

  pdf(paste0(figure_folder, peak_name,"_Epi_Count_nbDiffPeaks.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_Epi_names, function(i) nrow(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 0.05, lfc = 1, number=Inf))), names = contr_Epi_names, las = 2)# , ylim = c(0,800))
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_Epi_Fold_nbDiffPeaks.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_Epi_names, function(i) nrow(limma::topTable(fit_bayes_epi_fold, coef = i, adjust="BH", p.value = 0.05, lfc = 1, number=Inf))), names = contr_Epi_names, las = 2)
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_RNA_Input_nbDiffGenes.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_RNA_I_names, function(i) nrow(limma::topTable(fit_bayes_RNA_I, coef = gsub('_Input', '', i), adjust="BH", p.value = 0.05, lfc = 1, number=Inf))), names = contr_Epi_names, las = 2)
  dev.off()

  pdf(paste0(figure_folder, peak_name,"_Epi_Count_nbDiffPeaks_lfc2.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_Epi_names, function(i) nrow(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 0.05, lfc = 2, number=Inf))), names = contr_Epi_names, las = 2)# , ylim = c(0,800))
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_Epi_Fold_nbDiffPeaks_lfc2.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_Epi_names, function(i) nrow(limma::topTable(fit_bayes_epi_fold, coef = i, adjust="BH", p.value = 0.05, lfc = 2, number=Inf))), names = contr_Epi_names, las = 2)
  dev.off()
  pdf(paste0(figure_folder, peak_name,"_RNA_Input_nbDiffGenes_lfc2.pdf"), 10,10)
  par(mar = c(20,4,4,4));
  barplot(sapply(contr_RNA_I_names, function(i) nrow(limma::topTable(fit_bayes_RNA_I, coef = gsub('_Input', '', i), adjust="BH", p.value = 0.05, lfc = 2, number=Inf))), names = contr_Epi_names, las = 2)
  dev.off()

  # tfit <- treat(fit_bayes_epi, lfc=1)#, robust = TRUE)
  # dt <- decideTests(tfit)
  # summary(dt)
  # for (contr in contr_Epi_names) {
  #   plotMD(tfit, column=contr, status=dt[,contr], )
  # }

  assign("fold_epi",fold_epi, pos = globalenv())
  assign("fold_RNA",fold_RNA, pos = globalenv())
  assign("contr_Epi_names",contr_Epi_names, pos = globalenv())
  assign("fit_bayes_epi",fit_bayes_epi, pos = globalenv())
  assign("fit_bayes_epi_fold",fit_bayes_epi_fold, pos = globalenv())
  assign("fit_bayes_RNA_I",fit_bayes_RNA_I, pos = globalenv())

}
