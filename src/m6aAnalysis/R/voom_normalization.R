voom_normalisation <- function(){



  # fit models
  print("Voom normalization MeRIPSeq")
  dge_epi <- edgeR::DGEList(counts=counts_epi)
  dge_epi <- edgeR::calcNormFactors(dge_epi)
  dge_epi_voom <- limma::voom(dge_epi, design = model_epi_all, plot=TRUE)#, normalize="quantile")
  norm_voom_epi_all <- data.frame(dge_epi_voom, check.names = FALSE)

  #dge_epi <- dge_epi[, rownames(model_epi)]
  if (!all(rownames(model_epi) %in% colnames(dge_epi)))
    subids <- sub('-','.',rownames(model_epi))
  else
    subids <- rownames(model_epi)

  #norm_voom_epi <- data.frame(limma::voom(dge_epi[, subids], design = model_epi, plot=TRUE))

  dge_epiIP <- edgeR::DGEList(counts=counts_epi[, grep('_IP', colnames(counts_epi))])
  dge_epiIP <- edgeR::calcNormFactors(dge_epiIP)
  dge_epiIP_voom <- limma::voom(dge_epiIP, design = model_epi, plot=TRUE)#, normalize="quantile")
  norm_voom_epi <- data.frame(dge_epiIP_voom, check.names = FALSE)
  #
  # dge_epiInput <- edgeR::DGEList(counts=counts_epi[, grep('_.nput', colnames(counts_epi))])
  # dge_epiInput <- edgeR::calcNormFactors(dge_epiInput)
  # norm_voom_epiInput <- data.frame(limma::voom(dge_epiInput, design = model_epi_all[grep('_.nput', colnames(counts_epi)),], plot=TRUE))
  #
  # Transcriptomics
  print("Voom normalization RNASeq ")
  dge_RNA_I <- edgeR::DGEList(counts=counts_RNA_I)
  dge_RNA_I <- edgeR::calcNormFactors(dge_RNA_I)
  norm_voom_RNA_I <- data.frame(limma::voom(dge_RNA_I, design = model_RNA_I, plot=TRUE), check.names = FALSE)

  dge_RNA_IP <- edgeR::DGEList(counts=counts_RNA_IP)
  dge_RNA_IP <- edgeR::calcNormFactors(dge_RNA_IP)
  norm_voom_RNA_IP <- data.frame(limma::voom(dge_RNA_IP, design = model_RNA_IP, plot=TRUE), check.names = FALSE)


  print("Fit Fold")
  norm_fold_epi <- data.frame(norm_voom_epi_all[,as.character(target_epi_fold$IP)] -
                                norm_voom_epi_all[,as.character(target_epi_fold$Input)], #match(tolower(sub('-T', '.T', as.character(target_epi_fold$Input))),tolower(colnames(norm_voom_epi_all)))
                              check.names = FALSE)
  #norm_fold_epi <- data.frame((norm_voom_epiIP[,as.character(target_epi_fold$IP)] - norm_voom_epiInput[,as.character(target_epi_fold$Input)]))
  colnames(norm_fold_epi) <- target_epi_fold$DataName

  # fit with lmfit
  fit_epi <- lmFit(norm_voom_epi, design=model_epi)
  fit_epi_fold <- lmFit(norm_fold_epi, design=model_epi_fold)
  fit_RNA_I <- lmFit(norm_voom_RNA_I, design=model_RNA_I)

  assign("fit_epi",fit_epi,pos = globalenv())
  assign("fit_epi_fold",fit_epi_fold,pos = globalenv())
  assign("fit_RNA_I",fit_RNA_I,pos = globalenv())

  assign("norm_voom_epi_all",norm_voom_epi_all,pos = globalenv())
  assign("norm_voom_epi",norm_voom_epi,pos = globalenv())
  assign("norm_fold_epi",norm_fold_epi,pos = globalenv())
  assign("norm_voom_RNA_I",norm_voom_RNA_I,pos = globalenv())
  assign("norm_voom_RNA_IP",norm_voom_RNA_IP,pos = globalenv())


}
