voom_combat_normalization <- function(){


  # run ComBat batch correction
  if(combat_param){
    par_prior = TRUE
    mean_only = FALSE
  }else{
    par_prior = FALSE
    mean_only = TRUE
  }

  if(combat_all){
    # merge tables before running ComBat
    norm_voom_RNA = norm_voom_RNA_I
    norm_voom_RNA[,colnames(norm_voom_RNA_IP)] = norm_voom_RNA_IP
    norm_voom = norm_voom_epi
    norm_voom[rownames(norm_voom_RNA),] = norm_voom_RNA

    #SD = apply(norm_voom)
    SD = apply(norm_voom,1, sd, na.rm = TRUE)
    norm_voom$SD <- SD
    norm_voom <- norm_voom[norm_voom$SD > 0,]
    norm_voom$SD <- NULL

    for(i in 1:length(batch_remove)){
      batch_temp = batch_remove[[i]]
      print(paste("Remove batch: ",batch_temp," to combined count matrix",sep=""))
      norm_voom <- ComBat(dat=as.matrix(norm_voom), batch=as.numeric(target_epi[,batch_temp]), mod = model_epi, ref.batch = 1, par.prior = par_prior
                                , prior.plots = FALSE, mean.only = mean_only, BPPARAM = bpparam("SnowParam"))
      norm_voom <- data.frame(norm_voom)
    }
    norm_voom_epi <- norm_voom[rownames(norm_voom) %in% rownames(peaks),]
    norm_voom_RNA_I <- norm_voom[rownames(norm_voom) %in% genes$Gene_ID,]
    norm_voom_RNA_I <- norm_voom_RNA_I[,colnames(norm_voom) %in% colnames(counts_RNA_I),]

  }else{
    # for ComBat reduction only
    # merge tables before running ComBat
    SD = apply(norm_voom_epi,1, sd, na.rm = TRUE)
    norm_voom_epi$SD <- SD
    norm_voom_epi <- norm_voom_epi[norm_voom_epi$SD > 0,]
    norm_voom_epi$SD <- NULL
    SD = apply(norm_voom_RNA_I,1, sd, na.rm = TRUE)
    norm_voom_RNA_I$SD <- SD
    norm_voom_RNA_I <- norm_voom_RNA_I[norm_voom_RNA_I$SD > 0,]
    norm_voom_RNA_I$SD <- NULL

    for(i in 1:length(batch_remove)){
      batch_temp = batch_remove[[i]] #1 replaced by i by anne
      print(paste("Remove batch: ",batch_temp," to methyl site count matrix",sep=""))
      norm_voom_epi <- ComBat(dat=as.matrix(norm_voom_epi), batch=as.numeric(target_epi[,batch_temp]), mod = model_epi, ref.batch = 1, par.prior = par_prior
                          , prior.plots = FALSE, mean.only = mean_only, BPPARAM = bpparam("SnowParam"))
      norm_voom_epi <- data.frame(norm_voom_epi)
      print(paste("Remove batch: ",batch_temp," to gene count matrix",sep=""))
      norm_voom_RNA_I <- ComBat(dat=as.matrix(norm_voom_RNA_I), batch=as.numeric(target_RNA_I[,batch_temp]), mod = model_RNA_I, ref.batch = 1, par.prior = par_prior
                              , prior.plots = FALSE, mean.only = mean_only, BPPARAM = bpparam("SnowParam"))
      norm_voom_RNA_I <- data.frame(norm_voom_RNA_I)
    }
  }

  # get fold matrix  IP vs Input
  print("Fit Fold")
  print(head(target_epi_fold))
  norm_fold_epi <- data.frame((norm_voom_epi[,as.character(target_epi_fold$IP)] - norm_voom_epi[,as.character(target_epi_fold$Input)]))
  colnames(norm_fold_epi) <- target_epi_fold$DataName

  # fit with lmfit
  fit_epi <- lmFit(norm_voom_epi, design=model_epi)
  fit_epi_fold <- lmFit(norm_fold_epi, design=model_epi_fold)
  fit_RNA_I <- lmFit(norm_voom_RNA_I, design=model_RNA_I)

  assign("norm_voom_epi",norm_voom_epi,pos = globalenv())
  assign("fit_epi",fit_epi,pos = globalenv())

  assign("norm_fold_epi",norm_fold_epi,pos = globalenv())
  assign("fit_epi_fold",fit_epi_fold,pos = globalenv())

  assign("norm_voom_RNA_I",norm_voom_RNA_I,pos = globalenv())
  assign("fit_RNA_I",fit_RNA_I,pos = globalenv())

}
