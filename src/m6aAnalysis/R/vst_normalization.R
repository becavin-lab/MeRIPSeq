vst_normalisation <- function(){

  # VST normalization for MeRIPSeq
  dds <- DESeqDataSetFromMatrix(countData = counts_epi, colData = target_epi, design = formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)), collapse = " "), "\n")
  dds <- estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  dds <- estimateDispersions(dds, fitType = fitType)
  dds <- nbinomWaldTest(dds)
  norm_vst_epi <- data.frame(assay(varianceStabilizingTransformation(dds)))

  # VST for RNASeq Input
  dds <- DESeqDataSetFromMatrix(countData = counts_RNA_I, colData = target_RNA_I, design = formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)), collapse = " "), "\n")
  dds <- estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  dds <- estimateDispersions(dds, fitType = fitType)
  dds <- nbinomWaldTest(dds)
  norm_vst_RNA_I <- data.frame(assay(varianceStabilizingTransformation(dds)))

  # VST for RNASeq IP
  dds <- DESeqDataSetFromMatrix(countData = counts_RNA_IP, colData = target_RNA_IP, design = formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)), collapse = " "), "\n")
  dds <- estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  dds <- estimateDispersions(dds, fitType = fitType)
  dds <- nbinomWaldTest(dds)
  norm_vst_RNA_IP <- data.frame(assay(varianceStabilizingTransformation(dds)))

  # assign variables
  norm_voom_epi <- norm_vst_epi
  norm_voom_RNA_I <- norm_vst_RNA_I
  norm_voom_RNA_IP <- norm_vst_RNA_IP

  assign("norm_voom_epi",norm_voom_epi,pos = globalenv())
  assign("norm_voom_RNA_I",norm_voom_RNA_I,pos = globalenv())
  assign("norm_voom_RNA_IP",norm_voom_RNA_IP,pos = globalenv())

}
