voom_normalisation_exon <- function(){

  # Transcriptomics
  print("Voom normalization RNASes transcript")
  dge_Transc_I <- DGEList(counts=counts_Trans_I)
  dge_Transc_I <- calcNormFactors(dge_Transc_I)
  norm_voom_Trans_I <- data.frame(voom(dge_Transc_I, design = model_RNA_I, plot=TRUE))

  print("Voom normalization RNASes exon")
  dge_Exon_I <- DGEList(counts=counts_Exon_I)
  dge_Exon_I <- calcNormFactors(dge_Exon_I)
  norm_voom_Exon_I <- data.frame(voom(dge_Exon_I, design = model_RNA_I, plot=TRUE))

  # fit with lmfit
  fit_Trans_I <- lmFit(norm_voom_Trans_I, design=model_RNA_I)
  fit_Exon_I <- lmFit(norm_voom_Exon_I, design=model_RNA_I)

  assign("fit_Trans_I",fit_Trans_I,pos = globalenv())
  assign("fit_Exon_I",fit_Exon_I,pos = globalenv())

  assign("norm_voom_Trans_I",norm_voom_Trans_I,pos = globalenv())
  assign("norm_voom_Exon_I",norm_voom_Exon_I,pos = globalenv())

}
