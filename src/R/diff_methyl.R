###############################################################################
### R script to compare several MERIPSeq conditions with the SARTools packages
### Christophe Becavin, Hugo Varet, Anne Biton
### April 20th, 2018
################################################################################
################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
# local variables
#general_path <- "/Volumes/m6aAkker/"
# serveur variables
general_path <- "/pasteur/projets/policy01/m6aAkker/"

rscript_folder <- paste(general_path,"MeRIPSeq/src/m6aAnalysis/",sep="")
pkgload::load_all(path = rscript_folder, reset = TRUE)

args <- commandArgs(trailingOnly = TRUE)
exp_design_name <- args[1]
bed_name <- args[2]
peak_name <- paste(exp_design_name,"_",bed_name,sep="")

project_dir <- paste(general_path,"PeakDiffExpression/",peak_name,"/",sep="")
rdata_file <- paste(project_dir,peak_name,".RData",sep="")
# if(file.exists(rdata_file)){
#   print(paste("Load RData: ",rdata_file,sep = ""))
#   load(rdata_file)
# }else{

  annotation_name <- "gencode.vM13"
  exp_design_file = paste(project_dir,"/../../ExpDesign/",exp_design_name,"_exp_design.txt",sep="")

  print(paste("Run Rscript diff_methyl for",peak_name,"in",getwd()))

  featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                        "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                        "not_aligned", "too_low_aQual",
                        "N_noFeature","N_multimapping","N_unmapped","N_ambiguous")

  varInt <- c("BioCond")                                    # factor of interest
  batch <- c('dataset','Phase','Seq','Library','RNAprep')   # blocking factor: NULL (default) or "batch" for example
  locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors
  fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"

  vst_norm <- FALSE
  combat_all <- FALSE     # should Combat remove the batch to methyl site + Gene or separated
  combat_param <- TRUE    # Should Combat be run with paramtric or non parametric estimation
  batch_remove <- "Seq"   # list of batch to be removed by combat

  if(length(grep("Cec|Liv", exp_design_name)) > 0) {
    color_biocond <- c(Am_2019="seagreen4",Lp_2019="darkorange3",CONV_2019="gray0",GF_2019="dodgerblue",abx_2019='magenta',abxGF_2019='darkblue', Ec_2019="darkred",vanco_2019="DarkMagenta", 
                       abx_Liver_2019='magenta', Am_Liver_2019="seagreen4",Lp_Liver_2019="darkorange3",CONV_Liver_2019="gray0",GF_Liver_2019="dodgerblue",abxGF_Liver_2019='darkblue', 
                       Ec_Liver_2019="darkred",vanco_Liver_2019="DarkMagenta",
                       Am="seagreen4",Lp="darkorange3",CONV="gray0",ex_GF = 'gray51', GF="dodgerblue",abx='magenta',abxGF='darkblue',abx_GF='darkblue', Ec="darkred",vanco="DarkMagenta") #"gray51","tan4","tan3",
    color_biocond2 <-   color_biocond3 <-   color_biocond
    names(color_biocond2) <- paste0(names(color_biocond2), '.T', 13)
    names(color_biocond3) <- paste0(names(color_biocond3), '.T', 3)
    color_biocond <- c(color_biocond, color_biocond2, color_biocond3)
    
  } else {
    color_biocond <- c("seagreen4","seagreen3","gray0","gray51","tan4","tan3","firebrick1","DarkMagenta")
  }

  init_parameters()

  ##########################################
  #######      RUN ANALYSIS          #######
  ##########################################

  # load target, count, peaks, gene, annotation
  # include the design model you want to use in the load_target function for your experimental design
  load_target()

  load_peaks_genes()
  load_counts()

  # Plot some description figures to see the distributuon of raw counts in the datasets
  # Update : featuresToRemove in load_counts()
  plot_description()

  # VST normalization
  if(vst_norm){
    vst_normalisation()
    # Correct batch effect
    #voom_combat_normalization()
  }else{
    voom_normalisation()
  }

  # #Draw Clustering to verify data distribution
  #plot_clustering()

  # fit contrast
  load_contrast()
  calculate_contrast()
  save_contrast()
 
  # plot upset plots
  plot_upset('upset_plots.pdf')
 
  # Create diff summary
  diff_summary()

  # Extract peaks specific to a biocond and save in bed and fasta
  save_biocond_list()

  # Plot PCA and Heatmap of differentially methylated sites
  filter_table <- norm_voom_epi_all[diff_peaks,]
  plot_heatmap(filter_table,"")
  plot_pca(filter_table,"")

  #filter_table <- filter_table[,grep('_IP_',colnames(filter_table))]
  filter_table <- norm_voom_epi[diff_peaks,]
  plot_heatmap(filter_table,"_IP")
  plot_pca(filter_table,"_IP")

  filter_table <- norm_voom_epi_all[diff_peaks,]
  filter_table <- filter_table[,colnames(filter_table) %in% colnames(counts_RNA_I)]
  plot_heatmap(filter_table,"_Input")
  plot_pca(filter_table,"_Input")

  # filter_table <- data.frame(norm_voom_RNA_IP[diff_genes,])
  filter_table <- norm_voom_RNA_I[diff_genes,]
  plot_heatmap(filter_table,"_Gene")
  plot_pca(filter_table,"_Gene")

  #Plot PCA of the top 100 significantly differentially methylated sites with abs(logFC)>=1
  for (ntop in c(100)) {
    for (lfc in 1) {
    diff_peaks_lfc1 = unique(unlist(sapply(contr_Epi_names,
                                               function(i) {
                                                 head(rownames(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", 
                                                                               p.value = 0.05, lfc = lfc, number=Inf)),ntop)
                                               })))
  
    diff_genes_lfc1 = unique(unlist(sapply(sub('_Input','',contr_RNA_I_names),
                                           function(i) {
                                             head(rownames(limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", 
                                                                           p.value = 0.05, lfc = lfc, number=Inf)),ntop)
                                           })))
    
    filter_table <- data.frame(norm_voom_epi_all[diff_peaks_lfc1,])
    plot_heatmap(filter_table,paste0("_lfc",lfc,"_top",ntop))
    plot_pca(filter_table,paste0("_lfc",lfc,"_top",ntop))
    
    filter_table <- data.frame(norm_voom_epi[diff_peaks_lfc1,])
    plot_heatmap(filter_table,paste0("_IP_lfc", lfc,"_top",ntop))
    plot_pca(filter_table,paste0("_IP_lfc",lfc,"_top",ntop))

    filter_table <- data.frame(norm_voom_RNA_I[diff_genes_lfc1,])
    plot_heatmap(filter_table,paste0("_Input_Genes_lfc", lfc,"_top",ntop))
    plot_pca(filter_table,paste0("_Input_Genes_lfc",lfc,"_top",ntop))
    
    }
  }
  
  # correct data using limma::removeBatchEffect function for the plots
  # see advice of Gordon Smyth that i followed at https://support.bioconductor.org/p/83286/
  norm_voom_epi_allsave <- norm_voom_epi_all
  norm_voom_epi_allcorr <- removeBatchEffect(norm_voom_epi_all,  
                                             covariates = model_epi_all[, grep(paste0(batch,collapse='|'),colnames(model_epi_all),invert=F)], 
                                             design = model_epi_all[, setdiff(grep(paste0(c(levels(target_epi$BioCond)),collapse='|'),colnames(model_epi_all),invert=F), 
                                                                              c(which(colnames(model_epi_all) == "IP")))]) 
  norm_voom_epi_all <- norm_voom_epi_allcorr
  
  norm_voom_episave <- norm_voom_epi
  norm_voom_epicorr <- removeBatchEffect(norm_voom_epi, 
                                         covariates = model_epi[, grep(paste0(unique(target_epi$BioCond),collapse='|'),colnames(model_epi),invert=T)], 
                                         design = model_epi[, grep(paste0(unique(target_epi$BioCond),collapse='|'),colnames(model_epi),invert=F)])
  norm_voom_epi <- norm_voom_epicorr
  
  norm_voom_RNA_Isave <- norm_voom_RNA_I
  norm_voom_RNA_Icorr <- removeBatchEffect(norm_voom_RNA_I, 
                                           covariates = model_RNA_I[, grep(paste0(unique(target_RNA_I$BioCond),collapse='|'),colnames(model_RNA_I),invert=T)], 
                                           design = model_RNA_I[, grep(paste0(unique(target_RNA_I$BioCond),collapse='|'),colnames(model_RNA_I),invert=F)])
  norm_voom_RNA_I <- norm_voom_RNA_Icorr
  
  # Plot PCA and Heatmap of differentially methylated sites
  filter_table <- norm_voom_epi_all[diff_peaks,]
  plot_heatmap(filter_table,"Corrected")
  plot_pca(filter_table,"Corrected")
  
  #filter_table <- filter_table[,grep('_IP_',colnames(filter_table))]
  filter_table <- norm_voom_epi[diff_peaks,]
  plot_heatmap(filter_table,"_IPcorrected")
  plot_heatmap(filter_table,"_IPcorrected_scale", scaleRow=TRUE)
  plot_pca(filter_table,"_IPcorrected")
  
  filter_table <- norm_voom_epi_all[diff_peaks,]
  filter_table <- filter_table[,colnames(filter_table) %in% colnames(counts_RNA_I)]
  plot_heatmap(filter_table,"_InputCorrected")
  plot_heatmap(filter_table,"_InputCorrected_scale", scaleRow=TRUE)
  plot_pca(filter_table,"_InputCorrected")
  
  filter_table <- norm_voom_RNA_I[diff_genes,]
  plot_heatmap(filter_table,"_GeneCorrected")
  plot_pca(filter_table,"_GeneCorrected")
  
  #Plot PCA of the top 100 significant differentially methylated sites with abs(logFC)>=1
  for (ntop in c(100)) {
    for (lfc in 1) {
      diff_peaks_lfc1 = unique(unlist(sapply(contr_Epi_names,
                                             function(i) {
                                               head(rownames(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", 
                                                                             p.value = 0.05, lfc = lfc, number=Inf)),ntop)
                                             })))
      
      diff_genes_lfc1 = unique(unlist(sapply(sub('_Input','',contr_RNA_I_names),
                                             function(i) {
                                               head(rownames(limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", 
                                                                             p.value = 0.05, lfc = lfc, number=Inf)),ntop)
                                             })))
      
      filter_table <- data.frame(norm_voom_epi_all[diff_peaks_lfc1,])
      plot_heatmap(filter_table,paste0("Corrected_lfc",lfc,"_top",ntop))
      plot_pca(filter_table,paste0("Corrected_lfc",lfc,"_top",ntop))
      
      filter_table <- data.frame(norm_voom_epi[diff_peaks_lfc1,])
      plot_heatmap(filter_table,paste0("_IPcorrected_lfc", lfc,"_top",ntop))
      plot_pca(filter_table,paste0("_IPcorrected_lfc",lfc,"_top",ntop))

      filter_table <- data.frame(norm_voom_RNA_I[diff_genes_lfc1,])
      plot_heatmap(filter_table,paste0("_InputCorrected_Genes_lfc", lfc,"_top",ntop))
      plot_pca(filter_table,paste0("_InputCorrected_Genes_lfc",lfc,"_top",ntop))

    }
  }
  
  norm_voom_RNA_I <- norm_voom_RNA_Isave
  norm_voom_epi_all <- norm_voom_epi_allsave
  norm_voom_epi <- norm_voom_episave
  
  
  print("Save RData")
  print(rdata_file)
  save.image(rdata_file)
#}


## Search for transcript hyper and hypo methylated
# for(contr_names in contr_Epi_names){
#   diff_array <- read.delim(file = paste("DiffMethyl/",contr_names,"/",contr_names,".xls",sep=""),check.names = F)
#   hypo_dfarray <- diff_array[diff_array[,contr_names]<0,]
#   hyper_dfarray <- diff_array[diff_array[,contr_names]>0,]
#   gene_hypo_hyper <- Reduce(intersect,  list(hypo_dfarray$Gene_ID,hyper_dfarray$Gene_ID))
#   methyl_hypo_or_hyper <- Reduce(union,  list(hypo_dfarray$ID,hyper_dfarray$ID))
#   #write.table(gene_hypo_hyper, paste("DiffMethyl/",contr_names,"_hyper_hypo.txt",sep=""))
#   print(paste("Diff methyl for",contr_names," genes with hypo meth sites",length(hypo_dfarray$Gene_ID)))
#   print(paste("Diff methyl for",contr_names," genes with hyper meth sites",length(hyper_dfarray$Gene_ID)))
#   print(paste("Diff methyl for",contr_names," genes with hyper and hypo meth sites",length(gene_hypo_hyper)))
#   if(length(gene_hypo_hyper)>0){
#     filter_diff_array = diff_array[diff_array$ID %in% methyl_hypo_or_hyper,]
#     filter_diff_array = filter_diff_array[filter_diff_array$Gene_ID %in% gene_hypo_hyper,]
#     write.table(filter_diff_array, paste("DiffMethyl/",contr_names,"_diff_hyper_hypo.txt",sep=""), quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
#   }
# }
#
#
#
# print("Run save_biocond_list")
# # # # Play with gene lists
# bioconds <- unique(as.character(target_epi$BioCond))
# data_to_bioconds <- list()
# norm_median_epi <- norm_voom_epi
# for(i in 1:length(bioconds)){
#   data_biocond <- as.character(target_epi$labels[target_epi$BioCond==bioconds[i] & target_epi$TypeData=="IP"])
#   data_to_bioconds[[bioconds[i]]] = data_biocond
#   norm_median_epi[,bioconds[i]] = apply(norm_median_epi[,data_biocond], 1, function(x) median(x))
# }
#
# norm_median_epi <- norm_median_epi[,bioconds]
# norm_median_epi$ID <- rownames(norm_median_epi)
# norm_median_epi <- merge(norm_median_epi, peaks, by.all = "ID")
# rownames(norm_median_epi) <-norm_median_epi$ID
#
# cut_off = 2
# for(biocond in bioconds){
#   hypo_dfarray <- norm_median_epi[norm_median_epi[,biocond] < (-cut_off),]
#   hyper_dfarray <- norm_median_epi[norm_median_epi[,biocond] > cut_off,]
#   gene_hypo_hyper <- Reduce(intersect,  list(hypo_dfarray$Gene_ID,hyper_dfarray$Gene_ID))
#   methyl_hypo_or_hyper <- Reduce(union,  list(hypo_dfarray$ID,hyper_dfarray$ID))
#   #write.table(gene_hypo_hyper, paste("DiffMethyl/",contr_names,"_hyper_hypo.txt",sep=""))
#   print(paste(biocond,"genes with methyl sites hypo",length(hypo_dfarray$Gene_ID)))
#   print(paste(biocond,"genes with methyl sites hyper",length(hyper_dfarray$Gene_ID)))
#   print(paste(biocond,"gene with methyl sites hyper and hypo",length(gene_hypo_hyper)))
#   filter_norm_median_epi = norm_median_epi[norm_median_epi$ID %in% methyl_hypo_or_hyper,]
#   filter_norm_median_epi = filter_norm_median_epi[filter_norm_median_epi$Gene_ID %in% gene_hypo_hyper,]
#   write.table(filter_norm_median_epi, paste("DiffMethyl/",biocond,"_all_hyper_hypo.txt",sep=""), quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# }
#
# list_peaks <- diff_peaks
# norm_median_epi_filter <- norm_median_epi[list_peaks,]
# for(biocond in bioconds){
#    hypo_dfarray <- norm_median_epi_filter[norm_median_epi_filter[,biocond] < (-cut_off),]
#    hyper_dfarray <- norm_median_epi_filter[norm_median_epi_filter[,biocond] > cut_off,]
#    gene_hypo_hyper <- Reduce(intersect,  list(hypo_dfarray$Gene_ID,hyper_dfarray$Gene_ID))
#    methyl_hypo_or_hyper <- Reduce(union,  list(hypo_dfarray$ID,hyper_dfarray$ID))
#    print(paste(biocond,"genes with diff methyl sites hypo",length(hypo_dfarray$Gene_ID)))
#    print(paste(biocond,"genes with diff methyl sites hyper",length(hyper_dfarray$Gene_ID)))
#    print(paste(biocond,"gene with diff methyl sites hyper and hypo",length(gene_hypo_hyper)))
#    filter_norm_median_epi = norm_median_epi_filter[norm_median_epi_filter$ID %in% methyl_hypo_or_hyper,]
#    filter_norm_median_epi = filter_norm_median_epi[filter_norm_median_epi$Gene_ID %in% gene_hypo_hyper,]
#    write.table(filter_norm_median_epi, paste("DiffMethyl/",biocond,"_Meth_hyper_hypo.txt",sep=""), quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# }
