load_contrast <- function(){
  # Get all comparisons
  # You can create your own contrast HERE
  # contrasts = c("Cond2-Cond1","Cond3-Cond2")
  all_biocond <- unique(as.character(target_epi_fold$BioCond))
  contrasts <- combn(x = all_biocond, m = 2)
  contrasts <- paste(contrasts[2,],"-",contrasts[1,],sep="")
  if(exp_design_name == "CecAm"){
    contrasts <- c("CONV-GF", "CONV-abx", "CONV-ex_GF", "ex_GF-GF" , "ex_GF-abx" , "GF-abx", "CONV-vanco", "ex_GF-vanco", "vanco-abx", "vanco-GF", "CONV-Am", "ex_GF-Am", "GF-Am", "Am-abx", "vanco-Am")
  } else if (exp_design_name %in% c("Liver_2019", "Liver_all","LiverZT")) {
    # compare conditions within each time point
    contrasts <- c(sub('-', '-T3.', sub('^', 'T3.', contrasts)), sub('-', '-T13.', sub('^', 'T13.', contrasts)))
    # compare same condition across both time points
    contrasts <- c(contrasts, paste0('T13.',all_biocond,'-','T3.',all_biocond))
  }
  contr_Epi_names <- contrasts
  contr_RNA_I_names <- paste(contrasts, '_Input',sep= '')
  # # replace first constrast by appropriate names
  # if(exp_design_name == "CecAm"){
  #   contrasts <- c("CONV-GF", "CONV", "CONV-ex_GF", "ex_GF-GF" , "ex_GF" , "GF", "CONV-vanco", "ex_GF-vanco", "vanco", "vanco-GF", "CONV-Am", "ex_GF-Am", "GF-Am", "Am", "vanco-Am")
  # }else{
  #   for(index in 2:length(all_biocond)){
  #     contrasts[index-1] = all_biocond[index]
  #   }
  # }
  assign("contrasts",contrasts, pos = globalenv())
  assign("contr_Epi_names",contr_Epi_names, pos = globalenv())
  assign("contr_RNA_I_names",contr_RNA_I_names, pos = globalenv())
  assign("all_biocond",all_biocond, pos = globalenv())
}
