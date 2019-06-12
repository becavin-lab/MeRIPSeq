plot_pca <- function(filter_table,type,scale=FALSE){
  print("Plot PCA")
  print(nrow(filter_table))

  col_names <- colnames(filter_table)
  #filter_table = t(scale(t(filter_table), center = TRUE))
  rv = apply(filter_table, 1, var, na.rm = TRUE)
  # if(nrow(filter_table)>500){
  #   pca = prcomp(t(filter_table[order(rv, decreasing = TRUE),][1:500, ]))
  # }else{
    pca = prcomp(t(filter_table), scale. = scale)
  # }

  if (all(unique(as.character(target_epi$BioCond)) %in% names(color_biocond))) {
      colsidecolors <- color_biocond[as.character(target_epi[colnames(filter_table),'BioCond'])]
      #col_labels <- target_epi[colnames(filter_table),'new_name']
  } else {
      target_epi$BioCond <- as.factor(target_epi$BioCond)
      colsidecolors <- color_biocond[as.integer(target_epi[colnames(filter_table),'BioCond'])]
      #col_labels <- target_epi[colnames(filter_table),'new_name']
  }

  prp <- pca$sdev^2 * 100/sum(pca$sdev^2)
  prp <- round(prp[1:3], 2)
  png(filename=paste(figure_folder,"PCA_NoLabel_",peak_name,type,".png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- ggplot2::autoplot(pca,label = FALSE, size=4, colour = colsidecolors,
                              xlabs = paste0("PC1 (",prp[1], "%)"), ylabs = paste0("PC2 (", prp[2], "%)"),
                              main = "Principal Component Analysis - Axes 1 and 2")
  print(figure)
  dev.off()

  if (FALSE) {
  ## anne, to be removed when solved
    edfinal2$RNAprepdate <-   gsub('UTC', '', edfinal2$RNAprepdate)
    edfinal2$RNAprepdate_monthyear <-   gsub('-..$|UTC', '', edfinal2$RNAprepdate)
    edfinal2$RNAprepdate_year <-   gsub('-..-..$|UTC', '', edfinal2$RNAprepdate)
    edfinal2$mRNAprepdate <-   gsub(' .*', '', edfinal2$mRNAprepdate)
    edfinal2$mRNAprepdate_monthyear <-   gsub('-..$|UTC', '', edfinal2$mRNAprepdate)
    edfinal2$mRNAprepdate_year <-   gsub('-..-..$|UTC', '', edfinal2$mRNAprepdate)
    edfinal2$fragmentation <- gsub(' .*', '', edfinal2$fragmentation)
    edfinal2$fragmentation_monthyear <- gsub('-..$', '', gsub(' .*', '', edfinal2$fragmentation))
    edfinal2$fragmentation_year <- gsub('-..-..$', '', gsub(' .*', '', edfinal2$fragmentation))
    edfinal2$qbit_date <- gsub(' .*', '', edfinal2$qbit)
    edfinal2$qbit_monthyear <- gsub('-..$', '', gsub(' .*', '', edfinal2$qbit))
    edfinal2$qbit_year <- gsub('-..-..$', '', gsub(' .*', '',edfinal2$qbit))
    edfinal2$cage <- c("5G" = 'G', "6D" = 'D', "5" = 'X', "6" = 'X', "abxGF1" = 'abxGF', "abxGF2" = 'abxGF', "abxGF3" = 'abxGF', "abxGF4" = 'abxGF',
      "Am1" = 'Am', "Am10" = 'Am', "Am5" = 'Am', "2" = 'X', "C1" = 'C', "1D" = "D", "C2" = 'C', "C5" = 'C', "Ec5" = "Ec", "Ec6" = "Ec",
      "EC7" = "Ec", "GF3"= 'GF', "GF8"= 'GF', "GF1"= 'GF', "GF6"= 'GF', "Lp8"= 'Lp', "Lp9" = "Lp", "Lp3" = "Lp", "3G" = "G",
      "4G" = 'G', "3DD" = 'D', "4" = "X")[edfinal2$Mice]
    edfinal2$exp_RNAprep <-    paste(edfinal2$experiment, edfinal2$mRNAprepdate_year, sep = '_')
    edfinal2$RNA_mRNA_prep = paste(edfinal2$RNAPrep, edfinal2$mRNAprep, sep = '_')
    edfinal2$RNA_IP_prep = paste(edfinal2$RNAPrep, edfinal2$IP, sep = '_')

    target_epi[rownames(pca$x),]$mRNAprep[as.character(target_epi[rownames(pca$x),]$mRNAprep) == '3&5'] = '3'
    cols <- MineICA::annot2Color(edfinal2)
    names(color_biocond) <- unique(edfinal$BioCond)
    cols <- c(color_biocond, cols)
  pdf(file=paste(figure_folder,"PCA_diffGenes_Input_",peak_name,"_AllTechnicalInfo.pdf",sep=""), width = 9, height = 9)
  # for (annot in c("BioCond","Exp", "Library","seq", "RNAPrep", "RNAprepdate", "RNAprepdate_monthyear", "RNAprepdate_year", "mRNAprep", "mRNAprepdate",
  #                 "mRNAprepdate_monthyear","mRNAprepdate_year","qbit", 'qbit_date', 'qbit_monthyear', 'qbit_year', "IP", "libraryprepIP", "profileinput",
  #      "concentrationinput", "profileIP", "concentrationIP", "fragmentation", "fragmentation_monthyear", "fragmentation_year","cage","experiment","exp_RNAprep")) {
  for (annot in c("BioCond", "Time", "FlowCell", "Lane", "Index", "Mice", "seq", "RNAPrep", "mRNAprep", "IP", "RNA_mRNA_prep", "RNA_IP_prep", "fractionation")) {
  pcabis <- pca
  rownames(pcabis$x) <- edfinal2[match(gsub('.*_', '', rownames(pca$x)),edfinal2$Mice),][[annot]]
  colannot <- cols[as.character(edfinal2[match(gsub('.*_', '', rownames(pca$x)),paste(edfinal2$Mice,edfinal2$Time, sep='.')),][[annot]])]
  figure <- ggplot2::autoplot(pca,label = FALSE, size=4, colour = colannot,
                              x = 1, y = 2,
                              xlabs = paste0("PC1 (",prp[1], "%)"), ylabs = paste0("PC2 (", prp[2], "%)"),
                              main = paste0("Principal Component Analysis - Axes 1 and 2, ", annot))
  print(figure)
  figure <- ggplot2::autoplot(pca,label = TRUE, size=2, colour = colannot,
                              x = 1, y = 2,
                              xlabs = paste0("PC1 (",prp[1], "%)"), ylabs = paste0("PC2 (", prp[2], "%)"),
                              main = paste0("Principal Component Analysis - Axes 1 and 2, ", annot))
  print(figure)
  figure <- ggplot2::autoplot(pcabis,label = TRUE, size=2, colour = colannot,
                              x = 1, y = 2,
                              xlabs = paste0("PC1 (",prp[1], "%)"), ylabs = paste0("PC2 (", prp[2], "%)"),
                              main = paste0("Principal Component Analysis - Axes 1 and 2, ", annot))
  print(figure)
  figure <- ggplot2::autoplot(pca,label = FALSE, size=4, colour = colannot,
                              x = 3, y = 4,
                              xlabs = paste0("PC3 (",prp[3], "%)"), ylabs = paste0("PC4 (", prp[4], "%)"),
                              main = paste0("Principal Component Analysis - Axes 3 and 4, ", annot))
  print(figure)
  figure <- ggplot2::autoplot(pca,label = TRUE, size=2, colour = colannot,
                              x = 3, y = 4,
                              xlabs = paste0("PC3 (",prp[3], "%)"), ylabs = paste0("PC4 (", prp[4], "%)"),
                              main = paste0("Principal Component Analysis - Axes 3 and 4, ", annot))
  print(figure)
  figure <- ggplot2::autoplot(pcabis,label = TRUE, size=2, colour = colannot,
                              x = 3, y = 4,
                              xlabs = paste0("PC3 (",prp[3], "%)"), ylabs = paste0("PC4 (", prp[4], "%)"),
                              main = paste0("Principal Component Analysis - Axes 3 and 4, ", annot))
  print(figure)
  }
  dev.off()
 }
  png(filename=paste(figure_folder,"PCA_NoLabel_",peak_name,type,"_PC34.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- ggplot2::autoplot(pca,label = FALSE, size=4, colour = colsidecolors,
                              x = 3, y = 4,
                              xlabs = paste0("PC3 (",prp[3], "%)"), ylabs = paste0("PC4 (", prp[4], "%)"),
                              main = "Principal Component Analysis - Axes 3 and 4")
  print(figure)
  dev.off()

  # figure <- ggplot2::autoplot(pca,label = TRUE, size=4, colour = colsidecolors,
  #                             x = 5, y = 6,
  #                             xlabs = paste0("PC4 (",prp[5], "%)"), ylabs = paste0("PC6 (", prp[6], "%)"),
  #                             main = "Principal Component Analysis - Axes 5 and 6")
  # print(figure)

  png(filename=paste(figure_folder,"PCA_",peak_name,type,".png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- ggplot2::autoplot(pca,label = TRUE, size=2, label.size = 4, colour = colsidecolors,
                              xlabs = paste0("PC1 (",prp[1], "%)"), ylabs = paste0("PC2 (", prp[2], "%)"),
                              main = "Principal Component Analysis - Axes 1 and 2")
  print(figure)
  dev.off()

  png(filename=paste(figure_folder,"MDS_",peak_name,type,"_Dim12.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- plotMDS(filter_table,
          col = colsidecolors,
          dim.plot = c(1,2))
  print(figure)
  dev.off()

  png(filename=paste(figure_folder,"MDS_NoLabel_",peak_name,type,"_Dim12.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  filter_table_nocn <- filter_table
  colnames(filter_table_nocn) <- NULL
  figure <- plotMDS(filter_table_nocn,
                    col = colsidecolors,
                    pch = 19, cex = 2,
                    dim.plot = c(1,2))
  print(figure)
  dev.off()

  png(filename=paste(figure_folder,"MDS_",peak_name,type,"_Dim34.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- plotMDS(filter_table,
                    col = colsidecolors,
                    dim.plot = c(3,4))
  print(figure)
  dev.off()


  png(filename=paste(figure_folder,"PCA_",peak_name,type,"_PC34.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  figure <- ggplot2::autoplot(pca,label = TRUE, size=2, label.size = 4, colour = colsidecolors,
                              x = 3, y = 4,
                              xlabs = paste0("PC3 (",prp[3], "%)"), ylabs = paste0("PC4 (", prp[4], "%)"),
                              main = "Principal Component Analysis - Axes 3 and 4")
  print(figure)
  dev.off()

  # png(filename=paste(project_dir,"../All_Figures/PCA/PCA_",peak_name,type,".png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  # figure <- ggplot2::autoplot(pca,label = TRUE, size=2, label.size = 4,
  #                             xlabs = paste0("PC1 (",prp[1], "%)"), colour = colsidecolors,
  #                             ylabs = paste0("PC2 (", prp[2], "%)"),
  #                             main = "Principal Component Analysis - Axes 1 and 2")
  # print(figure)
  # dev.off()
  #
  # png(filename=paste(project_dir,"../All_Figures/PCA/PCA_",peak_name,type,"_PC34.png",sep=""), width = 480*1.5, height = 480*1.5, units = "px")
  # figure <- ggplot2::autoplot(pca,label = TRUE, size=2, label.size = 4,
  #                             xlabs = paste0("PC3 (",prp[3], "%)"), colour = colsidecolors,
  #                             ylabs = paste0("PC4 (", prp[4], "%)"),
  #                             x = 3, y = 4,
  #                             main = "Principal Component Analysis - Axes 1 and 2")
  # print(figure)
  # dev.off()
}
