plot_heatmap <- function(filter_table, type, scaleRow=TRUE){
  # Heatmap
  print("Plot HeatMap")
  print(nrow(filter_table))


  if(nrow(filter_table)>1000){
     rv = apply(filter_table, 1, var, na.rm = TRUE)
     rv = rv/rowMeans(filter_table)
     filter_table = filter_table[order(rv, decreasing = TRUE),][1:1000, ]
   }
  print(nrow(filter_table))
  filter_tableor <- filter_table
  filter_table = t(scale(t(filter_table), center = TRUE, scale = scaleRow))

  dist_peak_row <- dist(filter_table, method = "euclidean")
  dist_peak_col <- dist(t(filter_table), method = "euclidean")
  ordinal_row <- hclust(dist_peak_row, method = "ward.D2")$order
  ordinal_col <- hclust(dist_peak_col, method = "ward.D2")$order
  biocondInCol <- intersect(as.character(target_epi$BioCond), names(color_biocond))
  if (all(unique(as.character(target_epi$BioCond)) %in% names(color_biocond))) {
    colsidecolors <- color_biocond[as.character(target_epi[colnames(filter_table),'BioCond'])]
    col_labels <- target_epi[colnames(filter_table),'new_name']
  } else {
    target_epi$BioCond <- as.factor(target_epi$BioCond)
    colsidecolors <- color_biocond[as.integer(target_epi[colnames(filter_table),'BioCond'])]
    col_labels <- target_epi[colnames(filter_table),'new_name']
  }

  #myPalette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")), interpolate = 'spline')(50)
  myPalette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50)
  matrixOmics <- as.matrix(filter_table)
  #matrixOmics = t(scale(t(filter_tableor), center = TRUE, scale = scaleRow))


  # plot heatmap in correct folder
  pdf(file=paste(figure_folder,"HeatMap_",peak_name,type,".pdf",sep=""),width = 7, height = 7)
  HMAbsolute <- gplots::heatmap.2(t(matrixOmics),keysize = 1, labCol = NA, scale='none',
                                  margins = c(12, 12),
                                  Colv = ordinal_row,Rowv = ordinal_col,
                                  dendrogram = 'row',labRow = col_labels, col = myPalette,
                                  key=TRUE, density.info = 'density',trace = 'none',
                                  RowSideColors = colsidecolors,
                                  hclustfun = function(x) hclust(x, method = 'ward.D2'))
  dev.off()

  # create heatmap second time for general figure folder
  # pdf(file=paste(project_dir,"../All_Figures/HeatMap/HeatMap_",peak_name,type,".pdf",sep=""),width = 7, height = 7)
  # HMAbsolute <- gplots::heatmap.2(t(matrixOmics),keysize = 1, labCol = NA,
  #                                 Colv = ordinal_row,margins = c(10, 10),Rowv = ordinal_col,
  #                                 dendrogram = 'row',labRow = col_labels, col = myPalette,
  #                                 key=TRUE, density.info = 'density',trace = 'none',
  #                                 RowSideColors = colsidecolors,
  #                                 hclustfun = function(x) hclust(x, method = 'ward.D'))
  # dev.off()
}
