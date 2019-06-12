plot_clustering <- function(){
  
  # clustering without seq batch
  print("Plot Clustering")
  # display cluering for Epi, Input
  list_type <- list("Epi", "Input")
  norm_methyl <- data.frame(norm_voom_epi)
  I_gene <- data.frame(norm_voom_RNA_I)
  list_table <- list(norm_methyl,I_gene)
  #i=1
  for(i in 1:2){
    type = list_type[[i]]
    filter_table <- list_table[[i]]
    colsidecolors <- color_biocond[as.integer(target_RNA_I[colnames(filter_table),'Type'])]
    filter_table = t(scale(t(filter_table), center = TRUE))
    dist_peak_col <- dist(t(filter_table), method = "euclidean")
    nodePar <- list(col = colsidecolors)
    hc <- hclust(dist_peak_col, method = "ward.D2")
    hcd <- as.dendrogram(hc)
    png(filename = paste(figure_folder,"HCluster_",project_name,"_Genes_",type,".png",sep=""), width = 3000, height = 1000,res = 150)
    dend2 <- color_labels(hcd, labels = colnames(filter_table) , col = colsidecolors[hc$order]) 
    labels_cex(dend2)<- 0.6
    par(oma=c(1,1,1,1))
    plot(dend2,ylab = "Height", main = paste(project_name,"Clustering",type))
    dev.off()
    
    png(filename = paste(figure_folder,'../../HCluster_',project_name,"_Genes_",type,".png",sep=""), width = 3000, height = 1000,res = 150)
    dend2 <- color_labels(hcd, labels = colnames(filter_table) , col = colsidecolors[hc$order]) 
    labels_cex(dend2)<- 0.6
    par(oma=c(1,1,1,1))
    plot(dend2,ylab = "Height", main = paste(project_name,"Clustering",type))
    dev.off()
    
  }
  
}