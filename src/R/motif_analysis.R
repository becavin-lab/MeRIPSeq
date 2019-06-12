  #'   #' Title
  #'   #' Need to have motif search and then motif occurence analysiswith fimo
  #'   #' @return
  #'   #' @export
  #'   #'
  #'   #' @examples
args <- commandArgs(trailingOnly = TRUE)
exp_design_name <- args[1]
bed_name <- args[2]
peak_name <- paste(exp_design_name,"_",bed_name,sep="")

  #exp_design_name = "Cecum_all_withoutabxGF_MaxMaxValues_3"
  #exp_design_name = "Liver_all_MaxMaxValues_3"
  #exp_design_name = "Liver_all_T13_MaxMaxValues_3"
  motif_folder = "/Volumes/m6aAkker/PeakDiffExpression/Motif/"
  #type_motif = "_Diff"
  type_motif = ""
  motif_file = paste("Motif_", exp_design_name, type_motif, sep="")
  figure_motif_folder = paste(motif_folder,motif_file,"/", sep="")
  if (!dir.exists(figure_motif_folder)){
    dir.create(figure_motif_folder)
  }

  # create fasta
  # load(paste("/Volumes/m6aAkker-1/PeakDiffExpression/",exp_design_name,"/","Cecum_all_MaxMaxValues_3",".RData",sep=""))
  #   type_overlaps = c("3UTR","5UTR","CDS","Intron")
  # folderFasta <- paste("/Volumes/m6aAkker-1/PeakDiffExpression/",exp_design_name,"/Motif/fasta/",sep="")
  # for(type_overlap in type_overlaps){
  #   matrix_type = peaks[grep(type_overlap,peaks$Type_overlap),]
  #   seq_fasta <- fasta[names(fasta) %in% matrix_type$ID]
  #   print(paste(type_overlap,"Extract:",length(seq_fasta)))
  #   norm.filename = paste(folderFasta,exp_design_name,"_",type_overlap,".fasta",sep = "")
  #   seqinr::write.fasta(sequences=seq_fasta,names=names(seq_fasta),file.out=norm.filename)
  # }


motif_analysis <- function(motif_folder, exp_design_name, motif_file){
      # Graphic parameters
      width_hm = 10
      height_hm = 14

      # read occurence file
      count_filename = paste(motif_folder,motif_file,"_fasta_count.txt",sep="")
      motif_count <- read.table(count_filename, header=TRUE, check.names = F, sep="\t", row.names = 1)
      score_filename = paste(motif_folder,motif_file,"_fasta_score.txt",sep="")
      motif_score <- read.table(score_filename, header=TRUE, check.names = F, sep="\t", row.names = 1)
      print(paste(motif_folder,motif_file,".txt",sep=""))
      motif_array <- read.delim(paste(motif_folder,motif_file,".txt",sep=""), check.names = F, stringsAsFactors=F)

      # all filenames
      log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count.pdf",sep="")
      #log2count_table_filename = paste(motif_folder,motif_file,"_log2count.txt",sep="")
      log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_list.txt",sep="")
      fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore.pdf",sep="")
      #fimo_table_filename = paste(motif_folder,motif_file,"_fimoscore.txt",sep="")
      fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_list.txt",sep="")

      # pos_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Pos.pdf",sep="")
      # pos_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Pos_list.txt",sep="")
      # pos_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Pos.pdf",sep="")
      # pos_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Pos_list.txt",sep="")
      pos_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Pos_NoIntron.pdf",sep="")
      pos_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Pos_NoIntron_list.txt",sep="")
      pos_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Pos_NoIntron.pdf",sep="")
      pos_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Pos_NoIntron_list.txt",sep="")

      biocond_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Biocond.pdf",sep="")
      biocond_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Biocond_list.txt",sep="")
      biocond_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Biocond.pdf",sep="")
      biocond_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Biocond_list.txt",sep="")

      comp_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp.pdf",sep="")
      comp_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list.txt",sep="")
      comp_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp.pdf",sep="")
      comp_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list.txt",sep="")

      # Get all Meth list : contains(project_name)
      column_filter = grep(exp_design_name,colnames(motif_count))
      filter_motif_count = motif_count[,column_filter]
      filter_motif_score = motif_score[,column_filter]
      filter_motif_count[,paste(exp_design_name,"_CDS",sep="")] = NULL
      filter_motif_count[,paste(exp_design_name,"_3UTR",sep="")] = NULL
      filter_motif_count[,paste(exp_design_name,"_5UTR",sep="")] = NULL
      filter_motif_count[,paste(exp_design_name,"_Intron",sep="")] = NULL

      # log2count figure, table and list
      pdf(file=log2count_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
      #write.table(filter_motif_count[hm$tree_row$order,hm$tree_col$order],log2count_table_filename,sep="\t",col.names = T,row.names = T,quote = F)
      write.table(motif_array[hm$tree_row$order,],log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()
      # fimscore fiogure, table, and list
      pdf(file=fimo_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
      #write.table(filter_motif_count[hm$tree_row$order,hm$tree_col$order],fimo_table_filename,sep="\t",col.names = T,row.names = T,quote = F)
      write.table(motif_array[hm$tree_row$order,],fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()

      # Get all positions
      # Get all Meth list : contains(project_name)
      #col_filter = c('_Meth_CDS','_Meth_3UTR')
      #col_filter = c('_Meth_CDS','_Meth_Intron','_Meth_3UTR')
      col_filter = c('_CDS','_3UTR','_5UTR')
      col_filter = paste(exp_design_name,col_filter,sep='')
      filter_motif_count = motif_count[,col_filter]
      filter_motif_score = motif_score[,col_filter]


      # log2count figure, table and list
      pdf(file=pos_log2count_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
      write.table(motif_array[hm$tree_row$order,],pos_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()
      # fimscore fiogure, table, and list
      pdf(file=pos_fimo_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
      write.table(motif_array[hm$tree_row$order,],pos_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()

      # Get all Biocond :
      comp_list = grep("-",colnames(motif_count))
      filter_col <- colnames(motif_count)[!(1:length(colnames(motif_count)) %in% comp_list)]
      meth_list = grep("_Meth",filter_col)
      filter_col <- filter_col[!(1:length(filter_col) %in% meth_list)]
      exp_design = grep(exp_design_name,filter_col)
      filter_col <- filter_col[!(1:length(filter_col) %in% exp_design)]
      filter_motif_count = motif_count[,filter_col]
      filter_motif_score = motif_score[,filter_col]
      # Save log2count for BioCond
      pdf(file=biocond_log2count_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
      write.table(motif_array[hm$tree_row$order,],biocond_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()
      # Save fimoscore for BioCond
      pdf(file=biocond_fimo_figure_filename,width = width_hm, height = height_hm)
      hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
      write.table(motif_array[hm$tree_row$order,],biocond_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
      dev.off()

      # Get all Comparisons :
      comp_list = grep("-",colnames(motif_count))
      filter_col <- colnames(motif_count)[(1:length(colnames(motif_count)) %in% comp_list)]
      meth_list = grep("_Meth",filter_col)
      filter_col <- filter_col[!(1:length(filter_col) %in% meth_list)]

      # If both ZT13 and ZT3 are included
      if(length(grep("Liver_all_Max",exp_design_name))==1){
        comp_ZT3_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_T3.pdf",sep="")
        comp_ZT3_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list_T3.txt",sep="")
        comp_ZT3_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_T3.pdf",sep="")
        comp_ZT3_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list_T3.txt",sep="")
        comp_ZT13_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_T13.pdf",sep="")
        comp_ZT13_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list_T13.txt",sep="")
        comp_ZT13_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_T13.pdf",sep="")
        comp_ZT13_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list_T13.txt",sep="")
        comp_ZT3ZT13_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_T3T13.pdf",sep="")
        comp_ZT3ZT13_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list_T3T13.txt",sep="")
        comp_ZT3ZT13_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_T3T13.pdf",sep="")
        comp_ZT3ZT13_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list_T3T13.txt",sep="")

        filter_col_ZT3 = filter_col[grep("T3",filter_col)]
        filter_col_ZT13 = filter_col[grep("T13",filter_col)]
        filter_col_ZT3ZT13 = filter_col_ZT13[filter_col_ZT13 %in% filter_col_ZT3]
        filter_col_ZT3 = filter_col_ZT3[!(filter_col_ZT3 %in% filter_col_ZT13)]
        filter_col_ZT13 = filter_col_ZT13[!(filter_col_ZT13 %in% filter_col_ZT3ZT13)]
        filter_col_ZT13 [! filter_col_ZT13 %in% "T13.Lp-T13.Ec"]

        # ZT3 only comparisons
        filter_col = filter_col_ZT3
        comp_log2count_figure_filename = comp_ZT3_log2count_figure_filename
        comp_log2count_list_filename = comp_ZT3_log2count_list_filename
        comp_fimo_figure_filename = comp_ZT3_fimo_figure_filename
        comp_fimo_list_filename = comp_ZT3_fimo_list_filename
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        # Save log2count for Comparisons
        pdf(file=comp_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()


        # ZT13 only comparisons
        filter_col = filter_col_ZT13
        comp_log2count_figure_filename = comp_ZT13_log2count_figure_filename
        comp_log2count_list_filename = comp_ZT13_log2count_list_filename
        comp_fimo_figure_filename = comp_ZT13_fimo_figure_filename
        comp_fimo_list_filename = comp_ZT13_fimo_list_filename
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        # Save log2count for Comparisons
        pdf(file=comp_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()

        # ZT3 vs ZT13 only comparisons
        filter_col = filter_col_ZT3ZT13
        comp_log2count_figure_filename = comp_ZT3ZT13_log2count_figure_filename
        comp_log2count_list_filename = comp_ZT3ZT13_log2count_list_filename
        comp_fimo_figure_filename = comp_ZT3ZT13_fimo_figure_filename
        comp_fimo_list_filename = comp_ZT3ZT13_fimo_list_filename
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        # Save log2count for Comparisons
        pdf(file=comp_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()

      }else if(length(grep("Liver_all_T13_Max",exp_design_name))==1){
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        # Save log2count for Comparisons
        pdf(file=comp_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_array[hm$tree_row$order,],comp_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()

      }else{

        comp_CONV_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_CONV.pdf",sep="")
        comp_CONV_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list_CONV.txt",sep="")
        comp_CONV_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_CONV.pdf",sep="")
        comp_CONV_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list_CONV.txt",sep="")
        comp_GF_log2count_figure_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_GF.pdf",sep="")
        comp_GF_log2count_list_filename = paste(figure_motif_folder,motif_file,"_log2count_Comp_list_GF.txt",sep="")
        comp_GF_fimo_figure_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_GF.pdf",sep="")
        comp_GF_fimo_list_filename = paste(figure_motif_folder,motif_file,"_fimoscore_Comp_list_GF.txt",sep="")


        # Filter first the CONV datasets
        filter_col = c('GF-CONV', 'CONV-abx', 'GF-ex_GF', 'ex_GF-abx', 'CONV-Am', 'Lp-CONV')
        motif_filter = motif_array[c(grep("manually",motif_array$Data),(grep("CONV",motif_array$Data))),]
        rownames(motif_filter) = motif_filter$Motif
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        if(type_motif == "_Diff"){
          filter_motif_count = motif_count[motif_filter$Motif,filter_col]
          filter_motif_score = motif_score[motif_filter$Motif,filter_col]
        }
        # Save log2count for Comparisons
        pdf(file=comp_CONV_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_filter[hm$tree_row$order,],comp_CONV_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_CONV_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_filter[hm$tree_row$order,],comp_CONV_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()

        # Then filter GF datasets
        filter_col = c('GF-CONV', 'Lp-GF')
        motif_filter = motif_array[c(grep("manually",motif_array$Data),(grep("GF",motif_array$Data))),]
        rownames(motif_filter) = motif_filter$Motif
        filter_motif_count = motif_count[,filter_col]
        filter_motif_score = motif_score[,filter_col]
        if(type_motif == '_Diff'){
          filter_motif_count = motif_count[motif_filter$Motif,filter_col]
          filter_motif_score = motif_score[motif_filter$Motif,filter_col]
        }
        # Save log2count for Comparisons
        pdf(file=comp_GF_log2count_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(log2(filter_motif_count), scale = "row")
        write.table(motif_filter[hm$tree_row$order,],comp_GF_log2count_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()
        # Save fimoscore for Comparisons
        pdf(file=comp_GF_fimo_figure_filename,width = width_hm, height = height_hm)
        hm <- pheatmap::pheatmap(filter_motif_score, scale = "row")
        write.table(motif_filter[hm$tree_row$order,],comp_GF_fimo_list_filename,sep="\t",col.names = T,row.names = F,quote = F)
        dev.off()

      }


      # to show coorelation for the two group of motif
      #g <- factor(cutree(hm$tree_row, 2))
      #ggpairs(data.frame(filter_motif_score, g), columns = 1:12,
      #      ggplot2::aes(colour=g, alpha=0.5))

}

motif_analysis(motif_folder, exp_design_name, motif_file)

