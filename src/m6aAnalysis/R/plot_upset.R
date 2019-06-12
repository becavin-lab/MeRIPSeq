library(UpSetR)

# this function plots upsetR plots showing:
# the number of differentially methylated peaks and differentially expressed  genes, per comparison
# the number of differentially methylated peaks that are not associated withdifferentially expressed gene, per comparison

plot_upset <- function(filename) {

  # folder that will contain lists of peaks for the most frequent sets
  intersect_folder <- file.path(project_dir, 'intersections')
  dir.create(intersect_folder);

  # function to extract the 20 most frequent sets of upsetR plots
  writeTopSets <- function(isdiff_epi, annot, filename = '',  colId = 'ID', nbSets=40) {
    occ <- apply(isdiff_epi,1,paste0,collapse='')
    topint <- names(sort((table(occ)), decreasing = T)[2:nbSets]) # "0000000000000000000000000" is always first
    nbSets <- nbSets-(sum(is.na(topint)))
    topint <- topint[1:(nbSets-1)]
    freqs <- sort((table(occ)), decreasing = T)[2:nbSets]
    topint_peaks <- split(rownames(isdiff_epi), occ)[topint]
    topint_comp <- apply(isdiff_epi[match(topint, occ),], 1,
                         function(y) paste(colnames(isdiff_epi)[which(y == 1)], collapse=','))
    topint_nbcomp <- apply(isdiff_epi[match(topint, occ),], 1,
                         function(y) length(colnames(isdiff_epi)[which(y == 1)]))

    occ2peaks <- do.call(rbind,mapply(x=topint_peaks, y=topint_comp, freq=freqs, nb=topint_nbcomp,
                                      function(x,y,freq,nb) cbind(comparisons = y,
                                                               nb_comparisons = nb,
                                                               freq = freq,
                                                               annot[match(x, annot[[colId]]),]), SIMPLIFY = F))
    occ2peaks <- occ2peaks[order(occ2peaks$nb_comparisons, decreasing=T),]
    write.table(occ2peaks, file = filename, sep = '\t', quote = F, row.names = F)
  }


  pdf( file = file.path(figure_folder, filename), 10, 8)
  # for diff peaks, lfc 1
  isdiff_epi <- do.call(cbind,lapply(contr_Epi_names, function(i) {
    tt <- (limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf));
    as.numeric(abs(tt$logFC)>1 & tt$adj.P.Val < 0.05 )
  }))
  colnames(isdiff_epi) <- contr_Epi_names
  rownames(isdiff_epi) <- rownames(fit_bayes_epi)
  upset(as.data.frame(isdiff_epi), nsets = length(contr_Epi_names), point.size = 1,  order.by = c('freq'), mb.ratio = c(.5,.5), group.by = 'degree')
  #, decreasing = c(TRUE, FALSE))#, keep.order = TRUE, sets = contr_Epi_names)
  grid::grid.text('Differentially methylated peaks \n |logFC| > 1 & FDR < 5%',x = 0.7, y=.9)

  writeTopSets(isdiff_epi, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_lfc1.txt'),
               colId = 'ID')

  # for diff peaks, lfc 2
  isdiff_epi <- do.call(cbind,lapply(contr_Epi_names, function(i) {
    tt <- (limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf));
    as.numeric(abs(tt$logFC)>2 & tt$adj.P.Val < 0.05 )
  }))
  colnames(isdiff_epi) <- contr_Epi_names
  rownames(isdiff_epi) <- rownames(fit_bayes_epi)
  upset(as.data.frame(isdiff_epi), nsets = length(contr_Epi_names), point.size = 1,  order.by = c('freq'), mb.ratio = c(.5,.5), group.by = 'degree')#, decreasing = c(TRUE, FALSE))#, keep.order = TRUE, sets = contr_Epi_names)
  grid::grid.text('Differentially methylated peaks \n |logFC| > 2 & FDR < 5%',x = 0.7, y=.9)

  writeTopSets(isdiff_epi, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_lfc2.txt'),
               colId = 'ID')



  # for diff genes, lfc 1
  isdiff_epi <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
    tt <- (limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf));
    as.numeric(abs(tt$logFC)>1 & tt$adj.P.Val < 0.05 )
  }))
  colnames(isdiff_epi) <- gsub('_Input', '', contr_RNA_I_names)
  rownames(isdiff_epi) <- rownames(fit_bayes_RNA_I)
  upset(as.data.frame(isdiff_epi), nsets = length(contr_RNA_I_names), point.size = 1,  order.by = c('freq'), mb.ratio = c(.5,.5), group.by = 'degree')
  #, decreasing = c(TRUE, FALSE))#, keep.order = TRUE, sets = contr_Epi_names)
  grid::grid.text('Differentially expressed genes \n |logFC| > 1 & FDR < 5%',x = 0.7, y=.9)

  writeTopSets(isdiff_epi, annot = genes,
               filename = file.path(intersect_folder,'intersects_diff_genes_lfc1.txt'),
               colId = 'Gene_ID')

  # build a matrix which says for each DE peak whether its associated gene (using peaks_to_gene) is DE too
  isdiff_both <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
    ttp <- limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- ttg[peaks_to_gene[rownames(ttp),]$Gene_ID,]
    ttg[grep('NA',(rownames(ttg))),] <- 0
    as.numeric(abs(ttp$logFC)>1 & ttp$adj.P.Val < 0.05 & abs(ttg$logFC)>1 & ttg$adj.P.Val < 0.05)
  }))
  colnames(isdiff_both) <- gsub('_Input', '', contr_RNA_I_names)
  rownames(isdiff_both) <- rownames(fit_bayes_epi)

  writeTopSets(isdiff_both, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_and_genes_lfc1.txt'),
               colId = 'ID')

  isdiff_both_upup <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
    ttp <- limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- ttg[peaks_to_gene[rownames(ttp),]$Gene_ID,]
    ttg[grep('NA',(rownames(ttg))),] <- 0
    as.numeric((ttp$logFC)>1 & ttp$adj.P.Val < 0.05 & (ttg$logFC)>1 & ttg$adj.P.Val < 0.05)
  }))
  colnames(isdiff_both_upup) <- paste0(gsub('_Input', '', contr_RNA_I_names), '_', 'upup')
  rownames(isdiff_both_upup) <- rownames(fit_bayes_epi)

  writeTopSets(isdiff_both_upup, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_and_genes_lfc1_UpUp.txt'),
               colId = 'ID')

  isdiff_both_downdown <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
    ttp <- limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- ttg[peaks_to_gene[rownames(ttp),]$Gene_ID,]
    ttg[grep('NA',(rownames(ttg))),] <- 0
    as.numeric((ttp$logFC) < -1 & ttp$adj.P.Val < 0.05 & (ttg$logFC) < -1 & ttg$adj.P.Val < 0.05)
  }))
  colnames(isdiff_both_downdown) <- paste0(gsub('_Input', '', contr_RNA_I_names),'_','downdown')
  rownames(isdiff_both_downdown) <- rownames(fit_bayes_epi)

  writeTopSets(isdiff_both_downdown, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_and_genes_lfc1_DownDown.txt'),
               colId = 'ID')

  nb <- data.frame(sets=colnames(isdiff_both),
              "MethUp_GeneUp"=colSums(isdiff_both_upup, na.rm = T),
              "MethDown_GeneDown"=colSums(isdiff_both_downdown, na.rm = T))
  rownames(nb) <- nb$Name

  countNbupup <- function(row) {
   data <-  sum(row[grep('upup',names(row))] == 4) > 0
  }
  upset(as.data.frame(cbind(isdiff_both)),
                            #apply(isdiff_both_upup, 2, function(x) c('TRUE'=4,'FALSE'=5)[as.character(x)]),
                            #apply(isdiff_both_downdown, 2, function(x) c('TRUE'=4,'FALSE'=5)[as.character(x)]))),
              nsets = length(contr_Epi_names), #sets = contr_Epi_names,
              point.size = 1, order.by = c('freq'), mb.ratio = c(.5,.5),
              group.by = 'degree')#,
              # set.metadata = list(data = (nb),
              #                     # plots = list(list(type = "text", column = "MethUp_GeneUp", assign = 5),
              #                     #              list(type = "text", column = "MethDown_GeneDown", assign = 5))
              #                     plots = list(list(type = "hist", column = "MethUp_GeneUp", assign = 20),
              #                                  list(type = "hist", column = "MethDown_GeneDown", assign = 20)))
              #

        #queries = list(list(query = countNbupup, color = 'blue', active = T)) )
        #, decreasing = c(TRUE, FALSE))#, keep.order = TRUE, sets = contr_Epi_names)
  grid::grid.text('Differentially methylated peak associated \n with a differentially expressed gene \n |logFC| > 1 & FDR < 5%',x = 0.7, y=.9)


  colnames(isdiff_both_upup) <- gsub('_upup', '', colnames(isdiff_both_upup))
  upset(as.data.frame(apply(isdiff_both_upup,2,as.numeric)),
                            nsets = length(contr_Epi_names), #sets = contr_Epi_names,
                            point.size = 1, order.by = c('freq'), mb.ratio = c(.5,.5),
                            group.by = 'degree')
  grid::grid.text('Differentially methylated peak associated \n with a differentially expressed gene \n logFC > 1 & FDR < 5%',x = 0.7, y=.9)


  colnames(isdiff_both_downdown) <- gsub('_downdown', '', colnames(isdiff_both_downdown))
  upset(as.data.frame(apply(isdiff_both_downdown,2,as.numeric)),
                            nsets = length(contr_Epi_names), #sets = contr_Epi_names,
                            point.size = 1, order.by = c('freq'), mb.ratio = c(.5,.5),
                            group.by = 'degree')
  grid::grid.text('Differentially methylated peak associated \n with a differentially expressed gene \n logFC < -1 & FDR < 5%',x = 0.7, y=.9)

  # for diff peaks summarized at the gene level, lfc 1
  # e.g. each peak is transformed to its gene id, and repeats of the same genes are filtered out in the data matrix
  # build a matrix which says for each DE peak whether its associated gene (using peaks_to_gene) is DE too
  # isdiff_both <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
  #   ttp <- limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
  #   ttg <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
  #   ttg <- ttg[peaks_to_gene[rownames(ttp),]$Gene_ID,]
  #   ttg[grep('NA',(rownames(ttg))),] <- 0
  #   as.numeric(abs(ttp$logFC)>1 & ttp$adj.P.Val < 0.05 & abs(ttg$logFC)>1 & ttg$adj.P.Val < 0.05)
  # }))
  # colnames(isdiff_both) <- gsub('_Input', '', contr_RNA_I_names)
  # rownames(isdiff_both) <- rownames(fit_bayes_epi)
  #
  isdiff_notboth <- do.call(cbind,lapply(gsub('_Input', '', contr_RNA_I_names), function(i) {
    ttp <- limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf);
    ttg <- ttg[peaks_to_gene[rownames(ttp),]$Gene_ID,]
    ttg[grep('NA',(rownames(ttg))),] <- 0
    as.numeric(abs(ttp$logFC)>1 & ttp$adj.P.Val < 0.05 & !(abs(ttg$logFC)>1 & ttg$adj.P.Val < 0.05))
  }))
  colnames(isdiff_notboth) <- gsub('_Input', '', contr_RNA_I_names)
  rownames(isdiff_notboth) <- rownames(fit_bayes_epi)

  writeTopSets(isdiff_notboth, annot = peaks,
               filename = file.path(intersect_folder,'intersects_diff_peaks_and_NOTgenes_lfc1.txt'),
               colId = 'ID')

  upset(as.data.frame(apply(isdiff_notboth,2,as.numeric)),
        nsets = length(contr_Epi_names), #sets = contr_Epi_names,
        point.size = 1, order.by = c('freq'), mb.ratio = c(.5,.5),
        group.by = 'degree')
  grid::grid.text('Differentially methylated peak NOT associated \n with a differentially expressed gene \n |logFC| > 1 & FDR < 5%',x = 0.7, y=.9)


  dev.off()

  # diff peaks + diff genes
  #upset(as.data.frame(isdiff_epi[,grep("abx_GF|abxGF", colnames(isdiff_epi), invert = T)]), nsets = length(contr_RNA_I_names), point.size = 1,  order.by = c('freq'), mb.ratio = c(.5,.5), group.by = 'degree')

  #sort((table(apply(isdiff_epi,1,paste0,collapse=''))))
  #dir.create(file.path(project_dir, 'intersections')); write.table(peaks[names(which(apply(isdiff_epi,1,paste0,collapse='') == '101010101111111101010001110010111101')),], file = file.path(project_dir, 'intersections', '8peaks_commonTo23comp.txt'), sep = '\t', quote = F)
}
