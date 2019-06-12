load_counts <- function() {
  rawDirRNA <-
    "../../Expression/"   # path to the directory containing raw counts files
  rawDirEpi <- "../../Expression/"

  counts_epi <- loadCountData(target = target_epi,rawDir = rawDirEpi,featuresToRemove = featuresToRemove)
  counts_RNA_I <-
    loadCountData(target = target_RNA_I,
                  rawDir = rawDirRNA,
                  featuresToRemove = featuresToRemove)
  counts_RNA_IP <-
     loadCountData(target = target_RNA_IP,
                   rawDir = rawDirRNA,
                   featuresToRemove = featuresToRemove)

  counts_epi <- counts_epi[rownames(counts_epi) %in% peaks$ID, ]
  counts_RNA_I <-
    counts_RNA_I[rownames(counts_RNA_I) %in% genes$Gene_ID, ]
  counts_RNA_IP <-
     counts_RNA_IP[rownames(counts_RNA_IP) %in% genes$Gene_ID, ]

  # filter counts
  # keep genes with at least one cpm in at least 3 Input samples
  counts_RNA_I <- counts_RNA_I[which(rowSums(edgeR::cpm(counts_RNA_I)>=1) > 2),]
  counts_RNA_IP <- counts_RNA_IP[rownames(counts_RNA_I),]
  # keep peaks with at least 5 cpm in at least 3 IP and 3 Input samples
  counts_epi <- counts_epi[which(rowSums(edgeR::cpm(counts_epi[, grep('_IP',colnames(counts_epi))])>=5) > 2 &
         rowSums(edgeR::cpm(counts_epi[, grep('_.nput',colnames(counts_epi))])>=5) > 2),]

  # # for our merged dataset only:
  # #  selection of peaks expressed in at least half samples (or 2 samples) per condition and per dataset (2019 or before)
  # counts_epi <- counts_epi[(unique(unlist(lapply(split(colnames(counts_epi),target_epi$BioCond), function(cn) {
  #    ct = counts_epi[,cn];
  #    ctepi <- ct[, grep('_IP',colnames(ct))]
  #    ctinput <- ct[, grep('_.nput',colnames(ct))]
  #
  #    if (length(grep('2019',colnames(ct), invert=T)) > 0 & length(grep('2019',colnames(ct), invert=F)) > 0) {
  #      intersect(which(rowSums(edgeR::cpm(ctepi[, grep('2019',colnames(ctepi))]) >=5) >= 2 &
  #                        rowSums(edgeR::cpm(ctinput[, grep('2019',colnames(ctinput))]) >=5) >= 2),
  #                which(rowSums(edgeR::cpm(ctepi[, grep('2019',colnames(ctepi), invert = T)]) >=5) >= 2 &
  #                        rowSums(edgeR::cpm(ctinput[, grep('2019',colnames(ctinput), invert = T)]) >=5) >= 2 ))
  #    } else {
  #      which(rowSums(edgeR::cpm(ctepi) >=5) >= min(4, ncol(ctepi)) &
  #              rowSums(edgeR::cpm(ctinput) >=5) >= min(4, ncol(ctepi)))
  #    }
  #    }))))
  #  ,]

  assign("counts_epi", counts_epi, pos = globalenv())
  assign("counts_RNA_I", counts_RNA_I, pos = globalenv())
  assign("counts_RNA_IP", counts_RNA_IP, pos = globalenv())

}
