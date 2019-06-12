library(GenomicFeatures)
library(data.table)


plotVenns <- function(exp1='PeakDiffExpression/Cecum_all_withoutabxGF/Cecum_all_withoutabxGF.RData', 
                      exp2='PeakDiffExpression/Liver_all_T13/Liver_all_T13.RData',
                      bed_ref="Genome/pc_ep_mouse_mm10.bed",
                      filename="venn_LiverAll_CecumAll_ref.pdf",
                      select_counts=FALSE) {


  load(exp1)
  
  # load peak annotations of exp1
  peaksgr <- copy(peaks)
  setnames(peaksgr, c("chromo_peak",   "begin_peak",    "end_peak", 'chr'), c('chr', 'start', 'end','chr_gene'))
  peaksgr_cecum <- makeGRangesFromDataFrame(peaksgr, keep.extra.columns = TRUE, seqnames.field="chr")
  if (select_counts) peaksgr_cecum <- peaksgr_cecum[rownames(counts_epi)]
  names(peaksgr_cecum) <- paste0('Cecum_',names(peaksgr_cecum))
  
  # load MERIP reference sites
  m6aref <- fread(bed_ref)#sb_m6a_mouse_mm10.bed")
  setnames(m6aref, c('chr', 'start', 'end', 'type', 'id', 'strand'))
  m6arefgr <- makeGRangesFromDataFrame(as.data.frame(m6aref), keep.extra.columns = TRUE, start.field = "start", end.field = 'end', seqnames.field="chr")
  m6arefgr <- reduce(m6arefgr)
  names(m6arefgr) <- paste0('ref_',1:length(m6arefgr))
  
  # compute overlap between exp1 peaks and ref
  ov <- findOverlaps(query = peaksgr_cecum, subject = m6arefgr, minoverlap = 1)
  ov <- as.data.table(ov)
  ov[, peak := names(peaksgr_cecum)[queryHits]]
  ov[, ref := names(m6arefgr)[subjectHits]]
  ov[, ovWidth := width(pintersect(peaksgr_cecum[peak],m6arefgr[ref])), by = .I]
  # assign name of the ref peak that has the biggest overlap to the exp1 peak
  peaks_to_ref_cecum <- ov[, ref[which.max(ovWidth)], by = peak]
  names(peaksgr_cecum)[match(peaks_to_ref_cecum$peak,names(peaksgr_cecum))] <- peaks_to_ref_cecum$V1
  #ov[, refm6a := names(m6arefgr)[subjectHits]]
  
  # load peak annotations of exp2
  load(exp2)
  peaksgr <- copy(peaks)
  setnames(peaksgr, c("chromo_peak",   "begin_peak",  "end_peak", 'chr'), c('chr', 'start', 'end','chr_gene'))
  peaksgr_liver <- makeGRangesFromDataFrame(peaksgr, keep.extra.columns = TRUE, seqnames.field="chr")
  if (select_counts) peaksgr_liver <- peaksgr_liver[rownames(counts_epi)]
  names(peaksgr_liver) <- paste0('Liver_',names(peaksgr_liver))
  
  # compute overlap between exp2 peaks and ref
  ov_liverref <- findOverlaps(query = peaksgr_liver, subject = m6arefgr, minoverlap = 1)
  ov_liverref <- as.data.table(ov_liverref)
  ov_liverref[, peak :=names(peaksgr_liver)[queryHits]]
  ov_liverref[, ref := names(m6arefgr)[subjectHits]]
  ov_liverref[, ovWidth := width(pintersect(peaksgr_liver[peak],m6arefgr[ref])), by = .I]
  # assign name of the ref peak that has the biggest overlap to the exp2 peak
  peaks_to_ref_liver <- ov_liverref[, ref[which.max(ovWidth)], by = peak]
  names(peaksgr_liver)[match(peaks_to_ref_liver$peak,names(peaksgr_liver))] <- peaks_to_ref_liver$V1
  
  # compute overlap between exp1 and exp2
  ov_cecum_liver <- as.data.table(findOverlaps(query = peaksgr_liver, subject = peaksgr_cecum, minoverlap = 1))
  ov_cecum_liver[, peakCecum := names(peaksgr_cecum)[subjectHits]]
  ov_cecum_liver[, peakLiver := names(peaksgr_liver)[queryHits]]
  ov_cecum_liver[, ovWidth := width(pintersect(peaksgr_liver[peakLiver],peaksgr_cecum[peakCecum])), by = .I]
  # assign same name to the peak overlapping between exp1 and exp2, using the max overlap only
  peaks_to_cecum_liver <- ov_cecum_liver[, peakCecum[which.max(ovWidth)], by = peakLiver]
  names(peaksgr_liver)[match(peaks_to_cecum_liver$peakLiver,names(peaksgr_liver))] <- peaks_to_cecum_liver$V1
  
  if (!missing(filename) | !is.null(filename)) pdf(filename, 7,7)
  venn(data = list(Liver=names(peaksgr_liver), Cecum=names(peaksgr_cecum), "reference peaks"=names(m6arefgr)))
  if (!missing(filename) | !is.null(filename)) dev.off()
  
  #length(unique(ov$peak))/length(peaksgr_cecum)
  #length(unique(ov$peak));length(peaksgr_cecum)
  #length(unique(ov_liverref$peak))/length(peaksgr_liver)
  #length(unique(ov_cecum_liver$peakCecum))/length(peaksgr_cecum)
  #length(unique(ov_cecum_liver$peakLiver))/length(peaksgr_liver)
  #summary(as.numeric(table(ov_liverref$ref)))
}

plotVenns()

