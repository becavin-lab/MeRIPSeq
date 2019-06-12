load_counts_exon <- function() {
  rawDirTransExon <-
    "../../ExpressionExon/"   # path to the directory containing raw counts files

  target_Transc_I = target_RNA_I
  target_Exon_I = target_RNA_I
  target_Transc_I$files <- gsub("HTSeq_","HTSeq_Transcript_",target_RNA_I$files)
  target_Exon_I$files <- gsub("HTSeq_","HTSeq_Exon_",target_RNA_I$files)

  counts_Trans_I <-
    loadCountData(target = target_Transc_I,
                  rawDir = rawDirTransExon,
                  featuresToRemove = featuresToRemove)
  counts_Exon_I <-
    loadCountData(target = target_Exon_I,
                  rawDir = rawDirTransExon,
                  featuresToRemove = featuresToRemove)


  counts_Trans_I <-
    counts_Trans_I[rownames(counts_Trans_I) %in% transcripts$Transcript_ID, ]
  counts_Exon_I <-
    counts_Exon_I[rownames(counts_Exon_I) %in% unique(exons$Exon_ID), ]

  assign("counts_Trans_I", counts_Trans_I, pos = globalenv())
  assign("counts_Exon_I", counts_Exon_I, pos = globalenv())

}
