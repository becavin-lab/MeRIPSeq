load_exon <- function(){

  print("Load Transcripts")
  transcripts <- read.delim(file=paste("../../Genome/",annotation_name,".annotation.transcript.txt",sep=""), sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(transcripts),"transcripts loaded"))
  transcripts$Gene_ID = transcripts$gene_id
  transcripts$gene_id <- NULL

  print("Load Peaks for Transcripts")
  peaks_transcripts <- read.delim(file=paste("../../PeakDetection/Peaks/",peak_name,"_Peaks_Ref_Transcript_Relat.txt",sep = ""), sep = "\t", row.names = 1, stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(peaks_transcripts),"peaks loaded"))

  print("Load Peaks for Exons")
  peaks_exons <- read.delim(file=paste("../../PeakDetection/Peaks/",peak_name,"_Peaks_Ref_Transcript_Relat_Exon.txt",sep = ""), sep = "\t", row.names = 1, stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(peaks_exons),"peaks loaded"))

  print("Number of Transcripts per Gene")
  trans_stats <- read.delim(file=paste("../../Genome/",annotation_name,".annotation.transcript.stats.txt",sep=""), sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  print("Gene to transcript stats")
  summary(trans_stats$Nb_Transcript)
  value = length(trans_stats$Nb_Transcript[trans_stats$Nb_Transcript==1])
  print(paste("Number of Genes with 1 Transcript:",value," = ",(value/nrow(genes)),"%"))
  value = length(trans_stats$Nb_Transcript[trans_stats$Nb_Transcript==2])
  print(paste("Number of Genes with 2 Transcripts:",value," = ",(value/nrow(genes)),"%"))
  value = length(trans_stats$Nb_Transcript[trans_stats$Nb_Transcript>2])
  print(paste("Number of Genes with more than 2 Transcripts:",value," = ",(value/nrow(genes)),"%"))
  png(filename = paste(figure_folder,"TranscriptPerGenes.png",sep=""), width = 1000, height = 1000,res = 150)
  hist(trans_stats$Nb_Transcript, main = "Number of transcripts per gene", xlab = "Nb of transcripts per gene", breaks = 100)
  dev.off()

  print("Load Exons")
  exons <- read.delim(file=paste("../../Genome/",annotation_name,".annotation.exon.txt",sep=""), sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(exons),"exons by transcripts loaded"))
  transcripts$Gene_ID = transcripts$gene_id
  transcripts$gene_id <- NULL

  print("Load exon to gene and transcripts")
  gene_to_exon <- read.delim(file=paste("../../Genome/",annotation_name,".annotation.exon.gene.txt",sep=""), sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  trans_to_exon <- read.delim(file=paste("../../Genome/",annotation_name,".annotation.exon.transcr.txt",sep=""), sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  print("Gene to exons")
  summary(gene_to_exon$Length)
  png(filename = paste(figure_folder,"ExonsPerGene.png",sep=""), width = 1000, height = 1000,res = 150)
  hist(gene_to_exon$Length, main = "Number of exons per gene", xlab = "Nb of exons per gene", breaks = 100)
  dev.off()
  print("Transcript to exons")
  summary(trans_to_exon$Length)
  png(filename = paste(figure_folder,"ExonsPerTranscript.png",sep=""), width = 1000, height = 1000,res = 150)
  hist(trans_to_exon$Length, main = "Number of exons per transcript", xlab = "Nb of exons per transcript", breaks = 100)
  dev.off()

  assign("transcripts",transcripts,pos = globalenv())
  assign("exons",exons,pos = globalenv())
  assign("peaks_transcripts",peaks_transcripts,pos = globalenv())
  assign("peaks_exons",peaks_exons,pos = globalenv())

}
