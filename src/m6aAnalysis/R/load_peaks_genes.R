load_peaks_genes <- function(){

  # check peak distribution
  print("Load Peaks")
  peaks <- read.delim(file=paste("../../PeakDetection/Peaks/",peak_name,"_Peaks.txt",sep = ""), sep = "\t", row.names = 1, stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(peaks),"peaks loaded"))
  peaks$ID = row.names(peaks)
  peaks_to_gene = peaks[,c("ID","Gene_ID")]
  peaks_to_gene$ID = row.names(peaks_to_gene)

  # if MAX change peaks columns
  # peaks$max = peak_to_max[rownames(peaks),"MaxIndex"]
  # peaks$begin_peak = peaks$max - 50
  # peaks$end_peak = peaks$max + 50
  #peaks$max = NULL


  # Load annotation description
  print("Load Genes")
  gene_filename <- paste("../../Genome/",annotation_name,".annotation.gene.txt",sep="")
  print(gene_filename)
  genes <- read.delim(file=gene_filename, sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  print(paste(nrow(genes),"genes loaded"))
  genes$Gene_ID = genes$gene_id
  genes$gene_id <- NULL

  # Load peak sequence
  print("Load peak sequence")
  fasta_filename <- paste("../../PeakDetection/Peaks/",peak_name,".fasta",sep="")
  print(fasta_filename)
  fasta <- seqinr::read.fasta(fasta_filename)
  # Load fasta for peaks
  folderFasta <- paste("Motif/fasta/",sep="")
  if (!dir.exists(folderFasta)){
    dir.create(folderFasta)
  }
  norm.filename = paste(folderFasta,peak_name,".fasta",sep = "")
  seqinr::write.fasta(sequences=fasta,names=names(fasta),file.out=norm.filename)

  # Load peak bed file and transfrom it for guitar plot
  bed_filename <- paste("../../PeakDetection/Peaks/",peak_name,".bed",sep="")
  bed_file <- read.delim(bed_filename, header = FALSE)
  bed_file$V5 = 1
  new_bed_filename <- paste(path_guitar,"bed/",peak_name,".bed",sep = "")
  write.table(bed_file, new_bed_filename, col.names = F, row.names = F, quote = F, sep="\t")


  assign("fasta",fasta,pos = globalenv())
  assign("peaks",peaks,pos = globalenv())
  assign("peaks_to_gene",peaks_to_gene, pos = globalenv())
  assign("genes",genes,pos = globalenv())

}
