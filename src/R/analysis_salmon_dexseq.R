library(tximport)
library(limma)
library(edgeR)
library(GenomicFeatures)
library(data.table)
library(readr)
library(DEXSeq)
library(stageR)
library(foreach)
library(dplyr)
library(gridExtra)
library(DRIMSeq)
library(BiocParallel)

# load MeRIPSeq scripts
general_path <- "/pasteur/projets/policy01/m6aAkker/"
rscript_folder <- paste(general_path,"MeRIPSeq/src/m6aAnalysis/",sep="")
pkgload::load_all(path = rscript_folder, reset = TRUE)

# bed_name <- 'MaxMaxValues_3'
# cores <- 6
dir_salmon <- file.path(general_path,'Salmon/')
gtf <- file.path(general_path,'Genome/gencode.vM13.annotation.gtf')

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 | args[1] == '-h') {
  print('Usage: Rscript analysis_salmon.R <expDesignName> <dirSalmon> <cores>')
}
exp_design_name <- args[1]
dir_salmon <- args[2]
cores <- args[3]
gtf <- args[4]

if (missing(exp_design_name) | is.null(exp_design_name) | exp_design_name == '') {
  message('The name of the experimental design is missing, we will use Cecum_all as default')
}


BPPARAM=MulticoreParam(workers=cores)

# read sample annotations
exp_design_file <- file.path(general_path,'ExpDesign',paste0(exp_design_name,'_exp_design.txt'))
samps <- fread(exp_design_file)
samps[, sample_id:=Input_dataset]
samps[, dataset:='dt2018']
samps[grep('2019',DataName), dataset:='dt2019']

if (exp_design_name == 'Cecum_all') {
  # we reassign Libray 1a 1b and 1c to the library batch 1 because DEXSeq has to be ran by comparison 
  samps$Library <- as.character(samps$Library)
  samps$Library[grep('^1',samps$Library)] <- '1'
  design <- model.matrix(~0+BioCond+Library+seq+dataset+Phase, data = samps)
  design <- design[,-grep('seq2019',colnames(design))]
  design <- design[,-grep("datasetdt2019",colnames(design))]
  design <- design[,-grep("PhasePhase_2019",colnames(design))]
} else if (exp_design_name == 'Liver_all_T13') { 
  design <- model.matrix(~0+BioCond+dataset+seq+Library, data = samps)
} else if (exp_design_name == 'Liver_all') {
  samps$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', samps$BioCond)))
  design <- model.matrix(~0+Time:BioCond+dataset+seq+Library, data = samps)
  colnames(design) = gsub(":",".",colnames(design))
  colnames(design) = gsub("Time|BioCond","",colnames(design))
}

colnames(design) <- gsub('BioCond', '', colnames(design))

samps[, group:=BioCond]
samps[, condition:=BioCond]

# build contrast
target_epi_fold <-  samps
load_contrast()
# contrasts <- combn(x = unique(samps$BioCond), m = 2)
# contrasts <- paste(contrasts[2,],"-",contrasts[1,],sep="")
# contr_transcripts_names <- contrasts #gsub('-','-BioCond',paste0('BioCond',contrasts))
contr_transcripts_names <- contr_Epi_names
contr_matrix_transcripts <- limma::makeContrasts(contrasts = contr_transcripts_names, levels=colnames(design))

if (exp_design_name == 'Liver_all') {
  samps$BioCond <- paste(samps$Time,samps$BioCond,sep='.')
}


# read salmon counts using data.table
#dir_salmon <- file.path(general_path,'salmon/')
files <- list.files(dir_salmon, pattern = 'quant.sf', full.names=T, recursive = T)
names(files) <- basename(gsub('/quant.sf','',files))
counts <- lapply(files, fread)
# gencode
counts <- cbind(counts[[1]], do.call(cbind, lapply(counts[-1], function(x) x[,-c(1:4)])))
colnames(counts)[-c(1:4)] <- basename(gsub('/quant.sf','',files))
countsdt <- copy(counts)

# follows: https://bioconductor.org/packages/devel/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html

tx2gene <- data.table(TXNAME=countsdt$Name, 
                      GENEID=gsub("\\|.*", "", gsub('.*ENSMUSG', 'ENSMUSG', countsdt$Name)))
tx2gene[, ntx:=length(TXNAME), by=GENEID]
txdf <- tx2gene
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance="scaledTPM", txOut=TRUE)
names(txi)


cts <- txi$counts[,samps$sample_id]

counts <- data.table(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)

counts <- counts[, c('gene_id', 'feature_id', samps$sample_id), with=FALSE]
counts <- counts[which(rowSums(cts) > 0),]



library(DRIMSeq)
d <- dmDSdata(counts=as.data.frame(counts), samples=as.data.frame(samps))


#n <- nrow(samps)
n <- 3
n.small <- 3
dall <- dmFilter(d,
                 min_samps_feature_expr=n.small, min_feature_expr=10,
                 min_samps_feature_prop=n.small, min_feature_prop=0.1,
                 min_samps_gene_expr=n, min_gene_expr=10)

# how many isoforms have the genes that have been kept?
table(table(counts(dall)$gene_id))


  #### ======= 1. Run DEXSeq + stageR analysis of salmon scaledTPM per comparison of interest ====== ####
  
  resdtusave <- list()
  dxrsave <- list()
  dxr.gsave <- list()
  treated <- c()
  
  resdtu <- foreach(comp = contr_transcripts_names, .errorhandling = 'stop') %do% { #[10:length(contr_transcripts_names)]
    
    treated <- c(treated, comp)
    # full and reduced model for DEXSeq, had to be adjusted per comparison because batches can be present or not depending on the comparison
    if (exp_design_name == 'Cecum_all') {
      if (comp %in%  c("ex_GF-Am" )) {
        fullModel <- ~ sample + exon + condition:exon + Library:exon + dataset:exon + Phase:exon
        fullModel_drimseq <- ~  group + Library + dataset + Phase
        reducedModel <- ~ sample + exon + Library:exon + dataset:exon + Phase:exon
      } else if (comp %in% c("Ec-vanco")) {
        fullModel <- ~ sample + exon + condition:exon  + seq:exon + dataset:exon + Phase:exon
        fullModel_drimseq <- ~  group + seq + dataset + Phase
        reducedModel <- ~ sample + exon + seq:exon + dataset:exon + Phase:exon
      }else if (comp %in% c('Ec-ex_GF','Lp-ex_GF',"abx_GF-ex_GF")) {
        fullModel <- ~ sample + exon + condition:exon  
        fullModel_drimseq <- ~ group
        reducedModel <- ~ sample + exon 
      } else if (comp %in% c("Lp-Ec","abx_GF-Lp", "abx_GF-Ec") ) {
        fullModel <- ~ sample + exon + condition:exon + Library:exon 
        fullModel_drimseq <- ~ group + Library
        reducedModel <- ~ sample + exon + Library:exon
      } else {
        fullModel <- ~ sample + exon + condition:exon + Library:exon + seq:exon + dataset:exon + Phase:exon
        fullModel_drimseq <- ~  group + Library + seq + dataset + Phase
        reducedModel <- ~ sample + exon + Library:exon + seq:exon + dataset:exon + Phase:exon
      } 
    } else if (exp_design_name == 'Liver_all_T13') {
      if (comp == "Liver_Lp_T13-Liver_Ec_T13") {
        fullModel <- ~ sample + exon + condition:exon + Library:exon
        fullModel_drimseq <- ~ group  + Library
        reducedModel <- ~ sample + exon + Library:exon
      } else {
        fullModel <- ~ sample + exon + condition:exon + Library:exon + seq:exon + dataset:exon
        fullModel_drimseq <- ~ group  + Library + seq + dataset
        reducedModel <- ~ sample + exon + Library:exon + seq:exon + dataset:exon
      }
    } else if (exp_design_name == 'Liver_all') {
      # if (!limma::is.fullrank(design)) {
      #   design_full <- design[,-which(colnames(design) %in% limma::nonEstimable(design))]
      # }
      if (comp %in% c("T3.Lp-T3.Ec", "T13.Lp-T13.Ec", "T13.Lp-T3.Lp", "T13.Ec-T3.Ec")) {
        fullModel <- ~ sample + exon + condition:exon + Library:exon
        fullModel_drimseq <- ~ group  + Library 
        reducedModel <- ~ sample + exon + Library:exon 
      } else {
        fullModel <- ~ sample + exon + condition:exon + Library:exon + seq:exon #+ dataset:exon
        fullModel_drimseq <- ~ group  + Library + seq #+ dataset
        reducedModel <- ~ sample + exon + Library:exon + seq:exon #+ dataset:exon
      }
      
    } else {
      fullModel <- ~ sample + exon + condition:exon
      fullModel_drimseq <- ~ group
      reducedModel <- ~ sample + exon
    }
    
    print(comp)
    comp <- strsplit(comp, split='-')[[1]]
    
    # ====== DRIMSeq filtering =====
    dallsel <- dmDSdata(counts=as.data.frame(counts(dall)[,c('gene_id', 'feature_id', sub('-T3', '.T3',sub('-T13','.T13', as.character(subset(samps, BioCond %in% comp)$Input_dataset))))]), 
                  samples=as.data.frame(subset(dall@samples, BioCond %in% comp) %>% mutate(sample_id = sub('-T3', '.T3',sub('-T13','.T13', as.character(sample_id))) )))

    print(table(table(counts(dallsel)$gene_id)))
    

    n <- nrow(dallsel@samples)
    n.small <- min(table(dallsel@samples$condition))

    d <- dmFilter(dallsel,
                  min_samps_feature_expr=n.small, min_feature_expr=5,
                  #min_samps_feature_prop=n.small, min_feature_prop=0.1,
                  min_samps_gene_expr=n, min_gene_expr=5)
    
    # how many isoforms have the genes that have been kept?
    print(table(table(counts(d)$gene_id)))
    
    #### ==== DEXSeq ===== ###
    sample.data <- DRIMSeq::samples(d)
    sample.data <- subset(sample.data, BioCond %in% comp)
    sample.data$condition <- sample.data$BioCond
    count.data <- round(as.matrix(counts(d)[,gsub('-','.',as.character(sample.data$sample_id))]))
    count.data <- count.data #[1:1000,]
    
    
    #designForDEXSeq <- model.matrix(fullModel, data = sample.data)
    dxd <- DEXSeqDataSet(countData=count.data,
                         sampleData=sample.data,
                         design=fullModel,
                         featureID=counts(d)$feature_id,
                         groupID=counts(d)$gene_id)
    
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, quiet=TRUE, BPPARAM=BPPARAM)
    dxd <- testForDEU(dxd, 
                      reducedModel=reducedModel,
                      fullModel=fullModel, BPPARAM=BPPARAM)
    dxd <-  estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
    #https://support.bioconductor.org/p/116893/
    
    dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
    dxror <- dxr
    # compute a per-gene adjusted p-value, using the perGeneQValue function, 
    # which aggregates evidence from multiple tests within a gene to a single p-value for the gene and then corrects for multiple testing across genes. 
    qval <- perGeneQValue(dxr)
    dxr.g <- data.frame(gene=names(qval),qval)
    # For size consideration of the workflow R package, we reduce also the transcript-level results table to a simple data.frame:
    columns <- c("featureID","groupID","pvalue")
    dxr <- as.data.frame(dxr[,columns])
    dxr <- subset(dxr, !is.na(pvalue))
    print(head(dxr[order(dxr$pvalue),]))
    print(sum(dxr$pvalue < 0.05, na.rm = TRUE))  
    print(sum(dxr.g$qval < 0.05))  
    
    if (sum(dxr.g$qval < 0.05) > 0) {
      # ===== stageR ==== #
      # A typical analysis of differential transcript usage would involve asking first: "which genes contain any evidence of DTU?", 
      # and secondly, "which transcripts in the genes that contain some evidence may be participating in the DTU?" 
      # Note that a gene may pass the first stage without exhibiting enough evidence to identify one or more transcripts that are participating in the DTU. 
      # The stageR package is designed to allow for such two-stage testing procedures, where the first stage is called a screening stage and the second stage a confirmation stage12. The methods are general, and can also be applied to testing, for example, changes across a time series followed by investigation of individual time points, as shown in the stageR package vignette. 
      # We show below how stageR is used to detect DTU and how to interpret its output.
      
      # The code above is a conventional DEXSeq analysis for analysing differential transcript usage. 
      # It would proceed by either assessing the significant genes according to the gene-wise q-values or by assessing the significant transcripts according to the transcript-level p-values, 
      # after adjustment for multiple testing. Performing and interpreting both analyses does not provide appropriate FDR control and thus should be avoided. 
      # However, interpretation on the gene level combined with transcript-level results can provide useful biological insights and this can be achieved through stage-wise testing.
      
      
      # We must specify an alpha, which will be the overall false discovery rate target for the analysis. 
      # Unlike typical adjusted p-values or q-values, we cannot choose an arbitrary threshold later: 
      # after specifying alpha=0.05, we need to use 5% as the target in downstream steps. 
      # There are also convenience functions getSignificantGenes and getSignificantTx, which are demonstrated in the stageR vignette.
      
      
      strp <- function(x) gsub('\\.', 'p', gsub('\\|.*|', '', x)) #substr(x,1,15)
      pConfirmation <- matrix(dxr$pvalue,ncol=1)
      dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
      #dimnames(pConfirmation) <- list((dxr$featureID),"transcript")
      pScreen <- qval
      names(pScreen) <- strp(names(pScreen))
      tx2gene <- (dxr[,c("featureID", "groupID")])
      for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
      #rownames(tx2gene) <- NULL
      colnames(tx2gene) <- c('transcripts','genes')
      stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, 
                            pScreenAdjusted=TRUE, tx2gene=tx2gene)
      message(paste0('Number of NA p-value is:', sum(is.na(pConfirmation))))
      stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=FALSE)
      suppressWarnings({
        dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                       onlySignificantGenes=TRUE)
      })
      #getSignificantGenes( stageRObj )
      #getSignificantTx( stageRObj )
      
      # The final table with adjusted p-values summarizes the information from the two-stage analysis.
      # Only genes that passed the filter are included in the table, so the table already represents screened genes. 
      # The transcripts with values in the column, transcript, less than 0.05 pass the confirmation stage on a target 5% overall false discovery rate, or OFDR. 
      # This means that, in expectation, no more than 5% of the genes that pass screening will either (1) not contain any DTU, 
      # so be falsely screened genes, or (2) contain a falsely confirmed transcript. 
      # A falsely confirmed transcript is a transcript with an adjusted p-value less than 0.05 which does not exhibit differential usage across conditions. 
      # The stageR procedure allows us to look at both the genes that passed the screening stage and the transcripts with adjusted p-values less than our target alpha, 
      # and understand what kind of overall error rate this procedure entails. This cannot be said for 
      # an arbitrary procedure of looking at standard gene adjusted p-values and transcript adjusted p-values, 
      # where the adjustment was performed independently.
      print(head(dex.padj))
      dex.padj$featureID <- dxr$featureID[match(paste0(dex.padj$txID, '|',dex.padj$geneID), gsub('\\|OTT.*', '', dxr$featureID))]
      
      #resdtusave <- list(resdtusave,dex.padj)
      dxrsave <- list(dxrsave,dxror)
      dxr.gsave <- list(dxr.gsave,dxr.g)
      names(resdtusave)[length(names(resdtusave))] <- 
        names(dxrsave)[length(names(dxrsave))] <- 
        names(dxr.gsave)[length(names(dxr.gsave))] <- 
        paste(comp[1],comp[2],sep='-')
      
      print(dim(dex.padj))
      print(sum(dex.padj$gene < 0.05 & dex.padj$transcript < 0.05))
    } else {
      dex.padj <- NA
    }

    save(resdtusave, dxrsave, dxr.gsave, file = paste0(dir_salmon,exp_design_name,'_resdtu_temp.rda'))
    return(list(stager=dex.padj, dexseq=dxror, dexseq_qval=dxr.g)) 
  } 
  
  names(resdtu) <- contr_transcripts_names
  
  # save DTU results results 
  save(resdtu,  dall, file = paste0(dir_salmon,exp_design_name,'_resdtu_no2019sel.rda'))
  
  #load(paste0(dir_salmon,exp_design_name,'_resdtu_no2019sel.rda'))
 

  #### ======= 2. Add delta of isoform fractions to stageR output   ====== ####
  # compute isoform fractions and add both the average isoform fraction in each condition, 
  # and the delta of isoform fraction between conditions to the stageR output
  allifs <- foreach(comp=names(resdtu)) %do% {
      resdexseq=resdtu[[comp]]$dexseq 
      resstager=resdtu[[comp]]$stager
      #resdexseq@sampleData$condition <-gsub('T13.','',resdexseq@sampleData$condition)      
      
      if (!is.na(resstager)) {
        conds <- strsplit(comp,split='-')[[1]]
        resstager$geneID <- sub('p','.',resstager$geneID)
        resstager$txID <- sub('p','.',resstager$txID)
        # compute tx fractions
        dtc <- as.data.table(resdexseq$countData)
        setnames(dtc, colnames(resdexseq$countData), as.character(resdexseq@sampleData$DataName))
        
        dtc[, geneID := gsub(':.*', '', rownames(resdexseq$countData))]
        dtc[, txID := gsub('\\|.*|.*:', '', rownames(resdexseq$countData))]
        dtc[, ntx := length(txID), by=geneID]
        dtprops <- cbind(dtc[, lapply(.SD,function(x) x=x/sum(x)), by='geneID',.SDcols=as.character(resdexseq@sampleData$DataName)],txID=dtc$txID, ntx=dtc$ntx)
        dtprops[, mean_IF_cond1:=mean(unlist(.SD)), by = c('geneID','txID','ntx'), .SDcols=intersect(colnames(dtprops),subset(resdexseq@sampleData,condition %in% conds[1])$DataName)]
        dtprops[, mean_IF_cond2:=mean(unlist(.SD)), by = c('geneID','txID','ntx'), .SDcols=intersect(colnames(dtprops),subset(resdexseq@sampleData,condition %in% conds[2])$DataName)]
        dtprops[, dIF := mean_IF_cond1 - mean_IF_cond2]
        setnames(dtprops, 'mean_IF_cond1', paste0('mean_IF_', conds[1]))
        setnames(dtprops, 'mean_IF_cond2', paste0('mean_IF_', conds[2]))
        
        # update stageR output with tx fractions
        resdtu[[comp]]$stager <- left_join(resstager, dtprops[,c('txID',paste0('mean_IF_', conds[1]),paste0('mean_IF_', conds[2]), 'dIF'), with=F], by = 'txID')
      } 
      }
    
  
  #### ======= 3. Intersect DTU with peaks and add a threshold on isoform fractions  ====== ####
  dIF_threshold <- 0.1
  
  # load results of the MeRIPSeq pipeline for the experiment under consideration
  
  load(file.path(general_path, 'PeakDiffExpression', paste(exp_design_name, bed_name, sep = '_'),
                 paste0(paste(exp_design_name, bed_name, sep = '_'),'.RData')))

  # nb intersection between genes detected as DE and DTU
  lapply(resdtu, function(x) {xx = x$stager; if (!is.na(xx) ) length(intersect(unique(subset(xx,transcript < 0.05 & abs(dIF) >= dIF_threshold)$geneID), diff_genes)) })
  lapply(resdtu, function(x) {xx = x$stager; if (!is.na(xx) ) length(intersect(unique(subset(xx,gene < 0.05 & abs(dIF) >= dIF_threshold)$geneID), diff_genes)) })

  # load peak annotations
  peaksgr <- copy(peaks)
  setnames(peaksgr, c("chromo_peak",   "begin_peak",    "end_peak", 'chr'), c('chr', 'start', 'end','chr_gene'))
  peaksgr <- makeGRangesFromDataFrame(peaksgr, keep.extra.columns = TRUE, seqnames.field="chr")
  
  # load genes and transcript annotations
  txdb <- makeTxDbFromGFF(gtf, format = 'gtf')
  txgr <- transcripts(txdb)
  exonsgr <- exonsBy(txdb, by = 'tx')
  names(exonsgr) <- AnnotationDbi::select(txdb, keys = names(exonsgr), columns="TXNAME", keytype="TXID")$TXNAME
  
  # associate each peak to overlapping transcripts
  peaks_to_tx <- as.data.table(findOverlaps(query = exonsgr, subject = peaksgr, ignore.strand = TRUE, type = 'any'))
  peaks_to_tx[, peak := names(peaksgr)[subjectHits]]
  peaks_to_tx[, tx := names(exonsgr)[queryHits]]
  peak_to_tx <- peaks_to_tx
  
  # how many differentially methylated peaks are not overlapping a gencode transcript?
  length(setdiff(diff_peaks, peak_to_tx$peak))
  
  # for each comparison, we are going to associate the differentially methylated peaks with information 
  # regarding differential transcript usage of the transcripts they overlap.
  
  DMpeaksWithDTU <- 
    lapply(contr_Epi_names, function(i, peak_to_tx) {
      print(i)
      comp <- strsplit(i, split = '-')[[1]]
      comps <- c(paste(comp[1],comp[2],sep='-'), paste(comp[2],comp[1],sep='-'))
      # differentially expressed genes for this comparison (as found with Christophe, with voom, on read counts, not ideal it should have been done with tximport for consistency)
      resDEgenes <- limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = 1, lfc = 0, number=Inf)
      DEgenes <- rownames(limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf))
      # differentially methylated peaks for this comparison
      DMpeaks <- rownames(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf));
      # differentially used transcripts for this comparison
      # get contrast in case it was performed in the opposite direction (we don't use any log fold change yet so not a problem but could be later)
      compTx <- comps[which(comps %in% names(resdtu))] 
      
      resStageR <- resdtu[[compTx]]$stager
      if (!is.na(resStageR)) {
        
        # store featureIDs
        tx2featureID <- structure(.Data=rownames(resdtu[[compTx]]$dexseq),
                                  .Names=gsub('\\|.*|.*:', '', rownames(resdtu[[compTx]]$dexseq)))
        
        # reassign full featureIDs and substitute the 'p' that was introduced to solve stageR issues
        resStageR <- resStageR %>%  
          mutate(featureID=tx2featureID[txID]) %>%
          mutate(txID=gsub('p', '.',txID)) %>%
          mutate(geneID=gsub('p', '.',geneID)) 
        
        resStageR <- resStageR[abs(resStageR$dIF) >= dIF_threshold,]
        
        # extract differentially used transcripts and matching genes
        DTUtx <- subset(resStageR, gene < 0.05)$txID #transcript < 0.05)$txID
        DTUgenes <- unique(subset(resStageR, transcript < 0.05)$geneID)
        #select(txdb, keys = names(exonsgr), columns="TXNAME", keytype="TXID")$TXNAME
        
        # extract transcripts overlapping with differentially methylated peaks
        txinDMpeaks <- peak_to_tx[peak %in% DMpeaks,]$tx
        # subset results for those and add peak  and gene info
        txinDMpeaks_info <- resStageR[resStageR$txID %in% txinDMpeaks,]
        genesinDMpeaks_info <- genes[genes$Gene_ID %in% txinDMpeaks_info$geneID,] %>% rename(., geneID = Gene_ID)
        DTUinDMpeaks_info <- left_join(txinDMpeaks_info, genesinDMpeaks_info, by = 'geneID')
        
        # boxplots
        #boxplot(log(resdtu[[compTx]]$dexseq$countData[resStageR[resStageR$txID %in% txinDMpeaks,]$featureID[1],])~resdtu[[compTx]]$dexseq@sampleData$BioCond)
        # extract differentially methylated peaks associated with DTU and the overlapping transcripts for each peak
        DMpeaksWithDTU <- intersect(peak_to_tx[tx %in% as.character(DTUtx),]$peak, DMpeaks)
        if (length(DMpeaksWithDTU) > 0) {
          
          DMpeaksWithDTU_inDEgenes <- intersect(DMpeaksWithDTU, peaks_to_gene[peaks_to_gene$Gene_ID %in% DEgenes,]$ID)
          DMpeaksWithDTU_notInDEgenes <- setdiff(DMpeaksWithDTU, peaks_to_gene[peaks_to_gene$Gene_ID %in% DEgenes,]$ID)
          
          DTUinDMpeaks_info <- left_join(DTUinDMpeaks_info,
                                         peak_to_tx[peak %in% DMpeaksWithDTU & tx %in%txinDMpeaks_info$txID, list(peaks=list(peak)), by='tx'] %>% setnames(., 'tx', 'txID'),
                                         by = 'txID')
          DTUinDMpeaks_info <- DTUinDMpeaks_info %>% rename(DTU_gene_qval = gene) %>% rename(DTU_tx_OFDR = transcript)
          DTUinDMpeaks_info$DGE_gene_FDR = resDEgenes[DTUinDMpeaks_info$geneID,]$adj.P.Val
          
          # merge peak info if a transcript is associated with several peaks
          DTUinDMpeaks_info <- left_join(DTUinDMpeaks_info,
                                         as.data.frame(do.call(rbind,mapply(x=DTUinDMpeaks_info$peaks, y=DTUinDMpeaks_info$txID,
                                                                            function(x, y) {
                                                                              xcollapse <- apply(peaks[unlist(x),c("begin_peak", "end_peak", "length_peak",  "Type_overlap", "strand")] %>% 
                                                                                                   rename(peak_strand=strand),
                                                                                                 2,paste,collapse='|')
                                                                              
                                                                              xcollapse <- c(xcollapse,txID=y)
                                                                              return(xcollapse)
                                                                              
                                                                            },
                                                                            SIMPLIFY = FALSE))), 
                                         by = 'txID')
          
          
          DTUinDMpeaks_info$DGE_gene_status <- 'notDE'
          DTUinDMpeaks_info$DGE_gene_status[DTUinDMpeaks_info$DGE_gene_FDR < 0.05] <- 'DE'
          DTUinDMpeaks_info <- DTUinDMpeaks_info[order(DTUinDMpeaks_info$DTU_gene_qval,DTUinDMpeaks_info$DTU_tx_OFDR),]
          DTUinDMpeaks_info %>% mutate(DTU_dIF=dIF) 
          colnames(DTUinDMpeaks_info)[grep('mean_IF_', colnames(DTUinDMpeaks_info))] <- paste0('DTU_',colnames(DTUinDMpeaks_info)[grep('mean_IF_', colnames(DTUinDMpeaks_info))])
          DTUinDMpeaks_info$gene_name[is.na(DTUinDMpeaks_info$gene_name)] <- DTUinDMpeaks_info$geneID[is.na(DTUinDMpeaks_info$gene_name)]
          return(DTUinDMpeaks_info[,c("peaks", "geneID", "txID", "gene_name", "begin_peak", "end_peak",  
                                      colnames(DTUinDMpeaks_info)[grep('DTU_', colnames(DTUinDMpeaks_info))],
                                      #"DTU_", "DTU_", "DTU_", "DTU_gene_qval", "DTU_tx_OFDR", 
                                      "DGE_gene_FDR", "DGE_gene_status", "featureID", 
                                      "chr", "gene_begin", "gene_end", "strand", "type", 
                                      "EntrezID", "UniprotID", "UniprotIDs", "Description",
                                      "length_peak", "Type_overlap", "peak_strand")])      
        } else{
          return(NULL)
        }
      } else {
        return(NULL)
      }
      
      
    }, peak_to_tx = peak_to_tx)
  
  
  names(DMpeaksWithDTU) <- contr_Epi_names
  
  
  dir_dtu <- file.path(general_path, 'PeakDiffExpression', 
                       paste(exp_design_name, bed_name, sep = '_'),
                       'DTU_salmon')
  dir.create(dir_dtu)
  
    
  
  ### === Plot number of genes with evidence of DTU per comparison and the number of those that are oevrlapping a peak found to be differentially methylated === ###
  nbDTU_percond_overlappingDM <- unlist(lapply(DMpeaksWithDTU, function(x) if (!is.null(x)) length(unique(x$geneID)) else 0))
  nbDTU_percond_overlappingDM_notDE <- unlist(lapply(DMpeaksWithDTU, function(x) if (!is.null(x)) length(unique(subset(x, DGE_gene_status == 'notDE')$geneID)) else 0))
  nbDTU_percond_overlappingDM_notDE_genes <- unlist(lapply(DMpeaksWithDTU, function(x) if (!is.null(x)) paste(unique(subset(x, DGE_gene_status == 'notDE')$gene_name),collapse=',') else ""))
  nbDTU_percond_overlappingDM_DE <- unlist(lapply(DMpeaksWithDTU, function(x) if (!is.null(x)) length(unique(subset(x, DGE_gene_status == 'DE')$geneID)) else 0))

  # nb significant DTU
  nbdtu_percomp <- sapply(resdtu, function(x) {xx = x$stager; if (!is.na(xx) ) length(unique(subset(xx,transcript < 0.05 & abs(dIF) >= dIF_threshold)$geneID)) else 0})
  names(nbdtu_percomp)[which(! (names(nbdtu_percomp) %in% names(nbDTU_percond_overlappingDM_DE)))] <- sapply(names(nbdtu_percomp)[which(! (names(nbdtu_percomp) %in% names(nbDTU_percond_overlappingDM_DE)))], function(x) {conds=strsplit(x,split='-')[[1]]; paste(conds[2],conds[1],sep='-')})
  nbdtu_percomp <- unlist(nbdtu_percomp)
  nbdtu_percomp <- nbdtu_percomp[names(nbDTU_percond_overlappingDM_DE)]
  dtnbdtu <- data.table(comp=names(nbDTU_percond_overlappingDM_DE),
                        nbDTU_ovDMpeak_DEgene=nbDTU_percond_overlappingDM_DE, 
                        nbDTU_ovDMpeak_notDEgene=nbDTU_percond_overlappingDM_notDE) 
  dtnbdtu$nbDTU_notovDMpeak <- unlist(nbdtu_percomp)-(dtnbdtu$nbDTU_ovDMpeak_DEgene+dtnbdtu$nbDTU_ovDMpeak_notDEgene)
  # add
  dtnbdtu_summary <- data.table(comparison = names(nbdtu_percomp))
  dtnbdtu_summary$nbDTUgenes <- unlist(nbdtu_percomp)
  dtnbdtu_summary$nbDEgenes <- sapply(names(nbdtu_percomp), function(i) nrow(limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))
  dtnbdtu_summary$nbDMpeaks <- sapply(names(nbdtu_percomp), function(i) nrow(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))
  dtnbdtu_summary$nbGenes_withDMpeak <- sapply(names(nbdtu_percomp), function(i) length(unique(peaks_to_gene[rownames((limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf))),]$Gene_ID)))
  dtnbdtu_summary$nbDMpeaks_withDEgene <- sapply(names(nbdtu_percomp), 
                                                 function(i) {
                                                   dmpeaks <- unique(rownames(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))
                                                   sum(peaks_to_gene[match(dmpeaks, peaks_to_gene$ID),]$Gene_ID %in% (rownames((limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))))
                                                 })
  # peaks overlapping a DE gene, across conditions
  DMPeaks_withDEgene <- sapply(names(nbdtu_percomp), 
                                  function(i) {
                                    dmpeaks <- unique(rownames(limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))
                                    dmpeaks[(peaks_to_gene[match(dmpeaks, peaks_to_gene$ID),]$Gene_ID %in% unique(rownames((limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))))]
                                  })
  dtnbdtu_summary$nbDEgenes_withDMpeak <- sapply(names(nbdtu_percomp), 
                                          function(i) length(intersect(unique(peaks_to_gene[rownames((limma::topTable(fit_bayes_epi, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf))),]$Gene_ID),
                                                                       unique(rownames((limma::topTable(fit_bayes_RNA_I, coef = i, adjust="BH", p.value = .05, lfc = 1, number=Inf)))))))
  dtnbdtu_summary$nbDTUgenes_withDMpeak <- nbDTU_percond_overlappingDM
  dtnbdtu_summary$nbDTUgenes_withDMpeak_withDEgene <- nbDTU_percond_overlappingDM_DE
  dtnbdtu_summary$nbDTUgenes_withDMpeak_withoutDEgene <- nbDTU_percond_overlappingDM_notDE
  dtnbdtu_summary$nbDTUgenes_withDMpeak_withoutDEgene_list <- nbDTU_percond_overlappingDM_notDE_genes

  if (exp_design_name == 'Cecum_all') {
    dtnbdtu_summary <- dtnbdtu_summary[c(grep("^GF-CONV",dtnbdtu_summary$comparison), 
                                         grep("CONV-abx$",dtnbdtu_summary$comparison), 
                                         setdiff(grep("CONV",dtnbdtu_summary$comparison), grep("^GF-CONV|CONV-abx$",dtnbdtu_summary$comparison,invert=FALSE)), 
                                         setdiff(grep("^GF|-GF$",dtnbdtu_summary$comparison), grep("CONV|^GF-CONV|CONV-abx$",dtnbdtu_summary$comparison,invert=FALSE)),
                                         setdiff(order(dtnbdtu_summary$comparison), grep("-GF$|^GF|CONV|^GF-CONV|CONV-abx$",dtnbdtu_summary$comparison,invert=FALSE))),]
  } else if (exp_design_name == 'Liver_all') {
    dtnbdtu_summary <- dtnbdtu_summary[c(grep("GF.*CONV*|CONV*GF*",dtnbdtu_summary$comparison), 
                                         setdiff(grep("CONV",dtnbdtu_summary$comparison), grep("GF.*CONV*|CONV*GF*",dtnbdtu_summary$comparison,invert=FALSE)), 
                                         setdiff(grep("GF",dtnbdtu_summary$comparison), grep("CONV|GF.*CONV*|CONV*GF*",dtnbdtu_summary$comparison,invert=FALSE)),
                                         setdiff(order(dtnbdtu_summary$comparison), grep("GF|CONV|CONV.*abx.*|abx.*GF.*|GF.*CONV*|CONV*GF*",dtnbdtu_summary$comparison,invert=FALSE))),]
  }
    
  ggnbdtu <- dtnbdtu %>% melt(.) %>% ggplot(., aes(x=comp,y=value,fill=variable)) + 
    geom_bar(stat = 'identity', position='dodge') + 
    theme(axis.text.x = element_text(angle=90), legend.title = element_blank()) +
    ylab('#genes with DTU') + xlab('') 
  
  dir.create(file.path(dir_dtu, 'plots'))
  ggsave(plot=ggnbdtu, filename = file.path(dir_dtu, 'plots',paste0('nbDTUs','.pdf')), width=12, height=5)
  fwrite(x=dtnbdtu, file=file.path(dir_dtu, 'plots',paste0('nbDTUs','.txt')), sep='\t')
  
  # write DM peaks intersecting DTU
  DMpeaksWithDTU <- lapply(DMpeaksWithDTU, function(x) if(is.null(x)) return(data.frame(result=c('no intersection between differentially used transcripts and differentially methylated peaks',''))) else return(x))
  DMpeaksWithDTU[['DTU_summary']] <- dtnbdtu_summary
  WriteXLS::WriteXLS(DMpeaksWithDTU[c('DTU_summary', setdiff(names(DMpeaksWithDTU),'DTU_summary'))],#[!(sapply(DMpeaksWithDTU,is.null))], 
                     ExcelFileName =  file.path(dir_dtu,'DTUinDMpeak.xls'))
  
  
  
  #### ======= 3. Plot boxplot of tx proportions for DTU detected together with differentially methylated peaks  ====== ####

  marrangeGrob <- function (grobs, ..., ncol, nrow, layout_matrix = matrix(seq_len(nrow * ncol), nrow = nrow, ncol = ncol, byrow=TRUE), 
                            top = quote(paste("page", g, "of", npages))) 
  {
    n <- length(grobs)
    nlay <- max(layout_matrix, na.rm = TRUE)
    npages <- n%/%nlay + as.logical(n%%nlay)
    groups <- split(grobs, rep(seq_len(npages), each = nlay, 
                               length.out = n))
    pl <- vector(mode = "list", length = npages)
    for (g in seq_along(groups)) {
      params <- modifyList(list(...), list(top = eval(top), 
                                           layout_matrix = layout_matrix))
      pl[[g]] <- do.call(gridExtra::arrangeGrob, c(list(grobs = groups[[g]]), 
                                                   params))
    }
    class(pl) <- c("arrangelist", class(pl))
    pl
  }

  
  # for each comparison, plot proportions of tx of a gene where at least one tx overlaps with a differentially methylated peak

  foreach(comp=names(resdtu)) %do% {
    
    print(comp)
    resdexseq <- resdtu[[comp]]$dexseq
    if (!(comp %in% names(DMpeaksWithDTU))) compdtu <- paste(strsplit(comp,'-')[[1]][2],strsplit(comp,'-')[[1]][1],sep='-')  else  compdtu <- comp
    dtu <- DMpeaksWithDTU[[compdtu]]
    
    if (as.character(dtu[1,1]) != 'no intersection between differentially used transcripts and differentially methylated peaks') {
      # for each gene of the table, plot its tx proportions
      # i recomputed the props here bc i did that part before i decided to compute all isoform fractions 
      boxplots <- foreach( gene = unique(dtu$geneID), genename=unique(dtu$gene_name)) %do% {
        print(gene)
        txs <- gsub('.*:|\\|.*', '', rownames(resdexseq$countData[grep(paste(gene, collapse='|'),rownames(resdexseq$countData)),,drop=F]))
        props <- scale(resdexseq$countData[grep(gene,rownames(resdexseq$countData)),,drop=F], 
                       scale=colSums(resdexseq$countData[grep(gene,rownames(resdexseq$countData)),,drop=F]),
                       center=FALSE)
        colnames(props) <- resdexseq@sampleData$DataName
        rownames(props) <- gsub('.*:|\\|.*', '', rownames(props))
        melt(props) %>% rename(tx=Var1, sampleId=Var2, tx_prop=value) %>% mutate(condition=resdexseq@sampleData[match(sampleId,resdexseq@sampleData$DataName),]$BioCond) %>% 
          ggplot(., aes(x=condition, y=tx_prop, fill=tx)) + geom_boxplot() + ggtitle(paste(gene, genename, sep='\n')) + theme(legend.position = 'right') + ylab('tx_frac')
      }
      gge <- marrangeGrob(grobs=boxplots, ncol=3, nrow = 2)
      ggsave(plot=gge, filename = file.path(dir_dtu, 'plots',paste0(comp,'.pdf')), width = 14)
    }
  } 
    