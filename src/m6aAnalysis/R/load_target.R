  loadTargetFile <- function (targetFile, varInt, condRef, batch)
{
  target <- read.table(targetFile, header = TRUE, sep = "\t",
                       na.strings = "")
  if (!I(varInt %in% names(target)))
    stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target)))
    stop(paste("The batch effect", batch, "is not in the target file"))
  target[, varInt] <- as.factor(target[, varInt])
  if (!I(condRef %in% as.character(target[, varInt])))
    stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  target[, varInt] <- relevel(target[, varInt], ref = condRef)
  target <- target[order(target[, varInt]), ]
  rownames(target) <- as.character(target[, 1])
  if (min(table(target[, varInt])) < 2)
    stop(paste("The factor of interest", varInt, "has a level without replicates"))
  if (any(is.na(cbind(target[, c(varInt, batch)], target[,
                                                         1:2]))))
    stop("NA are present in the target file")
  if (!is.null(batch) && is.numeric(target[, batch]))
    warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  #if (any(grepl("[[:punct:]]", as.character(target[, varInt]))))
  #  stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  #cat("Target file:\n")
  #print(target)
  return(target)
}

load_target <- function(){

  # Create appropriate target file
  print("Create target file from exp_design file")
  print(exp_design_file)
  exp_design_table <- read.delim(file = exp_design_file, sep = "\t", stringsAsFactors = FALSE, header =  TRUE)
  labels = c(exp_design_table$IP_dataset,exp_design_table$Input_dataset)
  files = paste("HTSeq_",exp_design_name,"_",labels,"_",bed_name,"_peaks.txt",sep="")

  files_Input = paste("HTSeq_",exp_design_table$Input_dataset,".txt",sep="")
  files_IP = paste("HTSeq_",exp_design_table$IP_dataset,".txt",sep="")

  #List of parameters accessible in the model
  #list_batch=c("BioCond","seq","Seq_Phase","Library","Mice","Exp","Combined_batch","lane","Tissue")
  DataName = c(exp_design_table$DataName,exp_design_table$DataName)
  TypeData = c(rep("IP", nrow(exp_design_table)),rep("Input", nrow(exp_design_table)))
  BioCond = c(exp_design_table$BioCond,exp_design_table$BioCond)
  Seq = c(as.character(exp_design_table$seq),as.character(exp_design_table$seq))
  Phase=c(exp_design_table$Phase,exp_design_table$Phase)
  Library=c(exp_design_table$Library,exp_design_table$Library)
  Mice=c(exp_design_table$Mice,exp_design_table$Mice)
  Exp = c(exp_design_table$experiment,exp_design_table$experiment)
  Comb_batch=c(exp_design_table$Combined_batch,exp_design_table$Combined_batch)
  Lane = c(exp_design_table$lane,exp_design_table$lane)
  Tissue = c(exp_design_table$Tissue,exp_design_table$Tissue)
  RNAprep = c(exp_design_table$RNAPrep,exp_design_table$RNAPrep)
  mRNAprep = c(exp_design_table$mRNAprep,exp_design_table$mRNAprep)
  Library = as.character(c(exp_design_table$Library,exp_design_table$Library))
  LibraryIP = as.character(c(exp_design_table$libraryprepIP,exp_design_table$libraryprepIP))
  LibraryIPbis = c(exp_design_table$IP,exp_design_table$IP)
  #FlowCell = c(exp_design_table$FlowCell,exp_design_table$FlowCell)
  #FlowCell = c(sub('-.*', '', exp_design_table$IP_raw_dataset),sub('-.*', '', exp_design_table$IP_raw_dataset))
  profileIP <- c(exp_design_table$profileIP,exp_design_table$profileIP)
  RNAprepdate <- paste0('date_',gsub('-..T.*|-', '', c(exp_design_table$RNAprepdate,exp_design_table$RNAprepdate)))
  frag <- paste0('date_',gsub('-..T.*|-', '', c(exp_design_table$fragmentation,exp_design_table$fragmentation)))
  if ('Time' %in% colnames(exp_design_table)) {
    Time = c(exp_design_table$Time,exp_design_table$Time)
    target_table_epi <- data.frame(labels,files,DataName,BioCond,TypeData,Seq,Mice,Lane,Tissue,Time,LibraryIPbis,RNAprep)#, LibraryIPbis, RNAprep, mRNAprep)
  } else {
    target_table_epi <- data.frame(labels,files,DataName,BioCond,TypeData,Seq,Mice,Library,Phase)#Tissue, , Phase,RNAprepdate, mRNAprep)#,Phase)# LibraryIPbis, RNAprepdate, mRNAprep)#, Phase) #LibraryIPbis, LibraryIP, profileIP, RNAprepdate, frag, FlowCell, Exp)#LibraryIP,Exp,Library,Comb_batch,RNAprep,mRNAprep,
  }


  if (exp_design_name %in% c('Cecum_all','Liver_all', 'Liver_all_T13')) {
    # add dataset annotation
    target_table_epi$dataset <- ''
    target_table_epi[grep('2019',target_table_epi$labels, invert = TRUE),]$dataset <- 'dt1'
    target_table_epi[grep('2019',target_table_epi$labels),]$dataset <- 'dt2'
  }

  target_table_Input <- target_table_epi[target_table_epi$TypeData=="Input",]
  target_table_Input$files <- files_Input
  target_table_IP <- target_table_epi[target_table_epi$TypeData=="IP",]
  target_table_IP$files <- files_IP


  # save target files
  target_file_Epi <- paste(exp_design_name,"_target.txt",sep = "")
  target_file_RNA_Input <- paste(exp_design_name,"_target_Input.txt",sep = "")                           # path to the design/target file
  target_file_RNA_IP <- paste(exp_design_name,"_target_IP.txt",sep = "")                           # path to the design/target file
  write.table(target_table_epi, file=target_file_Epi, quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(target_table_Input, file=target_file_RNA_Input, quote=FALSE, row.names = FALSE, sep = "\t")
  write.table(target_table_IP, file=target_file_RNA_IP, quote=FALSE, row.names = FALSE, sep = "\t")

  condRef <- as.character(target_table_epi$BioCond[[1]])
  #batch <- "Seq"                                # blocking factor: NULL (default) or "batch" for example
  #varInt <- "BioCond"      # factor of interest

  # loading target file
  print("Load Target file for Peak expression and Gene expression")
  target_epi <- loadTargetFile(targetFile=target_file_Epi, varInt=varInt, condRef=condRef, batch=batch)

  target_epi_fold <- target_epi[target_epi$TypeData == "IP",]
  #target_epi_fold <- target_epi_fold[,c("labels","BioCond","Library","Phase","Seq","Exp","DataName")]
  #target_epi_fold <- target_epi_fold[,c("labels","BioCond","Seq","Batch","DataName")]
  target_epi_fold$IP <- target_epi_fold$labels
  target_epi_fold$Input <- sub("IP","Input",target_epi_fold$labels) #anne: replace Input by input
  target_epi_fold$labels <- NULL

  target_epiIP <- target_epi[grep('IP',rownames(target_epi)),]

  row.names(target_epi_fold) = target_epi_fold$DataName
  target_RNA_I <- loadTargetFile(targetFile=target_file_RNA_Input, varInt=varInt, condRef=condRef, batch=batch)
  target_RNA_IP <- loadTargetFile(targetFile=target_file_RNA_IP, varInt=varInt, condRef=condRef, batch=batch)


  print("Load models for MeRIPSeq and RNASeq")
  if(exp_design_name=="Cecum_all"){
    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond+target_epi$dataset+target_epi$Library+target_epi$Seq+target_epi$Phase)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0+target_epiIP$TypeData:target_epiIP$BioCond+target_epiIP$Library+target_epiIP$Seq+target_epiIP$dataset+target_epiIP$Phase)
    rownames(model_epi_bis) <- as.character(target_epiIP$labels)[grep('IP',as.character(target_epiIP$labels))]
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    model_epi_bis <- model_epi_bis[,c(grep('Mice',colnames(model_epi_bis), invert = TRUE),grep('Mice',colnames(model_epi_bis)))]
    colnames(model_epi_bis) = gsub("\\[.*|target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeDataIP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    model_epi_bis <- model_epi_bis[grep('IP',rownames(model_epi_bis)),]
    model_epi <- model_epi_bis
    model_epi <- model_epi[,ncol(model_epi):1]
    rownames(model_epi ) <- as.character(target_epiIP$labels)[grep('_IP',as.character(target_epiIP$labels))]

    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond+target_epi_fold$Library+target_epi_fold$Seq+target_epi_fold$dataset+target_epi_fold$Phase)
    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond+target_RNA_I$Library+target_RNA_I$Seq+target_RNA_I$dataset+target_RNA_I$Phase)
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond+target_RNA_IP$Library+target_RNA_IP$Seq+target_RNA_IP$dataset+target_RNA_IP$Phase)

  }else  if(exp_design_name=="Cecum_2019"){
    target_epi_fold$Input <- sub("input","Input",target_epi_fold$Input)
    target_epi$Seq <- paste0('Seq',target_epi$Seq)
    target_epi$RNAprep <- paste0('prep',target_epi$RNAprep)
    target_epi$RNAprepdate <- paste0('prep',target_epi$RNAprepdate)
    target_epi$LibraryIP <- paste0('prep',target_epi$LibraryIP)
    target_epiIP <- target_epi[grep('IP',rownames(target_epi)),]

    target_epi_fold$Seq <- paste0('Seq',target_epi_fold$Seq)
    target_epi_fold$RNAprep <- paste0('prep',target_epi_fold$RNAprep)
    target_RNA_I$Seq <- paste0('Seq',target_RNA_I$Seq)
    target_RNA_I$RNAprep <- paste0('prep',target_RNA_I$RNAprep)
    target_RNA_I$mRNAprep <- paste0('prep',target_RNA_I$mRNAprep)
    target_RNA_IP$RNAprep <- paste0('prep',target_RNA_I$RNAprep)
    target_RNA_IP$Seq <- paste0('Seq',target_RNA_I$Seq)


    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond+target_epi$LibraryIPbis+target_epi$RNAprepdate)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0+target_epiIP$TypeData:target_epiIP$BioCond+target_epiIP$LibraryIPbis+target_epiIP$RNAprepdate)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)[grep('IP',as.character(target_epiIP$labels))]
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    model_epi_bis <- model_epi_bis[,c(grep('Mice',colnames(model_epi_bis), invert = TRUE),grep('Mice',colnames(model_epi_bis)))]
    colnames(model_epi_bis) = gsub("\\[.*|target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeDataIP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    model_epi_bis <- model_epi_bis[grep('IP',rownames(model_epi_bis)),]
    model_epi <- model_epi_bis
    model_epi <- model_epi[,ncol(model_epi):1]

    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond+target_epi_fold$LibraryIPbis+target_epi_fold$RNAprepdate)
    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond+target_RNA_I$LibraryIPbis+target_RNA_I$RNAprepdate)
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond)

  }else if(exp_design_name %in% c("Liver_all")){

    target_epi$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', target_epi$BioCond)))
    target_epiIP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', target_epiIP$BioCond)))
    target_epi_fold$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', target_epi_fold$BioCond)))
    target_RNA_I$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', target_RNA_I$BioCond)))
    target_RNA_IP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', target_RNA_IP$BioCond)))

    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$Time:target_epi$BioCond+target_epi$dataset+target_epi$Seq+target_epi$Library)# + target_epi$Seq)# + target_epi$RNAprep)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("Time","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 +target_epiIP$TypeData:target_epiIP$Time:target_epiIP$BioCond+target_epiIP$dataset+target_epiIP$Seq+target_epiIP$Library)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("Time","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))

    model_epi_bis <- model_epi_bis[,-grep('seq2019',colnames(model_epi_bis))] # remove duplicates from the design matrix, the 2019 sequencing was done in one batch
    model_epi_bis <- model_epi_bis[,ncol(model_epi_bis):1]
    model_epi <- model_epi_bis

    model_epi_fold <- model.matrix(~0+target_epi_fold$Time:target_epi_fold$BioCond+target_epi_fold$dataset+target_epi_fold$Seq+target_epi_fold$Library)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))
    model_epi_fold <- model_epi_fold[,-grep('seq2019',colnames(model_epi_fold))]
    model_epi_fold <- model_epi_fold[,ncol(model_epi_fold):1]

    model_RNA_I <- model.matrix(~0+target_RNA_I$Time:target_RNA_I$BioCond+target_RNA_I$dataset+target_RNA_I$Seq+target_RNA_I$Library)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_I <- model_RNA_I[,-grep('seq2019',colnames(model_RNA_I))]

    model_RNA_IP <- model.matrix(~0+target_RNA_IP$Time:target_RNA_IP$BioCond+target_RNA_IP$dataset+target_RNA_I$Seq+target_RNA_I$Library)# + target_RNA_IP$Seq)# + target_RNA_IP$RNAprep)# + target_RNA_IP$Library)
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))
    model_RNA_IP <- model_RNA_IP[,-grep('seq2019',colnames(model_RNA_IP))]

    model_RNA_I <- model_RNA_I
    model_RNA_I <- model_RNA_I[,ncol(model_RNA_I):1]



  }else if(exp_design_name %in% c("Liver_all_T13")){

    target_epi$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epi$BioCond))))
    target_epiIP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epiIP$BioCond))))
    target_epi_fold$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epi_fold$BioCond))))
    target_RNA_I$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_RNA_I$BioCond))))
    target_RNA_IP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_RNA_IP$BioCond))))

    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond+target_epi$dataset+target_epi$Seq+target_epi$Library)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("Time","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 +target_epiIP$TypeData:target_epiIP$BioCond+target_epiIP$dataset+target_epiIP$Seq+target_epiIP$Library)#
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("Time","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))
    model_epi_bis <- model_epi_bis[,ncol(model_epi_bis):1]
    model_epi <- model_epi_bis

    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond+target_epi_fold$dataset+target_epi_fold$Seq+target_epi_fold$Library)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))

    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond+target_RNA_I$dataset+target_RNA_I$Seq+target_RNA_I$Library)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond+target_RNA_IP$dataset+target_RNA_IP$Seq+target_RNA_IP$Library)#
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))


  }else if(exp_design_name %in% c("Liver_2019")) {
    target_epiIP$RNAprep = as.factor(paste0('prep_',target_epiIP$RNAprep))
    target_epi$RNAprep = as.factor(paste0('prep',target_epi$RNAprep))
    target_epi_fold$RNAprep = as.factor(paste0('prep',target_epi_fold$RNAprep))
    target_epiIP$LibraryIPbis = as.factor(paste0('lib_',target_epiIP$LibraryIPbis))
    target_epi$LibraryIPbis = as.factor(paste0('lib_',target_epi$LibraryIPbis))
    target_epi_fold$LibraryIPbis = as.factor(paste0('lib_',target_epi_fold$LibraryIPbis))
    target_RNA_I$RNAprep = as.factor(paste0('prep_',target_RNA_I$RNAprep))
    target_RNA_I$mRNAprep = as.factor(paste0('prep',target_RNA_I$mRNAprep))
    target_RNA_I$LibraryIPbis = as.factor(paste0('prep',target_RNA_I$LibraryIPbis))
    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$Time:target_epi$BioCond +  target_epi$LibraryIPbis)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("Time","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 +target_epiIP$TypeData:target_epiIP$Time:target_epiIP$BioCond + target_epiIP$RNAprep)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("Time","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))

    model_epi <- model_epi_bis[,ncol(model_epi_bis):1]

    model_epi_fold <- model.matrix(~0+target_epi_fold$Time:target_epi_fold$BioCond + target_epi_fold$RNAprep)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))
    model_epi_fold <- model_epi_fold[,ncol(model_epi_fold):1]
    model_RNA_I <- model.matrix(~0+target_RNA_I$Time:target_RNA_I$BioCond + target_RNA_I$RNAprep)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_I <- model_RNA_I[,ncol(model_RNA_I):1]
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$Time:target_RNA_IP$BioCond)
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))



  }else if(exp_design_name %in% c("Liver_2019_T13")){

    target_epiIP$LibraryIPbis = as.factor(paste0('lib_',target_epiIP$LibraryIPbis))
    target_epi$LibraryIPbis = as.factor(paste0('lib_',target_epi$LibraryIPbis))
    target_RNA_I$RNAprep = as.factor(paste0('lib_',target_RNA_I$RNAprep))
    target_RNA_IP$RNAprep = as.factor(paste0('lib_',target_RNA_IP$RNAprep))
    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond + target_epi$LibraryIPbis)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 + target_epiIP$TypeData:target_epiIP$BioCond + target_epiIP$LibraryIPbis)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels) #[grep('IP',as.character(target_epiIP$labels))]
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))
    model_epi_bis <- model_epi_bis[,ncol(model_epi_bis):1]
    model_epi <- model_epi_bis

    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond+target_epi_fold$LibraryIPbis)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))

    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond + target_RNA_I$RNAprep)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond  + target_RNA_IP$RNAprep)
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))



  }else if(exp_design_name %in% c("LiverZT")){
    target_epi$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epi$BioCond))))
    target_epiIP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epiIP$BioCond))))
    target_epi_fold$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_epi_fold$BioCond))))
    target_RNA_I$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_RNA_I$BioCond))))
    target_RNA_IP$BioCond <- as.factor(gsub('_T.*','',gsub('Liver_', '', as.character(target_RNA_IP$BioCond))))

    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$Time:target_epi$BioCond +  target_epi$Seq)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("Time","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 +target_epiIP$TypeData:target_epiIP$Time:target_epiIP$BioCond + target_epiIP$Seq)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("Time","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))

    model_epi <- model_epi_bis[,ncol(model_epi_bis):1]

    model_epi_fold <- model.matrix(~0+target_epi_fold$Time:target_epi_fold$BioCond + target_epi_fold$Seq)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))
    model_epi_fold <- model_epi_fold[,ncol(model_epi_fold):1]
    model_RNA_I <- model.matrix(~0+target_RNA_I$Time:target_RNA_I$BioCond + target_RNA_I$Seq)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_I <- model_RNA_I[,ncol(model_RNA_I):1]
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$Time:target_RNA_IP$BioCond)
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))


  }else if(exp_design_name %in% c("Liver_2018_T13")){

    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond + target_epi$Seq + target_epi$Library)
    colnames(model_epi_all) = gsub("target_epi\\$","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("TypeData","",colnames(model_epi_all))
    colnames(model_epi_all) = gsub("BioCond","",colnames(model_epi_all))

    model_epi_bis <- model.matrix(~0 +target_epiIP$TypeData:target_epiIP$BioCond + target_epiIP$Seq + target_epiIP$Library)
    rownames(model_epi_bis ) <- as.character(target_epiIP$labels)
    model_epi_bis <- model_epi_bis[,!grepl("nput", colnames(model_epi_bis))]
    colnames(model_epi_bis) = gsub("target_epiIP\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("target_epi\\$","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("TypeData","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("BioCond","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub("IP:","",colnames(model_epi_bis))
    colnames(model_epi_bis) = gsub(":",".",colnames(model_epi_bis))
    model_epi <- model_epi_bis[,ncol(model_epi_bis):1]

    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond + target_epi_fold$Seq + target_epi_fold$Library)
    colnames(model_epi_fold) = gsub(":",".",colnames(model_epi_fold))
    model_epi_fold <- model_epi_fold[,ncol(model_epi_fold):1]
    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond + target_RNA_I$Seq + target_RNA_I$Library)
    colnames(model_RNA_I) = gsub(":",".",colnames(model_RNA_I))
    model_RNA_I <- model_RNA_I[,ncol(model_RNA_I):1]
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond + target_RNA_IP$Seq + target_RNA_IP$Library)
    colnames(model_RNA_IP) = gsub(":",".",colnames(model_RNA_IP))


  }else if(exp_design_name=="CecAm" ){
    model_epi_all <- model.matrix(~0+target_epi$Mice+target_epi$TypeData+target_epi$TypeData:target_epi$BioCond+target_epi$Seq+target_epi$Phase)
    model_epi <- model.matrix(~0+target_epiIP$BioCond + target_epiIP$Seq + target_epiIP$Phase)
    colnames(model_epi) <- gsub('.*\\$|BioCond', '', colnames(model_epi) )
    model_epi_fold <- model.matrix(~0+target_epi_fold$BioCond + target_epi_fold$Seq +target_epi_fold$Phase)
    model_RNA_I <- model.matrix(~0+target_RNA_I$BioCond + target_RNA_I$Seq + target_RNA_I$Phase)
    model_RNA_IP <- model.matrix(~0+target_RNA_IP$BioCond + target_RNA_IP$Seq + target_RNA_IP$Phase)
  }else if(exp_design_name=="LivOld" | exp_design_name == "LiverZT" | exp_design_name == "LiverZTall"){
    model_epi <- model.matrix(~1+target_epi$BioCond + target_epi$Seq)
    model_epi_fold <- model.matrix(~1+target_epi_fold$BioCond + target_epi_fold$Seq )
    model_RNA_I <- model.matrix(~1+target_RNA_I$BioCond + target_RNA_I$Seq )
    model_RNA_IP <- model.matrix(~1+target_RNA_IP$BioCond + target_RNA_IP$Seq )
  }else{
    model_epi <- model.matrix(~1+target_epi$BioCond)
    model_epi_fold <- model.matrix(~1+target_epi_fold$BioCond)
    model_RNA_I <- model.matrix(~1+target_RNA_I$BioCond )
    model_RNA_IP <- model.matrix(~1+target_RNA_IP$BioCond)
  }
  print("Model MeRIPSeq")
  print(model_epi)
  print("Model MeRIPSeq fold")
  print(model_epi_fold)
  print("Model RNASeq")
  print(model_RNA_I)


  colnames(model_epi) = gsub("target_epi\\$Time","",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$BioCond","",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$Seq","",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$Library","Library",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$Phase","",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$RNAprep","",colnames(model_epi))
  colnames(model_epi) = gsub("target_epi\\$dataset","",colnames(model_epi))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$Time","",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$BioCond","",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$Seq","",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$Library","Library",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$Phase","",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$RNAprep","",colnames(model_epi_fold))
  colnames(model_epi_fold) = gsub("target_epi_fold\\$dataset","",colnames(model_epi_fold))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$Time","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$BioCond","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$Seq","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$Library","Library",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$Phase","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$RNAprep","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$mRNAprep","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("&","",colnames(model_RNA_I))
  colnames(model_RNA_I) = gsub("target_RNA_I\\$dataset","",colnames(model_RNA_I))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$Time","",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$BioCond","",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$Seq","",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$Library","Library",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$Phase","",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$RNAprep","",colnames(model_RNA_IP))
  colnames(model_RNA_IP) = gsub("target_RNA_IP\\$dataset","",colnames(model_RNA_IP))

  assign("model_epi_all",model_epi_all,pos = globalenv())
  assign("model_epi",model_epi,pos = globalenv())
  assign("model_epi_fold",model_epi_fold,pos = globalenv())
  assign("model_RNA_I",model_RNA_I,pos = globalenv())
  assign("model_RNA_IP",model_RNA_IP,pos = globalenv())

  target_epiIP <- target_epi[grep('IP',rownames(target_epi)),]
  assign("target_epi",target_epi,pos = globalenv())
  assign("target_epiIP",target_epi,pos = globalenv())
  assign("target_epi_fold",target_epi_fold,pos = globalenv())
  assign("target_RNA_I",target_RNA_I,pos = globalenv())
  assign("target_RNA_IP",target_RNA_IP,pos = globalenv())
  assign("condRef",condRef,pos = globalenv())

}
