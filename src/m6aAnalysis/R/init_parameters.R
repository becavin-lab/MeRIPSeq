init_parameters <- function(){

  #print(project_dir)
  if (!dir.exists(project_dir)){
    dir.create(project_dir)
  }
  setwd(project_dir)

  figure_folder <- paste(project_dir,"/FigSummary/",sep = "")
  if (!dir.exists(figure_folder)){
    dir.create(figure_folder)
  }
  assign("figure_folder",figure_folder, pos = globalenv())

  motif_folder <- paste(project_dir,"/Motif/",sep = "")
  if (!dir.exists(motif_folder)){
    dir.create(motif_folder)
  }
  assign("motif_folder",motif_folder, pos = globalenv())

  path_guitar <- paste(project_dir,"/GUITARPlot/",sep="")
  if (!dir.exists(path_guitar)){
    dir.create(path_guitar)
  }
  folderBed <- "GUITARPlot/bed/"
  if (!dir.exists(folderBed)){
    dir.create(folderBed)
  }
  assign("path_guitar",path_guitar, pos = globalenv())



}
