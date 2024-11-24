# make folder with first level being the date and the second level being the experiment.title.  
# Default is 'experiment.title'

make.experiment.folder <- function (title = experiment.title) {
  folder.date <- paste0("/blue/weineid/output/", Sys.Date(),"/")
  folder.date.experiment <- paste0(folder.date, title,"/")
  
  if (!file.exists(folder.date)) {
    dir.create(folder.date)
  }
  if (!file.exists(folder.date.experiment)) {
    dir.create(folder.date.experiment)
  }
  message("The folder, '",folder.date.experiment, "', is now present")
  return(folder.date.experiment)
}
