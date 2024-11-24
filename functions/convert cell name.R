# convert cell name to useable name
convert.cell.name <-
  function(cell.name) {
    cell.name <- gsub("/", "-", cell.name) #remove a '/' that is sometimes present
    cell.name <- gsub("Proximal Tubule Epithelial Segment ", "PT-S", cell.name) # shorten name
    cell.name <- gsub("Distal Convoluted Tubule Type", "DCT", cell.name)
    cell.name <- gsub("Outer Medullary Collecting Duct", "OMCD", cell.name)
    cell.name <- gsub("Cortical Collecting Duct", "CCD", cell.name)
    cell.name <- gsub("Thick Ascending Limb", "TAL", cell.name)
    
    loop.ongoing <- TRUE
    while (loop.ongoing) {
      cell.name.length <- nchar(cell.name)
      if (substr(cell.name, cell.name.length, cell.name.length) == " ") {
        cell.name <- substr(cell.name, 1, cell.name.length - 1)
      } else {
        loop.ongoing <- FALSE
      }
    }
    return(cell.name)
  }
remove.dashes <- function(name) {
  new.name <- gsub("-", " ", name)
  new.name <- gsub("_", " ", new.name)
  return(new.name)
}