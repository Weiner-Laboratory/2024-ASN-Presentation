# receive a list of cell types and let user pick one of them
pick.a.cell <- function(cell.name.list) {
  source("/blue/weineid/functions/get input.R")
  cell.names <- as.list(cell.name.list)
  cell.names <- unique(cell.names)
  counter <- 1
  for (cell in cell.names){
    cell.details <- paste0(counter, ".  ", cell)
    message(cell.details)
    counter <- counter + 1
  }
  cell.chosen.number <- get.input("Enter cell number to choose")
  cell.chosen.number <- as.numeric(cell.chosen.number)
  cell.chosen <- cell.names[[cell.chosen.number]]
  message("The cell chosen was ", cell.chosen)
  return(cell.chosen)
}