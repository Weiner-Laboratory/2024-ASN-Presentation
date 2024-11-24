# receive a list of cell types and let user pick one of them
pick.a.cell <- function(cell.name.list) {
  cell.names <- unique(cell.name.list)
  pick.a.cell <- select.list(cell.names)
  pick.a.cell <- as.character(pick.a.cell)
#  pick.a.cell <- cell.names[pick.a.cell]
  return(pick.a.cell)
}


pick.a.cell.old <- function(cell.name.list) {
  cell.names <- unique(cell.name.list)
  for (counter in 1:length(cell.names)) {
    cell.details <- paste0(counter, ".  ", cell.names[counter])
    message(cell.details)
  }
  cell.chosen.number <- get.input("Enter 'all' for all cells, cell # to choose a cell, or 'exit' when done adding cells)")
  if(cell.chosen.number == "exit"){return("done")}
  cell.chosen.number <- as.numeric(cell.chosen.number)
  if (is.na(cell.chosen.number)) {
    message("Enter cell number, not name")
    return("error")
  }
  else {
    if (cell.chosen.number <= length(cell.names)) {
      cell.chosen <- as.character(cell.names[cell.chosen.number])
      return(cell.chosen)
    } else {
      message("Number needs to be in range of 1-", length(cell.names))
      return ("error")
    }
  }
}
