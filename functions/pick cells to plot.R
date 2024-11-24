pick.cells.to.plot <- function (data.set, how.many = 0) {
  source("/blue/weineid/version-controlled-software/functions/get input.R")
  # set up variables
  all.cells.unique <- unique(Idents(data.set))
  cells.to.plot <- NULL # start here for new set of cells
  get.cells <- TRUE
  # message user
  message("Select cells to plot, 0 for all cells or to exit")
  message("<Ctrl-2> to enter choices")
  while (get.cells) {
    new.cell <- select.list(all.cells.unique, multiple = T)
    if (length(new.cell) > 0) {
      for (counter in 1:length(new.cell)) {
        cells.to.plot <- c(cells.to.plot, as.character(new.cell[counter]))
        message(new.cell[counter], " added to list")
        if (how.many == 1) {
          get.cells <- FALSE
        }
      }
    } else {
      get.cells <- FALSE
    }
  }
  if (length(cells.to.plot) == 0) {
    cells.to.plot <- all.cells.unique
  } else {
    if (length(cells.to.plot) > 1) {
      cells.to.plot <- cells.to.plot[1]
    }
  }
  message("<Ctrl-1> to continue step-by-step")
  return(cells.to.plot)
}
