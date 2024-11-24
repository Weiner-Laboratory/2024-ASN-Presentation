# save plot as pdf file ----
save.plot <- function(plot,
                      file.name = plot.filename,
                      folder = folder.date.experiment,
                      orientation = "landscape") {
  library(ggplot2)
  if (!(substr(folder, nchar(folder), nchar(folder)) == "/")) {
    folder <- paste0(folder, "/")
  }
  file.with.folder <- paste0(folder, file.name)
  if (orientation == "landscape") {
    height.value <- 7.5
    width.value <- 10
  } else {if (tolower(orientation) == "portrait") {
      height.value <- 10
      width.value <- 7.5
    } else {
        message ("Can't interpret orientation option.  File not saved")
        return ("error")
      }
  }
  ggsave(
    plot = plot,
    file = file.with.folder,
    height = height.value,
    width = width.value,
    device = "pdf"
  )
  message("   saved to pdf file")
}
# create filename for a plot ----
plot.file <- function(plot.description, generic.name = title.with.genes) {
  plot.file <- paste0(generic.name," ", plot.description, ".pdf")
  return(plot.file)
}
# save plot to word document using officer ----
plot.to.word <- function(word.document, plot, size = "large") {
  size <- tolower(size)
  dimension <- NULL
  dimension <- matrix(nrow = 3, ncol = 4)
  dimension[1, ] <- c("large", 10, 6.5, 0.5)
  dimension[2, ] <- c("small", 3, 6.5, 0.75)
  dimension[3, ] <- c("medium", 4.5, 6.5, 1.0)
  dimension <- as.data.frame(dimension)
  rownames(dimension) <- dimension[, 1]
  colnames(dimension) <- c("size", "height", "width", "scaling")
  dimension$height <- as.numeric(dimension$height)
  dimension$width <- as.numeric(dimension$width)
  dimension$scaling <- as.numeric(dimension$scaling)
  
  word.document <- body_add_gg(
    word.document,
    plot,
    width = dimension[size, "width"],
    height = dimension[size, "height"],
    scale = dimension[size, "scaling"],
    pos = "after"
  )
  return(word.document)
}
