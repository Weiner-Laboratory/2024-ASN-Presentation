# convert to % value, with indicated number of digits to right of '.'
percent <- function(value, digits = 1) {
  value <- paste0(round(100*value,digits),"%")
  return(value)
}