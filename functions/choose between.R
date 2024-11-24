# let user choose between options
# enter the code 'source("/blue/weineid/functions/choose between.R")' to run this in primary code

choose.between <- function (options, text.line) {
  # make sure text.line has '-> ' at the end of the line
  if (!substr(text.line, nchar(text.line) - 3, nchar(text.line))
      == "-> ") {
    text.line <- paste(text.line, "-> ")
  }
  if (length(options) == 1) {
    # accept any response that is not a blank line
    choice <- ""
    message("Press Ctrl-2 to enter your response")
    while (nchar(choice) == 0) {
      choice <- readline(text.line)
    }
  } else {
    choice <- ""
    message("Press Ctrl-2 to enter your response")
    while (!choice %in% options) {
      choice <- readline(paste(text.line, "-> "))
      choice <- tolower(choice)
    }
  }
  message("Press Ctrl-1 to continue processing")
  return(choice)
}
get.input <- function(text.line, options = "") {
  # make sure text.line has '-> ' at the end of the line
  if (!substr(text.line, nchar(text.line) - 3, nchar(text.line))
      == "-> ") {
    text.line <- paste(text.line, "-> ")
  }
  if (length(options) == 1) {
    # accept any response that is not a blank line
    choice <- ""
    message("Press Ctrl-2 to enter your response")
    while (nchar(choice) == 0) {
      choice <- readline(text.line)
    }
  } else {
    choice <- ""
    message("Press Ctrl-2 to enter your response")
    while (!choice %in% options) {
      choice <- readline(text.line)
      choice <- tolower(choice)
    }
  }
  message("Press Ctrl-1 to continue processing")
  return(choice)
}
