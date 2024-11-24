# print variables used in a script to Word document
# first variable is the officer variable for the word document
# second variable is the list of variables to document
# function code ----
document.variables <- function(word.document, variables.to.document) {
  word.document <- body_add_par(word.document, "Variables used in analysis", style = "heading 2")
  for (variable in variables.to.document) {
    variable.data <- eval(as.name(variable))
    if (is.vector(variable.data)) {
      variable.string <- ""
      for (counter in 1:length(variable.data)) {
        variable.string <- paste(variable.string, variable.data[counter])
      }
      variable.string <- substr(variable.string, 2, nchar(variable.string))
      variable.data <- variable.string
    }
    word.document <- body_add_par(word.document, paste0(variable, ": ", variable.data))
  }
  # return(word.document)
  return("variables added to document")
}
