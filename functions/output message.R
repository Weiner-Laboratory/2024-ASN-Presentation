# Output 'the.text' to console and to word.document
output.message <- function(the.text, 
                           word.document = documentation.Word, 
                           message.delay = 2
                           ) {
  message(the.text)
  Sys.sleep(message.delay)
  word.document <- body_add_par(word.document,
                                the.text,
                                pos = "before")
  return(word.document)
}
# source("/blue/weineid/version-controlled-software/functions/output message.R")
