# let user enter their choice.
# Use second option to limit choices.
# Use third option to enter control-key options.  Default is TRUE (show Ctrl-x messages)
# Option 4 to input a number

get.input <-
  function(text.line,
           options = "",
           message.control.key = TRUE,
           numeric.only = FALSE) {
    message("") # enter blank line in the console
    # make sure text.line does not have '-> ' at the end of the line
    message(text.line)
    if (message.control.key) {
      message("Press Ctrl-2 to enter your response")
    }
    if (numeric.only) {
      choice <- select.list(1:99)
    } else {
      if (length(options) == 1) {
        # accept any response that is not a blank line
        choice <- ""
        while (nchar(choice) == 0) {
          choice <- readline("Any input is accepted -> ")
        }
      } else {
        # make sure choice is one of the options
        choice <- select.list(options)
      }
    }
    if (message.control.key) {
      message("Press Ctrl-1 to continue processing")
    }
    return(choice)
  }