message.with.delay <- function(documentation, delay.time = 2) {
  message(documentation)
  Sys.sleep(delay.time)
}
