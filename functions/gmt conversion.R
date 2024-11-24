convert.gmt <-
  function(gmt.to.convert) {
    my.test <- NULL
    my.test <- tibble("gs_name", "gene_symbol")
    colnames(my.test) <- c("gs_name", "gene_symbols")
    line.counter <- 1
    # for (name.counter in 1:length(gmt.to.convert[["geneset.names"]])) {
    for (name.counter in 1:25) {
      message(gmt.to.convert[["geneset.names"]][[name.counter]])
      for (genesets.counter in 1:length(gmt.to.convert[["genesets"]][[name.counter]])) {
        my.test[line.counter, 1] <- gmt.to.convert[["geneset.names"]][name.counter]
        my.test[line.counter, 2] <- gmt.to.convert[["genesets"]][[name.counter]][genesets.counter]
        line.counter <- line.counter + 1
      }
    }
    return(my.test)
  }
convert.entrez.gmt <- function(gmt.to.convert){
  temp <- gmt.to.convert[["genesets"]]
  names(temp) <- gmt.to.convert[["geneset.names"]]
  my.list <- temp
  return(my.list)
}