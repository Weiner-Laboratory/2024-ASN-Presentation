# identify gene to select cells ----
gene.select <- ""
ensembl.name <- F
while (!ensembl.name) {
  gene.select <- get.input("Gene to use to select cells")
  ensembl.name <- as.character(convert.to.ensembl(gene.select))
}
gene.cutoff <- 0
gene.cutoff <- get.input("Gene expression cutoff")
experiment.title <- paste(experiment.title.original,
                          gene.select,
                          gene.cutoff)
# enter gene.select as subset criteria below
if (gene.select == "all") {
  data.set <- data.set.original
} else {
  data.set <- subset(data.set.original, 
                     subset = (Car2 > gene.cutoff))
}
message(ncol(data.set), 
        " cells selected from ", 
        ncol(data.set.original), 
        " in original data set.")
