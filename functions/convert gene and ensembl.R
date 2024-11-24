convert.to.ensembl <- function (gene.name) {
  ensembl.ID <- getBM(
    filters = "external_gene_name",
    values = gene.name,
    attributes = "ensembl_gene_id",
    mart = ensembl
  )
  library(ensembldb)
  library(biomaRt)
  if (nrow(ensembl.ID) == 0) {
    message("Cannot match this gene in ensembl")
    documentation <- paste0("Gene '", gene.next, "' could not be matched to Ensembl")
    documentation.Word <- body_add_par(documentation.Word, documentation)
    return("")
  } else {
    ensembl.ID.length <- nchar(ensembl.ID)
    for (counter in (1:nrow(ensembl.ID))) {
      ensembl.ID.internal <- ensembl.ID[counter, 1]
      message(paste("Ensembl ID for", gene.name, "is", ensembl.ID.internal))
      if (is.element(ensembl.ID.internal, all.genes)) {
        return (ensembl.ID.internal)
      } else {
        message(paste(ensembl.ID.internal, "not in data set"))
        return("")
      }
    }
  }
}
convert.ensembl.to.gene <- function (ensembl.ID){ # takes a list and returns ensembl and gene names
  gene.name <- getBM(
    mart = ensembl,
    values = ensembl.ID,
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id","external_gene_name"),
    uniqueRows = T
  )
  return(gene.name)
}
# convert gene to entrez-id ----
convert.gene.to.entrez <- function (gene.list){
  gene.name <- getBM(
    mart = ensembl,
    values = gene.list,
    filters = "external_gene_name",
    attributes = c("ensembl_gene_id","external_gene_name", "entrezgene_id"),
    uniqueRows = T
  )
}