get.gene.name <- function() {
  gene.name.to.return <- NULL
  gene.to.plot <- NULL
  while (is.null(gene.name.to.return)) {
    gene.next <- get.input("Enter gene name", "", F) # nolint
    if (data.info[kidney.choice, "gene.case"] == "upper") {
      gene.next <- toupper(gene.next)
    }
    if (data.info[kidney.choice, "gene.case"] == "title") {
      gene.next <- str_to_title(gene.next)
    }
    if (data.info[kidney.choice, "gene.type"] == "gene") {
      if (is.element(gene.next, all.genes)) {
        gene.to.plot <- c(gene.to.plot, gene.next)
        gene.name.to.return <- c(gene.name.to.return, gene.next)
      } else {
        documentation <- paste(gene.next, "is not in data set")
        message(documentation)
        documentation.Word <- body_add_par(documentation.Word, documentation) # nolint
      }
    } else {
      ensembl.ID <- getBM(
        # nolint
        filters = "external_gene_name",
        values = gene.next,
        attributes = "ensembl_gene_id",
        mart = ensembl
      )
      if (nrow(ensembl.ID) == 0) {
        message("Cannot match this gene in ensembl")
        documentation <- paste0("Gene '", gene.next, "' could not be matched to Ensembl") # nolint
        documentation.Word <- body_add_par(documentation.Word, documentation) # nolint
      } else {
        ensembl.ID.length <- nchar(ensembl.ID) # nolint
        for (counter in (1:nrow(ensembl.ID))) {
          ensembl.ID.internal <- ensembl.ID[counter, 1]
          message(paste(
            "Ensembl ID for",
            gene.next,
            "is",
            ensembl.ID.internal
          ))
          if (is.element(ensembl.ID.internal, all.genes)) {
            gene.to.plot <- c(gene.to.plot, ensembl.ID.internal)
            gene.name.to.return <- c(gene.name.to.return, gene.next)
            documentation <- NULL
            documentation <- paste0("Gene '",
                                    gene.next,
                                    ".' Ensembl is ",
                                    ensembl.ID.internal,
                                    ".")
            documentation.Word <- body_add_par(documentation.Word, documentation, style = "Normal")
          } else {
            message(paste(ensembl.ID.internal, "not in data set"))
            documentation <- paste0(
              "Gene '",
              gene.next,
              "', Ensembl ID, ",
              ensembl.ID.internal,
              ", is not in this data set"
            )
            documentation.Word <- body_add_par(documentation.Word, documentation, style = "Normal")
          }
        }
      }
    }
  }
  return(gene.name.to.return)
}