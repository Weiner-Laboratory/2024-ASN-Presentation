get.dataset <- function(kidney.choice) {
  choices <- rownames(data.info)
  kidney.choice <- get.input ("Which data set to use?", choices)
  ensembl.dataset <- data.info[kidney.choice, "ensembl.dataset"]
  # create data file name and read the data and subset for normal if human data
  data.file.name <- paste0("/blue/weineid/datasets/", data.info[kidney.choice, "filename"])
  if (if.error(nchar(data.file.name.old),
               TRUE,!(data.file.name == data.file.name.old))) {
    message("Reading data")
    file.extension <- substr(data.file.name,
                             nchar(data.file.name) - 2,
                             nchar(data.file.name))
    if (file.extension == ".qs") {
      data.set <- qread(data.file.name) #uses qread for qs formatted files
    }
    if (file.extension == "rds") {
      data.set <- readRDS(data.file.name)
    }
    data.set.original <- data.set
  } else {
    message("Using same data file...")
    data.set <- data.set.original
  }
  Sys.sleep(delay.message)
  disease.subset <- data.info[kidney.choice, "disease.subset"]
  if (disease.subset == "yes") {
    choices <- unique(data.set$disease)
    kidney.subset <- get.input("Which subset to use?", choices)
    message("Subsetting for disease = ", kidney.subset)
    data.set <- subset(data.set, subset = (disease == kidney.subset))
  }
  # set 'idents'
  if (data.info[kidney.choice, "cell.ID"] == "azimuth.l2") {
    data.set$idw_celltype <- data.set$azimuth.l2
  }
  if (data.info[kidney.choice, "cell.ID"] == "azimuth.l3") {
    data.set$idw_celltype <- data.set$azimuth.l3
  }
  if (data.info[kidney.choice, "cell.ID"] == "cell_type") {
    data.set$idw_celltype <- data.set$cell_type
  }
  Idents(data.set) <- data.set$idw_celltype
  # load BiomaRT data ----
  gene.type <- data.info[kidney.choice, "gene.type"]
  get.ensembl <- F
  if (gene.type == "ensembl") {
    if (!is.error(ensembl)) {
      if (!ensembl@dataset == ensembl.dataset) {
        get.ensembl <- T
      }
    } else {
      get.ensembl <- T
    }
    if (get.ensembl) {
      ensembl <- useEnsembl(biomart = "genes", dataset = ensembl.dataset)
    }
  }
 return(data.set) 
}