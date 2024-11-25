# load libraries ----
library(Seurat)
library(officer)
library(qs) # used for faster seurat data save and read access
library(harmony)
library(ggplot2)
library(ensembldb)
library(biomaRt)
# load idw functions ----
source("/blue/weineid/version-controlled-software/functions/Make experiment folder.R") # nolint
source("/blue/weineid/version-controlled-software/functions/get input.R")
source("/blue/weineid/version-controlled-software/functions/word documentation.R") # nolint
source("/blue/weineid/version-controlled-software/functions/data set information.R")
source("/blue/weineid/version-controlled-software/functions/convert gene and ensembl.R")
# get experiment details ----
experiment.title <- get.input("Enter experiment title")
experiment.title.original <- experiment.title
experiment.subtitle <- Sys.Date()
folder.date.experiment <- make.experiment.folder(experiment.title)
r.script.file <- rstudioapi::getSourceEditorContext()$path
data.info <- get.data.info()
time.start <- Sys.time()
# start Word documentation ----
documentation.Word <- NULL # nolint
documentation.Word <- read_docx() # nolint
documentation.Word <- body_add_par(documentation.Word,
                                   # nolint
                                   paste("Experiment:", experiment.title),
                                   style = "heading 1")
documentation.Word <- body_add_par(documentation.Word,
                                   # nolint
                                   paste("Analysis date:", experiment.subtitle),
                                   style = "heading 2")
documentation.Word <- body_add_par(documentation.Word, "Analysis parameters", style = "heading 2") # nolint
# load data set ----
choices <- rownames(data.info)
data.set.choice <- get.input ("Which data set to use?", choices)
# create data file name and read the data and subset for normal if human data
data.file.name <- paste0("/blue/weineid/datasets/", data.info[data.set.choice, "filename"])
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
if (data.info[data.set.choice, "disease.subset"] == "yes") {
  message("  Subsetting for normal data")
  data.set <- subset(data.set, subset = (disease == "normal"))
}
# set 'idents'
if (data.info[data.set.choice, "cell.ID"] == "subclass.l2") {
  Idents(data.set) <- data.set$subclass.l2
}
if (data.info[data.set.choice, "cell.ID"] == "cell_type") {
  Idents(data.set) <- data.set$cell_type
}
# normalize data ----
message("log normalize, etc ...")
data.set <- NormalizeData(data.set,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
data.set <- FindVariableFeatures(data.set, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data.set), 10)
all.genes <- rownames(data.set)
data.set <- ScaleData(data.set, features = all.genes)
data.set <- RunPCA(data.set, features = VariableFeatures(object = data.set))
DimHeatmap(data.set,
           dims = 1:15,
           cells = 500,
           balanced = TRUE)
DimHeatmap(data.set,
           dims = 1:2,
           cells = 5000,
           balanced = TRUE)
FeatureScatter(data.set, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- VariableFeaturePlot(data.set)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)
plot1 + plot2
VizDimLoadings(data.set, dims = 1:2, reduction = "pca")
# JackStraw analysis ----
message("JackStraw analysis ...")
data.set <- JackStraw(data.set, num.replicate = 5) # originally 100
data.set <- ScoreJackStraw(data.set, dims = 1:20)
JackStrawPlot(data.set, dims = 1:20)
ElbowPlot(data.set)
# get number of dimensions for further reduction ----
dimensions <- get.input("Enter number of dimensions to use", , , TRUE)
message("Find neighbors, clusters, etc ...")
data.set <- FindNeighbors(data.set, dims = 1:dimensions)
data.set <- FindClusters(data.set, resolution = 1.0)
data.set <- RunUMAP(data.set, dims = 1:dimensions)

DimPlot(data.set,
        pt.size = 0.1,
        raster = F,
        reduction = 'umap')
# Harmony batch correction ----
message("Harmony batch correction ...")
data.set <- RunHarmony(
  data.set,
  group.by.vars = data.info[data.set.choice, "group.by"],
  # reduction = 'PCA',
  assay.use = 'RNA',
  project.dim = FALSE
)
data.set <- RunUMAP(data.set, dims = 1:dimensions, # dims = 2:30, original from KH
                    reduction = "harmony")
data.set <- FindNeighbors(data.set, dims = 1:dimensions)
data.set <- FindClusters(data.set, resolution = 1.0)
p2 <- DimPlot(data.set,
              group.by = data.info[data.set.choice, "group.by"],
              pt.size = 0.2,
              raster = F) +
  ggplot2::ggtitle("Harmony integration")
p2
p3 <- DimPlot(data.set,
              label = TRUE,
              pt.size = 0.2,
              raster = F)
p3

DimPlot(
  data.set,
  label = FALSE,
  pt.size = 0.2,
  raster = F,
  group.by = data.info[data.set.choice, "group.by"]
)
DimPlot(data.set, label = FALSE, pt.size = 0.1)

# Join layers
# data.set <- JoinLayers(data.set)
# save the normalized and other adjusted dataset
path.normalized <- "/blue/weineid/datasets/normalized/"
filename.normalized <- paste0(path.normalized, data.info[data.set.choice, "filename"])
qsave (data.set, filename.normalized)
# convert top10 to gene if in ensembl code ----
if (data.info[data.set.choice, "gene.type"] == "ensembl") {
  ensembl.dataset <- data.info[data.set.choice, "ensembl.dataset"]
  ensembl <- useEnsembl(biomart = "genes", dataset = ensembl.dataset)
  top10.table <- convert.ensembl.to.gene(top10)
  top10.gene <- top10.table[1,2]
  for (counter in 2:10) {
    top10.gene <- paste(top10.gene, top10.table[counter,2], sep = ", ")
  }
}
time.end <- Sys.time()
# Document the work done in Word ----
variables.used <- c(
  "data.file.name",
  "data.set.choice",
  "filename.normalized",
  "r.script.file",
  "top10",
  "dimensions",
  "time.start",
  "time.end"
)
if (data.info[data.set.choice, "gene.type"] == "ensembl") {
  variables.used <- c(variables.used,
                      "top10.gene")
}
document.variables (documentation.Word, variables.used)
print(documentation.Word,
      target = paste0(folder.date.experiment, experiment.title, ".docx"))

# done
message("Done!")
