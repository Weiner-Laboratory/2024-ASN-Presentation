# Load libraries ----
library(clusterProfiler)
library(Seurat)#v4.3.1
library(dplyr)
library(msigdbr)
library(harmony)
library(Azimuth)
library(openxlsx)
library(officer)
library(qs)
# load idw functions ----
source("/blue/weineid/version-controlled-software/functions/Make experiment folder.R") # nolint
source("/blue/weineid/version-controlled-software/functions/get input.R")
source("/blue/weineid/version-controlled-software/functions/word documentation.R") # nolint
source("/blue/weineid/version-controlled-software/functions/output message.R")
source("/blue/weineid/version-controlled-software/functions/convert cell name.R")
# Set experiment title ----
experiment.title <- get.input("Enter experiment title")
output.data.name <- get.input("Enter name for output dataset")
experiment.title.original <- experiment.title
experiment.subtitle <- Sys.Date()
folder.date.experiment <- make.experiment.folder(experiment.title)
r.script.file <- rstudioapi::getSourceEditorContext()$path
# set variables ----
min.cells <- 5 # 3 in original KH analysis
nFeature.min <- 500 # number of genes in a cell; 500 in original KH analysis;
nFeature.max <- 9999999 # From Susztak paper - 3000
nCountRNA.min <- 2000 # 2000 in original KH analysis; 200 in Susztak paper
nCountRNA.max <- 9999999 # 3000 in Susztak paper
max.mt <- 2 # 2 in original KH analysis; 50 in Susztak paper
p.val.cutoff <- 0.01
log2FC.Cutoff <- 1.0
minimum.pct <- 0.25
de.test <- "wilcox" # default setting is wilcox
StartTime <- Sys.time()
enter.key <- "\n" #blank line to format output
data.directory <- "/blue/weineid/datasets/AR-KO/"
# Import data from idw computer ----
message("Importing data ...")
sample100 <- Read10X(data.dir = paste0(data.directory, "100_filtered_feature_bc_matrix/"))
sample102 <- Read10X(data.dir = paste0(data.directory, "102_filtered_feature_bc_matrix/"))
#create seurat objects for each sample
all.samples <- c("sample100", "sample102")
for (data.name in all.samples) {
  message(paste("Converting", data.name, "..."))
  temp <- eval(as.name(data.name))
  temp <- CreateSeuratObject(
    counts = eval(as.name(data.name)),
    project = data.name,
    min.cells = min.cells,
    min.features = nFeature.min
  )
  temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern = "^mt-")
  temp <- subset(temp, subset = (percent.mt < max.mt))
  assign(data.name, temp)
}
#add genotype and sex to each file
sample100$groupid <- "wt"
sample102$groupid <- "wt"
sample100$sex <- "male"
sample102$sex <- "female"

#merge to create one object
message("Merging files ...")
all.data <- merge(x = sample100, y = list(sample102))

all.data

#log normalize, PCA etc
message("log normalize, etc ...")
all.data <- NormalizeData(all.data,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
all.data <- FindVariableFeatures(all.data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(all.data), 10)
all.genes <- rownames(all.data)
all.data <- ScaleData(all.data, features = all.genes)
all.data <- RunPCA(all.data, features = VariableFeatures(object = all.data))
DimHeatmap(all.data,
           dims = 1:15,
           cells = 500,
           balanced = TRUE)
DimHeatmap(all.data,
           dims = 1:2,
           cells = 5000,
           balanced = TRUE)

FeatureScatter(all.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- VariableFeaturePlot(all.data)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)
plot1 + plot2
VizDimLoadings(all.data, dims = 1:2, reduction = "pca")

message("Find neighbors, clusters, etc ...")
all.data <- FindNeighbors(all.data, dims = 1:20)
all.data <- FindClusters(all.data, resolution = 1.0)
all.data <- RunUMAP(all.data, dims = 1:20)
DimPlot(all.data, pt.size = 0.1, reduction = 'umap')
DimPlot(all.data,
        pt.size = 0.1,
        reduction = 'umap',
        split.by = "sex")
DimPlot(all.data,
        pt.size = 0.1,
        reduction = 'umap',
        group.by = "groupid")

# Harmony batch correction ----
message("Harmony batch correction ...")
all.data <- RunHarmony(
  all.data,
  group.by.vars = 'sex',
  # reduction = 'PCA',
  assay.use = 'RNA',
  project.dim = FALSE
)
all.data <- RunUMAP(all.data, dims = 1:30, # dims = 2:30, original from KH
                    reduction = 'harmony')
all.data <- FindNeighbors(all.data, dims = 1:20)
all.data <- FindClusters(all.data, resolution = 1.0)

p2 <- DimPlot(all.data, group.by = 'groupid', pt.size = 0.1) +
  ggplot2::ggtitle("Harmony integration")
p2
p3 <- DimPlot(all.data, label = TRUE, pt.size = 0.1)
p3

DimPlot(all.data,
        label = FALSE,
        pt.size = 0.1,
        group.by = "groupid")
DimPlot(all.data,
        label = FALSE,
        pt.size = 0.1,
        group.by = "sex")
DimPlot(all.data, label = FALSE, pt.size = 0.1)
# Join layers
all.data <- JoinLayers(all.data)
# run azimuth to identify cell types----
all.data.azimuth <- RunAzimuth(all.data, "kidneyref")
all.data$azimuth.l2 <- all.data.azimuth$predicted.annotation.l2
all.data$azimuth.l2.score <- all.data.azimuth$predicted.annotation.l2.score
all.data$azimuth.l3 <- all.data.azimuth$predicted.annotation.l3
all.data$azimuth.l3.score <- all.data.azimuth$predicted.annotation.l3.score
# save RDS file ----
save.file.name <- paste0(folder.date.experiment, output.data.name)
saveRDS(all.data,
        file = paste0(save.file.name, ".rds")
        )
qsave(all.data, paste0(save.file.name, ".qs"))
all.data.meta.data <- all.data@meta.data
# find celltype-specif degs ----
Idents(all.data) <- all.data$azimuth.l3
try(rm(deg.workbook))
deg.workbook <- createWorkbook()
# get number of cells of each cell type
cell.type.number <- as.data.frame(table(all.data$azimuth.l3))
addWorksheet(deg.workbook, sheetName = "cells")
writeData(deg.workbook, sheet = "cells", cell.type.number)
# get list of cell types with at least min.cells in the group
cell.type.number <- as.data.frame(table(all.data$azimuth.l3))
colnames(cell.type.number) <- c("cell", "N")
cell.type.number.minimum <- cell.type.number %>% filter(N > min.cells)
all.cell.types <- cell.type.number.minimum[, 1]
column.order <- c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")
# go through all.cell.types and find cell-type specific deg's
counter <- 1
message("Find celltype-specific deg's")
for (cell in all.cell.types) {
  message(paste(
    "  ",
    counter,
    "of",
    length(all.cell.types),
    format(Sys.time(), '%H:%M:%S'),
    cell
  ))
  cell.markers <- FindMarkers(
    all.data,
    ident.1 = cell,
    logfc.threshold = log2FC.Cutoff,
    verbose = T,
    min.cells.feature = min.cells,
    min.cells.group = min.cells,
    min.pct = minimum.pct
  )
  cell.markers <- subset(cell.markers, subset = (p_val_adj < p.val.cutoff))
  counter <- counter + 1
  cell.markers$gene <- rownames(cell.markers)
  cell.markers <- cell.markers[, column.order]
  # modify cell name use as tab in Excel sheet
  cell.name <- convert.cell.name(cell)
  sheet.name <- substr(cell.name, 1, 29)
  # add degs to the spreadsheet
  addWorksheet(deg.workbook, sheetName = sheet.name)
  writeData(deg.workbook, sheet = sheet.name, cell.markers)
}
# save excel documents ----
excel.file.name <- paste0(folder.date.experiment, "celltype-specific degs for ",experiment.title, ".xlsx")
saveWorkbook(deg.workbook, file = excel.file.name, overwrite = T)
# document experiment's variables
variables.used <- c(
  "min.cells",
  "minimum.pct",
  "nCountRNA.max",
  "nCountRNA.min",
  "nFeature.min",
  "nFeature.max",
  "max.mt",
  "p.val.cutoff",
  "excel.file.name",
  "r.script.file",
  "all.samples",
  "output.data.name"
)
try(rm(documentation.Word))
documentation.Word <- read_docx()
documentation.Word <-
  body_add_par(documentation.Word,
               paste("Experiment:", experiment.title),
               style = "heading 1")
documentation.Word <-
  body_add_par(documentation.Word,
               paste(
                 "Analysis date:",
                 experiment.subtitle,
                 format(Sys.time(), '%H:%M:%S')
               ),
               style = "heading 2")
documentation.Word <-
  body_add_par(documentation.Word, "")

document.variables (documentation.Word, variables.used)
print(documentation.Word,
      target = paste0(folder.date.experiment, "Experiment summary.docx"))
message("Done")
