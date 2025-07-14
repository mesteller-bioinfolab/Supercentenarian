library(Seurat)
library(tidyverse)  
library(writexl)
library(scDblFinder)
library(DoubletFinder)
library(BiocParallel)

dataset <- "M116"
step <-  "03_norm_clust_annot"

projectRoot <- "~/Documents/PROJECTS/07.2_Supercentenaria/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "02_doublets", "data.RDS")
data <- readRDS(data_path)

# Scale
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)

# Linear dimensional reduction
data <- RunPCA(data, features = VariableFeatures(object = data))
pca <- DimPlot(data, reduction = "pca", group.by = "sample")
ggsave(plot = pca, filename = here::here("plots", analysis, "pca.pdf"), width = 8, height = 6)

ElbowPlot(data)
ggsave(file = here::here("plots", analysis, "elbow_plot.pdf"), width = 8, height = 6)

# Clustering
data <- FindNeighbors(data, dims = 1:15)

# Non-linear dimensional - visualization
data <- RunUMAP(data, reduction.use = "pca", dims = 1:15)

data <- FindClusters(data, resolution = 0.5)


plot <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = T)
ggsave(plot, filename = file.path(savePlots, "RNA_snn_res.0.5.pdf"), width = 9, height = 7)

saveRDS(data, file.path(saveProcessedData, "data.RDS"))

Idents(data) <- "RNA_snn_res.0.5"
markers <- FindAllMarkers(data, only.pos = T)
saveRDS(markers, file.path(saveResults, "RNA_snn_res.0.5.RDS"))
writexl::write_xlsx(markers, file.path(saveResults, "RNA_snn_res.0.5.xlsx"))

