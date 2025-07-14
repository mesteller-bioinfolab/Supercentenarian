library(zellkonverter)
library(Seurat)

dataset <- "M116"
step <-  "04.1.create_input_for_CellTypist"

projectRoot <- "~/Documents/PROJECTS/07.2_Supercentenaria/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "03_norm_clust_annot", "data.RDS")
data <- readRDS(data_path)

# Convert to h5ad
# NOTE: celltypist need log normalized data
sce_obj <- as.SingleCellExperiment(data, assay = c("RNA"))
zellkonverter::writeH5AD(sce_obj, file.path(saveProcessedData, "logcounts.h5ad"), X_name = 'logcounts')
