library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "M116"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/07.2_Supercentenaria/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset)

# Load data
expression_matrix <- Read10X(data.dir = here::here(raw_data_path, "filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = expression_matrix, min.cells = 3)

data
# An object of class Seurat 
# 23394 features across 8953 samples within 1 assay 
# Active assay: RNA (23394 features, 0 variable features)
# 1 layer present: counts

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    36    1925    2334    2436    2817    7590 


# QC METRICS
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# Visualize QC metrics before filtering
# number of cells before filtering: 384
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_before_filtering.pdf"), width = 8, height = 8)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot = plot1 + plot2, filename = file.path(savePlots, "QC_scatter_plot.pdf"), width = 12, height = 6)

pdf(file = file.path(savePlots, "n_feature_density.pdf"))
plot(density(data$nFeature_RNA))
abline(v = 300)
dev.off()

print(summary(data$nFeature_RNA))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  36    1925    2334    2436    2817    7590 



# Visualize the number UMIs/transcripts per cell
plot <- data@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)
ggsave(plot = plot, filename = file.path(savePlots, "n_feature_density.pdf"), width = 12, height = 6)

data <- subset(data, subset = nFeature_RNA > 300 & percent.mt < 20)

plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_filtering.pdf"), width = 8, height = 8)

saveRDS(data, file.path(saveProcessedData, "data.RDS"))

