library(tidyverse)
library(Seurat)

dataset <- "M116"
step <-  "04.3.add_CellTyspist_pred"

projectRoot <- "~/Documents/PROJECTS/07.2_Supercentenaria/"
source(file.path(projectRoot, "misc", "colors.R"))

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "03_norm_clust_annot", "data.RDS")
data <- readRDS(data_path)


# Read predictions 
celltypist_low <- read.csv(file.path(projectRoot, "results", dataset, "04.2.run_celltypist", "celltypist_pred_Immune_All_Low.csv"))
celltypist_high <- read.csv(file.path(projectRoot, "results", dataset, "04.2.run_celltypist", "celltypist_pred_Immune_All_High.csv"))

# Clean format
rownames(celltypist_low) <- celltypist_low$X
rownames(celltypist_high) <- celltypist_high$X
celltypist_low$X <- NULL
celltypist_high$X <- NULL
colnames(celltypist_low) <- paste0(colnames(celltypist_low), "_Immune_All_Low")
colnames(celltypist_high) <- paste0(colnames(celltypist_high), "_Immune_All_High")

# Incorporate them to Seurat object
data <- data %>% AddMetaData(celltypist_low)
data <- data %>% AddMetaData(celltypist_high)

# Replace "/" by " or "
data$majority_voting_Immune_All_High <- gsub("/", " or ", data$majority_voting_Immune_All_High)
data$predicted_labels_Immune_All_High <- gsub("/", " or ", data$predicted_labels_Immune_All_High)
data$majority_voting_Immune_All_Low <- gsub("/", " or ", data$majority_voting_Immune_All_Low)
data$predicted_labels_Immune_All_Low <- gsub("/", " or ", data$predicted_labels_Immune_All_Low)

# Save object
saveRDS(data, file = file.path(saveProcessedData, "data.RDS"))

data$majority_voting_Immune_All_Low <- factor(data$majority_voting_Immune_All_Low, levels = names(cols_majority_voting_Immune_All_Low))
data$majority_voting_Immune_All_High <- factor(data$majority_voting_Immune_All_High, levels = names(cols_majority_voting_Immune_All_High))

plot <- DimPlot(data, group.by = "majority_voting_Immune_All_High", label = F, cols = cols_majority_voting_Immune_All_High)
ggsave(plot, filename = file.path(savePlots, "majority_voting_Immune_All_High.pdf"), width = 9, height = 7)

plot <- DimPlot(data, group.by = "majority_voting_Immune_All_High", label = T, cols = cols_majority_voting_Immune_All_High)
ggsave(plot, filename = file.path(savePlots, "majority_voting_Immune_All_High_labeled.pdf"), width = 9, height = 7)

plot <- DimPlot(data, group.by = "majority_voting_Immune_All_Low", label = F, cols = cols_majority_voting_Immune_All_Low)
ggsave(plot, filename = file.path(savePlots, "majority_voting_Immune_All_Low.pdf"), width = 9, height = 7)

plot <- DimPlot(data, group.by = "majority_voting_Immune_All_Low", label = T, cols = cols_majority_voting_Immune_All_Low)
ggsave(plot, filename = file.path(savePlots, "majority_voting_Immune_All_Low_labeled.pdf"), width = 9, height = 7)


# Differencial expression analysis
Idents(data) <- "majority_voting_Immune_All_High"
markers_Immune_All_High <- FindAllMarkers(data, only.pos = T)
saveRDS(markers_Immune_All_High, file.path(saveResults, "DEGs_majority_voting_Immune_All_High.RDS"))

markers_split <- list()
markers_split <- markers_Immune_All_High %>%
  split(.$cluster) %>%
  map(~ .x %>% filter(p_val_adj<0.01) %>% 
        arrange(p_val_adj, desc(avg_log2FC)))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, file.path(saveResults, "DEGs_majority_voting_Immune_All_High.xlsx"))


Idents(data) <- "majority_voting_Immune_All_Low"
markers_Immune_All_Low <- FindAllMarkers(data, only.pos = T)
saveRDS(markers_Immune_All_Low, file.path(saveResults, "DEGs_majority_voting_Immune_All_Low.RDS"))

markers_split <- list()
markers_split <- markers_Immune_All_Low %>%
  split(.$cluster) %>%
  map(~ .x %>% filter(p_val_adj<0.01) %>% 
        arrange(p_val_adj, desc(avg_log2FC)))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, file.path(saveResults, "DEGs_majority_voting_Immune_All_Low.xlsx"))


# Heatmaps
markers_Immune_All_High %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
plot <- DoHeatmap(data, features = top10$gene, group.by = "majority_voting_Immune_All_High", group.colors = cols_majority_voting_Immune_All_High) + NoLegend()
ggsave(plot, filename = file.path(savePlots, "heatmap_majority_voting_Immune_All_High.png"), width = 12, height = 9)

markers_Immune_All_Low %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
plot <- DoHeatmap(data, features = top10$gene, group.by = "majority_voting_Immune_All_Low", group.colors = cols_majority_voting_Immune_All_Low) + NoLegend()
ggsave(plot, filename = file.path(savePlots, "heatmap_majority_voting_Immune_All_Low.png"), width = 15, height = 17)


