#Necessary packages for SPATIAL
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
#install.packages("hdf5r")
library(hdf5r)
library(SeuratDisk)


#-Load M116 data------------------------------------------------------------------------------
M116_data <- readRDS("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/processed_data/M116/04.3.add_CellTyspist_pred/data.RDS")

DimPlot(M116_data, reduction = "umap", group.by = "majority_voting_Immune_All_Low", label = TRUE)

#-PLOTS
FeaturePlot(M116_data, features = "ACSM3", split.by = "group") & theme(plot.title = element_text(size=10))
FeaturePlot(M116_data, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

VlnPlot(M116_data, features = "MCM6", pt.size = 0.1, group.by = "predicted.celltype.l2", split.by = "group", cols = c("palegreen1", "steelblue1", "tan2", "indianred1")) + NoLegend()

count_matrix <- as.matrix(GetAssayData(M116_data, assay = "RNA", slot = "counts"))
write.csv(count_matrix, file = "M116_data_raw_counts_matrix.csv", quote = FALSE)

# Load predicted labels
predicted_labels <- read.csv("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/Celltypist_ALL_dataset/M116/predicted_labels.csv", row.names = 1)

# Check the structure
head(predicted_labels)

all(rownames(predicted_labels) %in% colnames(M116_data))

# Add predicted labels to metadata
M116_data <- AddMetaData(M116_data, metadata = predicted_labels)

M116_data$celltypist <- predicted_labels$cell_type

DimPlot(M116_data, group.by = "predicted_labels", label = TRUE)

# Define cell types to remove for too little cells
cells_to_remove <- c("CD16- NK cells", "Mast cells", "Mono-mac", "B cells", "DC", "Pro-B cells", "Transitional B cells", "NK cells", "CRTAM+ gamma-delta T cells", "Epithelial cells", "Cycling T cells", "Monocytes", "Follicular helper T cells", "Double-negative thymocytes", "Plasmablasts", "Alveolar macrophages", "MNP", "Double-positive thymocytes", "DC1", "Macrophages", "Type 17 helper T cells", "CD8a/a", "ILC3", "Germinal center B cells", "Megakaryocyte precursor", "Migratory DCs", "Fibroblasts", "CD8a/b(entry)", "Treg(diff)", "Tem/Effector helper T cells PD1+", "Type 1 helper T cells")

# Subset the object to KEEP only desired cells
M116_filtered <- subset(M116_data, subset = predicted_labels %in% cells_to_remove == FALSE)

DimPlot(M116_filtered, group.by = "predicted_labels", label = TRUE)

#-Subsetting B cells for ABC specific markers--------------------------------------------------
# Method 1: Using subset()
Bcells_M116_data <- subset(
  M116_data,
  subset = majority_voting_Immune_All_Low %in% c(
    "Age-associated B cells",
    "Naive B cells",
    "Memory B cells"
  )
)

DimPlot(Bcells_M116_data, reduction = "umap", group.by = "majority_voting_Immune_All_Low", label = TRUE)

# Extract the current metadata
metadata <- Bcells_M116_data@meta.data

# Rename "Naive B cells" and "Memory B cells" to "Other B cells"
metadata$majority_voting_Immune_All_Low <- ifelse(
  metadata$majority_voting_Immune_All_Low %in% c("Naive B cells", "Memory B cells"),
  "Other B cells",
  metadata$majority_voting_Immune_All_Low
)

# Update the Seurat object's metadata
Bcells_M116_data@meta.data <- metadata

# Extract metadata
metadata <- Bcells_M116_data@meta.data

# Rename "1" to "Other B cells" (assuming "1" = Naive/Memory B cells)
metadata$majority_voting_Immune_All_Low <- ifelse(
  metadata$majority_voting_Immune_All_Low == "1",
  "Age-associated B cells",
  metadata$majority_voting_Immune_All_Low
)

# Update the Seurat object
Bcells_M116_data@meta.data <- metadata

table(Bcells_M116_data@meta.data$majority_voting_Immune_All_Low)

DimPlot(Bcells_M116_data, reduction = "umap", group.by = "majority_voting_Immune_All_Low", label = TRUE)

#-Markers 2 cell types-----------------------------
# Using 'majority_voting_Immune_All_Low' as identities:
Idents(Bcells_M116_data) <- "majority_voting_Immune_All_Low"
table(Idents(Bcells_M116_data))  

# Find markers for all clusters (using Wilcoxon rank-sum test by default)
markers <- FindAllMarkers(
  Bcells_M116_data,
  min.pct = 0.4,           
  logfc.threshold = 0.25,  
  only.pos = TRUE,         
  assay = "RNA"            
)

#-Save as Excel file
write.xlsx(
  markers, 
  file = "B cells only markers min.pct0.4.xlsx", 
  rowNames = FALSE,  # Don't write row names since we have gene column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")  # Make header bold
)

#-GO enrichment ABCs with B cells isolated--------------------------------------------------------

head(Genes_GOenrich_ABCcluster_top100)

library(clusterProfiler)
library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse
library(enrichplot)

org_db <- org.Hs.eg.db

# Initialize list to store results
pathway_results <- list()

# Loop through each cell type
for (ct in names(Genes_GOenrich_ABCcluster_top100)) {
  gene_list <- Genes_GOenrich_ABCcluster_top100[[ct]]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org_db)
  
  # Run KEGG or GO enrichment
  enrich_res <- enrichGO(gene         = entrez_ids$ENTREZID,
                         OrgDb        = org_db,
                         keyType      = "ENTREZID",
                         ont          = "BP",     # Biological Process
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  # Store result
  pathway_results[[ct]] <- enrich_res
}

# View top pathways for a cell type
head(pathway_results[["gdT"]])

# Dotplot
dotplot(pathway_results[["CD16 Mono"]], showCategory = 10) + ggtitle("CD14 Mono Pathway Enrichment")

# Barplot
barplot(enrich_res, showCategory = 10)


#-PATHWAYS for UPreg and Downreg genes in ABC-------------------------------------------------------------------

head(Genes_GOenrich_ABCcluster)

library(clusterProfiler)
library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse
library(enrichplot)

org_db <- org.Hs.eg.db

# Initialize list to store results
pathway_results <- list()

# Loop through each cell type
for (ct in names(Genes_GOenrich_ABCcluster)) {
  gene_list <- Genes_GOenrich_ABCcluster[[ct]]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org_db)
  
  # Run KEGG or GO enrichment
  enrich_res <- enrichGO(gene         = entrez_ids$ENTREZID,
                         OrgDb        = org_db,
                         keyType      = "ENTREZID",
                         ont          = "MF",     # Biological Process
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  # Store result
  pathway_results[[ct]] <- enrich_res
}

# View top pathways for a cell type
head(pathway_results[["gdT"]])

# Dotplot
dotplot(pathway_results[["CD16 Mono"]], showCategory = 10) + ggtitle("CD14 Mono Pathway Enrichment")

# Barplot
barplot(enrich_res, showCategory = 10)

#-KEGG new way for GO enrich-------------------------------------------------------------------------------------

library(dplyr)

# Keep significant upregulated genes
genes_up <- Table_SX_KEGG %>%
  dplyr::filter(FC > 0 & p_val_adj < 0.05) %>%
  pull(gene)

library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db if you're working with mouse data

entrez_up <- bitr(genes_up, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# KEGG enrichment analysis
kegg_res <- enrichKEGG(gene = entrez_up$ENTREZID,
                       organism = 'hsa',  # 'mmu' for mouse
                       pvalueCutoff = 0.05)

# View top pathways
head(kegg_res)

library(ReactomePA)

reactome_res <- enrichPathway(gene = entrez_up$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)

# View enriched pathways
head(reactome_res)

# KEGG
barplot(kegg_res, showCategory = 10, title = "KEGG Pathways")

# Reactome
barplot(reactome_res, showCategory = 10, title = "Reactome Pathways")

# Filter descriptions containing your terms
terms_of_interest <- c("immune", "cytokine", "aging", "lifespan", "senescence", "autophagy")

# Get matching descriptions
filtered_terms <- as.data.frame(reactome_res) %>%
  dplyr::filter(grepl(paste(terms_of_interest, collapse = "|"), Description, ignore.case = TRUE))

# Plot only filtered pathways manually (custom plot)
library(ggplot2)
ggplot(filtered_terms, aes(x = reorder(Description, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Filtered Reactome Pathways",
       x = NULL,
       y = "Gene Count") +
  theme_minimal()

ggplot(filtered_terms, aes(x = reorder(Description, -Count), y = Count, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log10(adj p-value)") +
  labs(title = "Filtered Reactome Pathways",
       x = NULL,
       y = "Gene Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14))

#-For GSEA--------------------------------------------------------------------------------
# Named vector of logFCs
gene_list <- Table_SX_KEGG$avg_log2FC
names(gene_list) <- Table_SX_KEGG$gene

# Convert gene names to Entrez
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list) <- gene_list_entrez$ENTREZID

# GSEA using KEGG
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = 'hsa',
                     nPerm = 1000,
                     minGSSize = 10,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

dotplot(gsea_kegg, showCategory = 10)

#-Load Hashimoto data (other supercentenarians SCs)--------------------------------------------------------------------
data_Hashimoto <- readRDS("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/07.2_Supercentenarian/processed_data/Hashimoto_et_al/04.3.add_CellTyspist_pred/data.RDS")

DimPlot(data_Hashimoto, reduction = "umap", group.by = "predicted_labels", label = TRUE)

# Rename the column in the Seurat object's metadata
colnames(data_Hashimoto@meta.data)[colnames(data_Hashimoto@meta.data) == "majority_voting_Immune_All_Low"] <- "predicted_labels"

# Replace "or" with "/" in the predicted_labels column
data_Hashimoto@meta.data$predicted_labels <- gsub(
  pattern = " or ", 
  replacement = "/", 
  x = data_Hashimoto@meta.data$predicted_labels
)

DimPlot(data_Hashimoto, group.by = "predicted_labels", label = TRUE)

# Define cell types to remove (e.g., "Fibroblasts" and "Doublets")
cells_to_remove <- c("CD16- NK cells", "Mast cells", "Mono-mac", "B cells", "DC", "Pro-B cells", "Transitional B cells", "NK cells", "CRTAM+ gamma-delta T cells", "Epithelial cells", "Cycling T cells", "Monocytes", "Follicular helper T cells", "Double-negative thymocytes", "Plasmablasts", "Alveolar macrophages", "MNP", "Double-positive thymocytes", "DC1", "Macrophages", "Type 17 helper T cells", "CD8a/a", "ILC3", "Germinal center B cells", "Megakaryocyte precursor", "Migratory DCs", "Fibroblasts", "CD8a/b(entry)", "Treg(diff)", "Tem/Effector helper T cells PD1+", "Type 1 helper T cells")

# Subset the object to KEEP only desired cells
data_Hashimoto_filtered <- subset(data_Hashimoto, subset = predicted_labels %in% cells_to_remove == FALSE)

DimPlot(data_Hashimoto_filtered, group.by = "predicted_labels", label = TRUE)

#-PLOTS
FeaturePlot(data_Hashimoto, features = "ACSM3", split.by = "group") & theme(plot.title = element_text(size=10))
FeaturePlot(data_Hashimoto, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

VlnPlot(data_Hashimoto, features = "MCM6", pt.size = 0.1, group.by = "predicted.celltype.l2", split.by = "group", cols = c("palegreen1", "steelblue1", "tan2", "indianred1")) + NoLegend()

#-Load Terakhova data (youngers)----------------------------------------------------------------

# load h5ad as Seurat
Terekhova_data = schard::h5ad2seurat('all_pbmcs_rna_subsampled_10000_by_ag.h5ad')

count_matrix <- as.matrix(GetAssayData(Terekhova_data, assay = "RNA", slot = "counts"))
write.csv(count_matrix, file = "Terekhova_data_raw_counts_matrix.csv", quote = FALSE)

# Load predicted labels
predicted_labels <- read.csv("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/Celltypist_ALL_dataset/Terekhova/predicted_labels.csv", row.names = 1)

# Check the structure
head(predicted_labels)

all(rownames(predicted_labels) %in% colnames(Terekhova_data))

# Add predicted labels to metadata
Terekhova_data <- AddMetaData(Terekhova_data, metadata = predicted_labels)

Terekhova_data$celltypist <- predicted_labels$cell_type

DimPlot(Terekhova_data, group.by = "predicted_labels", reduction = "umap", label = TRUE)

# Define cell types to remove (e.g., "Fibroblasts" and "Doublets")
cells_to_remove <- c("CD16- NK cells", "Mast cells", "Mono-mac", "B cells", "DC", "Pro-B cells", "Transitional B cells", "NK cells", "CRTAM+ gamma-delta T cells", "Epithelial cells", "Cycling T cells", "Monocytes", "Follicular helper T cells", "Double-negative thymocytes", "Plasmablasts", "Alveolar macrophages", "MNP", "Double-positive thymocytes", "DC1", "Macrophages", "Type 17 helper T cells", "CD8a/a", "ILC3", "Germinal center B cells", "Megakaryocyte precursor", "Migratory DCs", "Fibroblasts", "CD8a/b(entry)", "Treg(diff)", "Tem/Effector helper T cells PD1+", "Type 1 helper T cells", "HSC/MPP", "NKT cells", "Intermediate macrophages", "gamma-delta T cells", "ILC", "Proliferative germinal center B cells", "ETP", "Transitional NK", "Early lymphoid/T lymphoid", "Neutrophil-myeloid progenitor", "MEMP", "Monocyte precursor", "GMP", "Kupffer cells", "DC precursor", "Cycling NK cells", "DC3", "Megakaryocyte-erythroid-mast cell progenitor", "Small pre-B cells", "T(agonist)", "Early erythroid", "Memory CD4+ cytotoxic T cells", "Late erythroid")

# Subset the object to KEEP only desired cells
Terekhova_data_filtered <- subset(Terekhova_data, subset = predicted_labels %in% cells_to_remove == FALSE)

DimPlot(Terekhova_data_filtered, group.by = "predicted_labels", label = TRUE)

#---------------------------------------------------------------------------------------------------------------
#-Integration of SCs and Younger datasets-----------------------------------------------------------------------

# Rename the column directly in the metadata
names(Terekhova_data_filtered@meta.data)[names(Terekhova_data_filtered@meta.data) == "Age_group"] <- "sample"

# Verify the change
head(Terekhova_data_filtered@meta.data$sample)  # Should show your sample data
colnames(Terekhova_data_filtered@meta.data)     # Should show "sample" instead of "Age_group"

# Set Age_range for each object
Terekhova_data_filtered$Age_general <- "Young"
data_Hashimoto$Age_general <- "SCs"
M116_filtered$Age_general <- "SCs"

# Get intersection of features across all objects
common_features <- Reduce(intersect, list(
  rownames(M116_filtered),
  rownames(data_Hashimoto),
  rownames(Terekhova_data_filtered)
))

# Subset all to common features
M116_filtered <- subset(M116_filtered, features = common_features)
data_Hashimoto <- subset(data_Hashimoto, features = common_features)
Terekhova_data_filtered <- subset(Terekhova_data_filtered, features = common_features)

Terekhova_data_filtered[["RNA"]] <- Terekhova_data_filtered[["originalexp"]]
DefaultAssay(Terekhova_data_filtered) <- "RNA"
Terekhova_data_filtered[["originalexp"]] <- NULL

Terekhova_data_filtered <- NormalizeData(Terekhova_data_filtered)
Terekhova_data_filtered <- FindVariableFeatures(Terekhova_data_filtered, selection.method = "dispersion")

# Add metadata to track sample origin
M116_filtered$dataset <- "M116"
data_Hashimoto$dataset <- "Hashimoto"
Terekhova_data_filtered$dataset <- "Terekhova"

# Create list
obj.list <- list(M116_filtered, data_Hashimoto, Terekhova_data_filtered)

# Normalize and find variable features (for all objects)
obj.list <- lapply(obj.list, NormalizeData)
obj.list <- lapply(obj.list, FindVariableFeatures)

Terekhova_data <- subset(Terekhova_data, cells = sample(colnames(Terekhova_data), 20000))

Terekhova_data <- FindVariableFeatures(Terekhova_data, selection.method = "dispersion")

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list)
Young_SCs_integrated_Celltype <- IntegrateData(anchorset = anchors)

# Run standard workflow on integrated data
Young_SCs_integrated_Celltype <- ScaleData(Young_SCs_integrated_Celltype)
Young_SCs_integrated_Celltype <- RunPCA(Young_SCs_integrated_Celltype)
Young_SCs_integrated_Celltype <- RunUMAP(Young_SCs_integrated_Celltype, dims = 1:30)
# Cluster cells (adjust resolution here!)
Young_SCs_integrated_Celltype <- FindNeighbors(Young_SCs_integrated_Celltype, dims = 1:30) # Use PC range from ElbowPlot
Young_SCs_integrated_Celltype <- FindClusters(
  Young_SCs_integrated_Celltype,
  resolution = 0.8)

DimPlot(Young_SCs_integrated_Celltype, reduction = "umap", group.by = "predicted_labels", split.by = "dataset", label = TRUE)

Young_SCs_integrated_Celltype <- subset(
  Young_SCs_integrated_Celltype,
  subset = !is.na(Young_SCs_integrated_Celltype$predicted_labels)
)

#-RDS objects
saveRDS(Young_SCs_integrated_Celltype, file="Young_SCs_integrated_Celltype.rds")

Young_SCs_integrated_CellType$Age_range <- as.factor(ordered(Young_SCs_integrated$Age_range, levels=c("Young","SCs")))
Young_SCs_integrated$dataset <- as.factor(ordered(Young_SCs_integrated$dataset, levels=c("Terekhova","Hashimoto","M116")))

#-PLOTS
DefaultAssay(Young_SCs_integrated) <- "RNA"
FeaturePlot(Young_SCs_integrated, features = "S100A9", split.by = "Age_range") & theme(plot.title = element_text(size=10))
FeaturePlot(Young_SCs_integrated, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

VlnPlot(Young_SCs_integrated, features = "MCM6", pt.size = 0.1, group.by = "predicted.celltype.l2", split.by = "group", cols = c("palegreen1", "steelblue1", "tan2", "indianred1")) + NoLegend()

#-DEGs between Young and SCs and Young and M116----------------------------------------------

#Terekhova vs. Hashimoto
degs_tere_vs_hash <- FindMarkers(
  Young_SCs_integrated,
  ident.1 = "Terekhova",
  ident.2 = "Hashimoto",
  group.by = "dataset",       # Compare datasets (not celltypes)
  test.use = "wilcox",        # Default test
  logfc.threshold = 0.25,     # Minimum log2 fold change
  min.pct = 0.1,              # Gene detected in ≥10% of cells
  only.pos = FALSE            # Return all DEGs
)

# Add gene column and comparison label
degs_tere_vs_hash <- degs_tere_vs_hash %>%
  tibble::rownames_to_column("gene") %>%
  mutate(comparison = "Terekhova_vs_Hashimoto")

#M116 vs. Terekhova
degs_m116_vs_tere <- FindMarkers(
  Young_SCs_integrated_Azimuth,
  ident.1 = "M116",
  ident.2 = "Terekhova",
  group.by = "dataset",
  test.use = "MAST",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(comparison = "M116_vs_Terekhova")

# M116 vs. Hashimoto
degs_m116_vs_hashimoto <- FindMarkers(
  object = Young_SCs_integrated_Azimuth,
  ident.1 = "M116",
  ident.2 = "Hashimoto",
  group.by = "dataset",
  test.use = "MAST",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.cells.group = 3,  # Minimum cells per group
  only.pos = FALSE,
  verbose = TRUE
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(comparison = "M116_vs_Hashimoto")  # Fixed comparison label

# Filter out any potential NA results
degs_m116_vs_hashimoto <- degs_m116_vs_hashimoto %>%
  dplyr::filter(!is.na(p_val_adj))

# View top differentially expressed genes
head(degs_m116_vs_hashimoto %>% arrange(p_val_adj))

#Combine results
all_degs <- bind_rows(degs_m116_vs_tere, degs_m116_vs_hashimoto)

library(openxlsx)
# Save as Excel file
write.xlsx(
  all_degs, 
  file = "all_degs_M116vsAll.xlsx", 
  rowNames = FALSE,  # Don't write row names since we have gene column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")  # Make header bold
)

#-DEG between Young and SCs ALL-------------------------------------------------------------

# Duplicate the 'dataset' column to 'dataset2'
Young_SCs_integrated_Celltype@meta.data$dataset2 <- Young_SCs_integrated_Celltype@meta.data$dataset

# Change values in dataset2 column
Young_SCs_integrated_Celltype@meta.data$dataset2 <- ifelse(
  Young_SCs_integrated_Celltype@meta.data$dataset %in% c("Hashimoto", "Terekhova"),
  "NotM116",
  Young_SCs_integrated_Celltype@meta.data$dataset2  # Keep original value if not matched
)

degs_young_vs_scs <- FindMarkers(
  object = Young_SCs_integrated_Celltype,  # Your Seurat object
  ident.1 = "M116",                      # First group (Young)
  ident.2 = "NotM116",                        # Second group (SCs)
  group.by = "dataset2",                 # Metadata column containing the groups
  test.use = "MAST",                    # Wilcoxon rank-sum test (default)
  logfc.threshold = 0.25,                 # Minimum log2 fold change
  min.pct = 0.1,                          # Gene detected in ≥10% of cells
  min.cells.group = 3,  # Minimum cells per group
  only.pos = FALSE                        # Return both up/downregulated genes
) %>%
  tibble::rownames_to_column("gene") %>%  # Convert rownames to a column
  mutate(comparison = "M116_vs_NotM116")      # Label the comparison

library(openxlsx)
write.xlsx(
  degs_young_vs_scs,
  file = "DEGs_M116_vs_NotM166.xlsx",
  rowNames = FALSE,  # Gene names are already in a column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

#By celltype DEGs-------------------------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(openxlsx)

# Get all cell types
cell_types <- unique(Young_SCs_integrated_Celltype$predicted_labels)
all_degs <- list()

for (celltype in cell_types) {
  # Subset to current cell type
  cells_subset <- subset(Young_SCs_integrated_Celltype, 
                         subset = predicted_labels == celltype)
  
  # Skip if too few cells
  if (length(unique(cells_subset$dataset)) < 2) next
  
  cat("\nProcessing cell type:", celltype, "\n")
  
  # Function to create valid sheet names
  create_sheet_name <- function(celltype, comparison) {
    # Shorten cell type name if needed
    short_celltype <- substr(gsub("[^[:alnum:]]", "", celltype), 1, 15)
    short_comparison <- ifelse(comparison == "M116_vs_Terekhova", "vsTere",
                               ifelse(comparison == "M116_vs_Hashimoto", "vsHash", comparison))
    
    # Combine and ensure <= 31 chars
    sheet_name <- paste0(short_celltype, "_", short_comparison)
    substr(sheet_name, 1, 31)
  }
  
  # M116 vs Terekhova
  degs_tere <- tryCatch({
    res <- FindMarkers(
      cells_subset,
      ident.1 = "M116",
      ident.2 = "Terekhova",
      group.by = "dataset",
      test.use = "MAST",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      min.cells.group = 3,
      only.pos = FALSE,
      verbose = FALSE
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(comparison = "M116_vs_Terekhova",
             celltype = celltype)
    
    # Add to results with valid sheet name
    sheet_name <- create_sheet_name(celltype, "M116_vs_Terekhova")
    all_degs[[sheet_name]] <- res
    res
  }, error = function(e) {
    message(paste("Error in", celltype, "vs Terekhova:", e$message))
    NULL
  })
  
  # M116 vs Hashimoto
  degs_hash <- tryCatch({
    res <- FindMarkers(
      cells_subset,
      ident.1 = "M116",
      ident.2 = "Hashimoto",
      group.by = "dataset",
      test.use = "MAST",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      min.cells.group = 3,
      only.pos = FALSE,
      verbose = FALSE
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(comparison = "M116_vs_Hashimoto",
             celltype = celltype)
    
    # Add to results with valid sheet name
    sheet_name <- create_sheet_name(celltype, "M116_vs_Hashimoto")
    all_degs[[sheet_name]] <- res
    res
  }, error = function(e) {
    message(paste("Error in", celltype, "vs Hashimoto:", e$message))
    NULL
  })
}

# Remove any NULL entries
all_degs <- all_degs[!sapply(all_degs, is.null)]

# Save to Excel
if (length(all_degs) > 0) {
  write.xlsx(
    all_degs,
    file = "DEGs_by_celltype_correct.xlsx",
    rowNames = FALSE,
    firstRow = TRUE,
    headerStyle = createStyle(textDecoration = "bold"),
    colWidths = "auto"
  )
  cat("Successfully saved DEG results to DEGs_by_celltype_comparisons.xlsx\n")
} else {
  cat("No DEG results to save.\n")
}

#---Azimuth annotation-------------------------------------------------------------------------------------------

devtools::install_github("satijalab/seurat", "seurat5")
devtools::install_github("satijalab/seurat-data", "seurat5")
devtools::install_github("satijalab/azimuth", "seurat5")
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(Matrix)
library(irlba)

library(Matrix)
options(future.globals.maxSize= 1000000000)
BiocManager::install("TFBSTools", type = "source", force = TRUE)
colSums(Young_SCs_integrated_Azimuth)
TFBSTools::colSums(Young_SCs_integrated_Azimuth)
colSums(Young_SCs_integrated_Azimuth)

DefaultAssay(Young_SCs_integrated_Azimuth) <- "RNA"
Young_SCs_integrated_Azimuth <- JoinLayers(Young_SCs_integrated_Azimuth)
DefaultAssay(Young_SCs_integrated_Azimuth) <- "SCT"
Young_SCs_integrated_Azimuth <- RunAzimuth(Young_SCs_integrated_Azimuth, reference = "pbmcref")

DimPlot(Young_SCs_integrated, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3)
DimPlot(Young_SCs_integrated, group.by = "predicted.celltype.l2", split.by = "group", label = FALSE, label.size = 3, raster=FALSE)
DimPlot(Young_SCs_integrated, group.by = "predicted.celltype.l2", split.by = "group", label = FALSE, label.size = 3, raster=FALSE) + scale_color_manual(values = c("#15A390","#CE5755","#c22f87","#5E3203","#F37612", "#F5ADC5", "#719CB9", "#ad8f57", "#A9D39F", "cornsilk4", "#E3E500", "#C7F7FF", "#3957c4", "#C994F2", "lightcyan3", "cyan4", "darksalmon", "#8A6C09", "#A8D5E2", "darkorange", "#FF715B", "#5398BE", "#F2CD5D", "yellowgreen", "darkolivegreen", "peru", "indianred4", "lightskyblue", "plum4", "bisque2"))

#-RDS objects
saveRDS(Young_SCs_integrated, file="Young_SCs_integrated_Azimuth.rds")

count_matrix <- as.matrix(GetAssayData(Young_SCs_integrated_Azimuth, assay = "RNA", slot = "counts"))
write.csv(count_matrix, file = "Young_SCs_integrated_Azimuth_data_raw_counts_matrix.csv", quote = FALSE)

library(Seurat)
library(dplyr)

# Extract average expression (log-normalized recommended)
avg_expr <- AverageExpression(
  Young_SCs_integrated,
  assays = "RNA",                # Use "integrated" if harmonized
  group.by = "sample",          # Column with your 3 groups
  layer = "data",                # "data" for normalized, "counts" for raw
  return.seurat = FALSE
)$RNA                            # $RNA for RNA assay (or $integrated)

# 2. Convert to dataframe with gene names as a column
avg_expr_df <- as.data.frame(avg_expr) %>%
  tibble::rownames_to_column(var = "Gene")  # Moves gene names from row names to column

library(openxlsx)
# Save as Excel file
write.xlsx(
  avg_expr_df, 
  file = "avg_expr_ALLsamples.xlsx", 
  rowNames = FALSE,  # Don't write row names since we have gene column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")  # Make header bold
)

#-Subset only T cells from annotated ALL datasets-------------------------------------------------------
t_subtypes <- c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "dnT", "gdT", "Treg")
t_cells <- subset(Young_SCs_integrated_Azimuth, 
                  subset = predicted.celltype.l2 %in% t_subtypes)

# Visualize
DimPlot(t_cells, group.by = "predicted.celltype.l2", label = TRUE, cols = colour_order) + NoLegend()
DimPlot(t_cells, group.by = "predicted.celltype.l2", label = TRUE, cols = colour_order, split.by = "Age_range")

t_cells$Age_range <- as.factor(ordered(t_cells$Age_range, levels=c("Young","SCs")))
t_cells$dataset <- as.factor(ordered(t_cells$dataset, levels=c("Terekhova","Hashimoto","M116")))

# Re-run analysis on just T cells -> only if necessary 
#t_cells <- NormalizeData(t_cells)
#t_cells <- FindVariableFeatures(t_cells)
#t_cells <- ScaleData(t_cells)
#t_cells <- RunPCA(t_cells)
#t_cells <- RunUMAP(t_cells, dims = 1:30)

# Find T cell markers
DimPlot(t_cells, group.by = "predicted.celltype.l2", label = TRUE) + NoLegend()

#---------Cell frequency table and number of cells per cluster------------

colour_order <- c("#15A390","#CE5755","#c22f87","#5E3203","#F37612", "#F5ADC5", "#719CB9", "#ad8f57", "#A9D39F", "cornsilk4", "#E3E500", "#C7F7FF")

# Create a frequency table using Azimuth cell types and group
freq_table_sample <- as.data.frame(
  prop.table(
    x = table(t_cells$predicted.celltype.l2, t_cells@meta.data$Age_range),
    margin = 2
  )
)

# Rename columns for clarity
colnames(freq_table_sample) <- c("CellType", "Age_range", "Proportion")

# Plot the data
ggplot(data = freq_table_sample, aes(x = Age_range, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Age_range", y = "Proportion of Cells", fill = "Cell Type") +
  theme_minimal() +
  scale_fill_manual(values=colour_order)

# Create a table of the number of cells per cell type and group
cellsPerCluster <- as.data.frame(
  table(t_cells$predicted.celltype.l2, t_cells@meta.data$Age_range)
)

# Rename columns for clarity
colnames(cellsPerCluster) <- c("CellType", "Age_range", "CellCount")

# View the cell count table
View(cellsPerCluster)

#-Graphs for T cells naïve and senescent percentages------------------------------------------------------------
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(dplyr)

# 1. Calculate percentage of naive cells
plot_data <- t_cells@meta.data %>%
  group_by(dataset, sample) %>%
  summarise(
    naive_pct = 100 * sum(predicted.celltype.l2 %in% c("CD4 Naive", "CD8 Naive")) / n(),
    .groups = 'drop'
  ) %>%
  mutate(dataset = factor(dataset))

# 2. Calculate mean and SEM
summary_data <- plot_data %>%
  group_by(dataset) %>%
  summarise(
    mean_pct = mean(naive_pct),
    sem = sd(naive_pct)/sqrt(n()),
    y_position = mean_pct + sem + 2,
    .groups = 'drop'
  )

# Custom colors
custom_colors <- c("#91BFDB", "gold", "#FE9000")

# 3. Create base plot
p <- ggplot() +
  geom_bar(
    data = summary_data,
    aes(x = dataset, y = mean_pct, fill = dataset),
    stat = "identity", width = 0.6, alpha = 0.7
  ) +
  scale_fill_manual(values = custom_colors) +
  geom_errorbar(
    data = summary_data,
    aes(x = dataset, ymin = mean_pct - sem, ymax = mean_pct + sem),
    width = 0.2
  ) +
  geom_quasirandom(
    data = plot_data,
    aes(x = dataset, y = naive_pct),
    color = "black", size = 3, width = 0.2
  ) +
  geom_text(
    data = summary_data,
    aes(x = dataset, y = y_position, label = sprintf("%.1f%%", mean_pct)),
    size = 4, fontface = "bold"
  ) +
  labs(
    title = "Percentage of Naive T cells (CD4+CD8+)",
    x = NULL, y = "Percentage of total T cells"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

# 4. Add SINGLE p-value between first two groups
if(length(levels(plot_data$dataset)) >= 2) {
  groups <- levels(plot_data$dataset)
  group1 <- groups[1]
  group2 <- groups[2]
  
  # Get data for first two groups
  data1 <- plot_data$naive_pct[plot_data$dataset == group1]
  data2 <- plot_data$naive_pct[plot_data$dataset == group2]
  
  if(length(data1) >= 2 && length(data2) >= 2) {
    # Calculate p-value
    p_val <- t.test(data1, data2)$p.value
    p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
    
    # Calculate position (1.4x max height)
    y_max <- max(summary_data$mean_pct[1:2] + summary_data$sem[1:2])
    y_pos <- y_max * 1.4
    
    # Add to plot - ONLY ONCE
    p <- p + 
      # Horizontal line
      annotate("segment",
               x = 1, xend = 2,
               y = y_pos, yend = y_pos,
               color = "black",
               linewidth = 0.8) +
      
      # Text label (ONLY ONCE)
      annotate("text",
               x = 1.5, y = y_pos * 1.05,  # Slightly above line
               label = p_text,
               size = 5,
               fontface = "bold") +
      
      # Adjust y-axis
      coord_cartesian(ylim = c(0, y_pos * 1.2))
  }
}

# Print the plot
print(p)

#-Senescent cells

# Define markers (replace if necessary)
marker_pos <- "B3GAT1"  # Marker you want to be positive
marker_neg <- "CD28"  # Marker you want to be negative

# Create a new metadata column for the condition
t_cells$pos_pos_neg_neg <- ifelse(
  GetAssayData(t_cells, assay = "RNA")[marker_pos, ] > 0 &  # Positive for marker_pos
    GetAssayData(t_cells, assay = "RNA")[marker_neg, ] == 0,  # Negative for marker_neg
  "Pos_Pos_Neg_Neg",
  "Other"
)

# 1. First ensure your metadata contains the group column
head(t_cells@meta.data[c("dataset", "pos_pos_neg_neg")])

# 2. Create the split plot
split_plot <- DimPlot(t_cells,
                      group.by = "pos_pos_neg_neg",
                      split.by = "dataset",  # Column containing your experimental groups
                      cols = c("gray90", "red"),  # Other cells in gray, target cells in red
                      order = TRUE) +      # Plot positive cells on top
  ggtitle(paste(marker_pos, "+ /", marker_neg, "- cells by Group")) +
  theme(legend.position = "top")

# 3. Adjust layout if needed
split_plot + plot_layout(guides = "collect")

library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)

# Define your custom color palette for bars
custom_colors <- c("#91BFDB", "gold", "#FE9000")

# 1. Calculate percentages for each sample
plot_data <- t_cells@meta.data %>%
  group_by(dataset, sample) %>%
  summarize(
    pos_count = sum(pos_pos_neg_neg == "Pos_Pos_Neg_Neg"),
    total_cells = n(),
    pos_pct = (pos_count / total_cells) * 100
  ) %>%
  ungroup()

# 2. Prepare errorbar data (only for datasets with >1 sample)
errorbar_data <- plot_data %>%
  group_by(dataset) %>%
  filter(n() > 1) %>%  # Only include datasets with multiple samples
  ungroup()

# 3. Create base plot with statistical comparisons
p <- ggplot(plot_data, aes(x = dataset, y = pos_pct)) +
  # Mean bars with custom colors
  stat_summary(
    fun = mean, 
    geom = "bar", 
    aes(fill = dataset),
    width = 0.6, 
    alpha = 0.7
  ) +
  scale_fill_manual(values = custom_colors) +
  # Error bars
  geom_errorbar(
    data = errorbar_data,
    stat = "summary",
    fun.data = mean_se,
    width = 0.2
  ) +
  # Individual samples in black
  geom_quasirandom(
    color = "black",  # Changed to black
    size = 3, 
    width = 0.2, 
    show.legend = FALSE
  ) +
  # Mean percentage labels
  stat_summary(
    aes(label = sprintf("%.1f%%", ..y..)),
    fun = mean,
    geom = "text",
    vjust = -1.5,
    size = 4,
    fontface = "bold"
  ) +
  # Add p-value comparison between first two datasets
  stat_compare_means(
    comparisons = list(c(levels(factor(plot_data$dataset))[1], 
                         levels(factor(plot_data$dataset))[2])),
    method = "t.test",
    label = "p.format",
    tip.length = 0.01,
    vjust = -0.5,
    bracket.size = 0.7,
    label.size = 4,
    y.position = max(plot_data$pos_pct) * 1.1  # Position above highest bar
  ) +
  labs(
    title = paste(marker_pos, "+", marker_neg, "- cells"),
    x = NULL, 
    y = "Percentage of cells"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Remove legend
  ) +
  theme(
    axis.line.x = element_line(linewidth = 1.2, color = "black"),
    axis.line.y = element_line(linewidth = 1.2, color = "black")
  )

# 4. Print the plot
print(p)

# Define your markers (replace with your actual gene/peak names)
marker_pos1 <- "B3GAT1"  # First marker you want to be positive
marker_pos2 <- "KLRG1"   # Second marker you want to be positive
marker_neg <- "CD28"     # Marker you want to be negative

# Create a new metadata column for the condition
t_cells$pos_pos_neg_neg <- ifelse(
  GetAssayData(t_cells, assay = "RNA")[marker_pos1, ] > 0 &  # Positive for first marker
    GetAssayData(t_cells, assay = "RNA")[marker_pos2, ] > 0 &  # Positive for second marker
    GetAssayData(t_cells, assay = "RNA")[marker_neg, ] == 0,   # Negative for marker_neg
  "Pos_Pos_Neg_Neg",
  "Other"
)

# 1. First ensure your metadata contains the group column
head(t_cells@meta.data[c("dataset", "pos_pos_neg_neg")])

# 2. Create the split plot
split_plot <- DimPlot(t_cells,
                      group.by = "pos_pos_neg_neg",
                      split.by = "dataset",  # Column containing your experimental groups
                      cols = c("gray90", "red"),  # Other cells in gray, target cells in red
                      order = TRUE) +      # Plot positive cells on top
  ggtitle(paste(marker_pos1, "+", marker_pos2, "+ /", marker_neg, "- cells by Group")) +
  theme(legend.position = "top")

# 3. Adjust layout if needed
split_plot + plot_layout(guides = "collect")

library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)

# Define markers
marker_pos1 <- "B3GAT1"
marker_pos2 <- "KLRG1"
marker_neg <- "CD28"

# Create cell population identifier
t_cells$target_population <- ifelse(
  GetAssayData(t_cells, assay = "RNA")[marker_pos1, ] > 0 &
    GetAssayData(t_cells, assay = "RNA")[marker_pos2, ] > 0 &
    GetAssayData(t_cells, assay = "RNA")[marker_neg, ] == 0,
  "Double_Pos_Neg",
  "Other"
)

# Define your custom color palette for bars
custom_colors <- c("#91BFDB", "gold", "#FE9000")

# Calculate percentages
plot_data <- t_cells@meta.data %>%
  group_by(dataset, sample) %>%
  summarize(
    pop_pct = mean(target_population == "Double_Pos_Neg") * 100,
    .groups = 'drop'
  ) %>%
  mutate(dataset = factor(dataset))

# Prepare statistical comparisons
sample_counts <- plot_data %>% 
  group_by(dataset) %>% 
  tally()

comparable_pairs <- combn(
  sample_counts %>% filter(n >= 2) %>% pull(dataset),
  2,
  simplify = FALSE
)

# Create plot with squared background
p <- ggplot(plot_data, aes(x = dataset, y = pop_pct)) +
  # Background with squares
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank()
  ) +
  
  # Mean bars with custom colors
  stat_summary(
    fun = mean, 
    geom = "bar", 
    aes(fill = dataset),
    width = 0.6, 
    alpha = 0.7
  ) +
  scale_fill_manual(values = custom_colors) +
  
  # Error bars
  geom_errorbar(
    data = filter(plot_data, dataset %in% sample_counts$dataset[sample_counts$n >= 2]),
    stat = "summary", 
    fun.data = mean_se, 
    width = 0.2
  ) +
  
  # Individual samples in black
  geom_quasirandom(
    color = "black",
    size = 3, 
    width = 0.2, 
    show.legend = FALSE
  ) +
  
  # Mean percentage labels
  stat_summary(
    aes(label = sprintf("%.1f%%", ..y..)),
    fun = mean, 
    geom = "text",
    vjust = -1.5, 
    size = 4, 
    fontface = "bold"
  ) +
  
  # Dynamic p-value additions
  {
    if(length(comparable_pairs) > 0) {
      stat_compare_means(
        comparisons = comparable_pairs,
        method = "t.test",
        label = "p.format",
        tip.length = 0.02,
        bracket.size = 0.6,
        step.increase = 0.1,
        size = 3.5
      )
    }
  } +
  
  # Labels and theme
  labs(
    title = paste0(marker_pos1, "+", marker_pos2, "+/", marker_neg, "- population"),
    subtitle = "Percentage across datasets",
    x = NULL,
    y = "Percentage of cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Add axis lines
  theme(
    axis.line.x = element_line(linewidth = 1.2, color = "black"),
    axis.line.y = element_line(linewidth = 1.2, color = "black")
  )

print(p)

# 8. Create the split DimPlot as in your example
split_plot <- DimPlot(t_cells,
                      group.by = "pos_pos_neg_neg",
                      split.by = "dataset",
                      cols = c("gray90", "red"),
                      order = TRUE) +
  ggtitle(paste(marker_pos1, "+", marker_pos2, "+ /", marker_neg, "- cells by Dataset")) +
  theme(legend.position = "top")

# Combine both plots if desired
combined_plots <- (p / split_plot) + plot_layout(heights = c(1, 2))
print(combined_plots)


# Define markers
cd28_neg <- "CD28"
cd57_pos <- "B3GAT1"

# Calculate cell states
t_cells$cd28_neg <- GetAssayData(t_cells, assay = "RNA")[cd28_neg, ] == 0
t_cells$cd57_pos <- GetAssayData(t_cells, assay = "RNA")[cd57_pos, ] > 0
t_cells$cd28_neg_cd57_pos <- t_cells$cd28_neg & t_cells$cd57_pos

# Get counts per dataset
count_data <- t_cells@meta.data %>%
  group_by(dataset) %>%
  summarise(
    total_cells = n(),
    cd28_neg_count = sum(cd28_neg),
    cd57_pos_count = sum(cd57_pos),
    cd28_neg_cd57_pos_count = sum(cd28_neg_cd57_pos),
    .groups = 'drop'
  ) %>%
  mutate(
    cd28_neg_pct = 100 * cd28_neg_count / total_cells,
    cd57_pos_pct = 100 * cd57_pos_count / total_cells,
    cd28_neg_cd57_pos_pct = 100 * cd28_neg_cd57_pos_count / total_cells
  )

# View results
print(count_data)

library(ggplot2)

ggplot(count_data, aes(x = dataset, y = cd28_neg_cd57_pos_count, fill = dataset)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cd28_neg_cd57_pos_count), vjust = -0.5) +
  labs(title = "Absolute Count of CD28- CD57+ Cells",
       y = "Number of Cells",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(count_data, aes(x = dataset, y = cd28_neg_cd57_pos_pct, fill = dataset)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", cd28_neg_cd57_pos_pct)), vjust = -0.5) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Percentage of CD28- CD57+ Cells",
       y = "% of Total Cells",
       x = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


# Assuming your count_data has:
# - dataset: categorical variable for x-axis
# - cd28_neg_cd57_pos_pct: percentage of CD28-CD57+ cells
# You'll need to calculate the "Other T cells" percentage (100 - cd28_neg_cd57_pos_pct)

# Reshape data to long format
plot_data <- count_data %>%
  mutate(Other_T_cells = 100 - cd28_neg_cd57_pos_pct) %>%
  pivot_longer(cols = c(cd28_neg_cd57_pos_pct, Other_T_cells),
               names_to = "population",
               values_to = "percentage")

# Create plot without percentage labels
ggplot(plot_data, aes(x = dataset, y = percentage, fill = population)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("cd28_neg_cd57_pos_pct" = "red", 
                               "Other_T_cells" = "gray70"),
                    labels = c("CD28- CD57+", "Other T cells")) +
  labs(title = "Percentage of CD28- CD57+ Cells Among T Cells",
       y = "% of Total Cells",
       x = "",
       fill = "Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

#-Percentage graph with 3 markers 

# Define markers
cd28_neg <- "CD28"
cd57_pos <- "B3GAT1"
klrg1_pos <- "KLRG1"

# Calculate cell states
t_cells$cd28_neg <- GetAssayData(t_cells, assay = "RNA")[cd28_neg, ] == 0
t_cells$cd57_pos <- GetAssayData(t_cells, assay = "RNA")[cd57_pos, ] > 0
t_cells$klrg1_pos <- GetAssayData(t_cells, assay = "RNA")[klrg1_pos, ] > 0

# Define triple-positive population (CD28- CD57+ KLRG1+)
t_cells$cd28_neg_cd57_pos_klrg1_pos <- t_cells$cd28_neg & t_cells$cd57_pos & t_cells$klrg1_pos

# Get counts per dataset
count_data <- t_cells@meta.data %>%
  group_by(dataset) %>%
  summarise(
    total_cells = n(),
    cd28_neg_count = sum(cd28_neg),
    cd57_pos_count = sum(cd57_pos),
    klrg1_pos_count = sum(klrg1_pos),
    cd28_neg_cd57_pos_klrg1_pos_count = sum(cd28_neg_cd57_pos_klrg1_pos),
    .groups = 'drop'
  ) %>%
  mutate(
    cd28_neg_pct = 100 * cd28_neg_count / total_cells,
    cd57_pos_pct = 100 * cd57_pos_count / total_cells,
    klrg1_pos_pct = 100 * klrg1_pos_count / total_cells,
    cd28_neg_cd57_pos_klrg1_pos_pct = 100 * cd28_neg_cd57_pos_klrg1_pos_count / total_cells
  )

# View results
print(count_data)

# Reshape data to long format
plot_data <- count_data %>%
  mutate(Other_T_cells = 100 - cd28_neg_cd57_pos_klrg1_pos_pct) %>%
  pivot_longer(
    cols = c(cd28_neg_cd57_pos_klrg1_pos_pct, Other_T_cells),
    names_to = "population",
    values_to = "percentage"
  )

# Create plot
ggplot(plot_data, aes(x = dataset, y = percentage, fill = population)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(
    values = c("cd28_neg_cd57_pos_klrg1_pos_pct" = "red", "Other_T_cells" = "gray70"),
    labels = c("CD28- CD57+ KLRG1+", "Other T cells")
  ) +
  labs(
    title = "Percentage of CD28- CD57+ KLRG1+ Cells Among T Cells",
    y = "% of Total Cells",
    x = "",
    fill = "Population"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

#---------------------------------------------------------------------------------------------------------------
#-T cells DEGs---------------------------------------------------------------------------------------------------

library(Seurat)
library(MAST)
library(purrr)
library(dplyr)
library(openxlsx)

# Define cell types
cell_types <- unique(t_cells$predicted.celltype.l2)
cell_types <- cell_types[!is.na(cell_types)]

# Run MAST for each cell type
all_mast_degs <- map_dfr(cell_types, function(ct) {
  subset_obj <- subset(t_cells, predicted.celltype.l2 == ct)
  
  # Check if "SCs" and "Young" exist in this subset
  if (!all(c("SCs", "Young") %in% unique(subset_obj$Age_range))) {
    message(sprintf("Skipping cell type %s: missing comparison groups", ct))
    return(NULL)
  }
  
  tryCatch({
    FindMarkers(
      subset_obj,
      ident.1 = "SCs",
      ident.2 = "Young",
      group.by = "Age_range",
      test.use = "MAST",
      latent.vars = "nCount_RNA",
      min.pct = 0.1,
      logfc.threshold = 0.25
    ) %>%
      mutate(
        celltype = ct,
        gene = rownames(.)  # Add gene names as a new column
      )
  }, error = function(e) {
    message(sprintf("Error in cell type %s: %s", ct, e$message))
    return(NULL)
  })
})

# Filter significant DEGs
significant_degs <- all_mast_degs %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  select(gene, everything())  # Put gene column first

# Save as Excel file
write.xlsx(
  significant_degs, 
  file = "MAST_DEGs_AgeRange_by_CellType.xlsx", 
  rowNames = FALSE,  # Don't write row names since we have gene column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")  # Make header bold
)

#-PATHWAYS for UPreg genes in T cells SCs------------------------------------------------------------------------

head(Genes_GOenrich_Downreg_Tcells)

library(clusterProfiler)
library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse
library(enrichplot)

org_db <- org.Hs.eg.db

# Initialize list to store results
pathway_results <- list()

# Loop through each cell type
for (ct in names(Genes_GOenrich_Downreg_Tcells)) {
  gene_list <- Genes_GOenrich_Downreg_Tcells[[ct]]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org_db)
  
  # Run KEGG or GO enrichment
  enrich_res <- enrichGO(gene         = entrez_ids$ENTREZID,
                         OrgDb        = org_db,
                         keyType      = "ENTREZID",
                         ont          = "BP",     # Biological Process
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  # Store result
  pathway_results[[ct]] <- enrich_res
}

# View top pathways for a cell type
head(pathway_results[["gdT"]])

# Dotplot
dotplot(pathway_results[["CD16 Mono"]], showCategory = 10) + ggtitle("CD14 Mono Pathway Enrichment")

# Barplot
barplot(enrich_res, showCategory = 10)

#---------------------------------------------------------------------------------------------------
#-Others hallmarks of Aging-------------------------------------------------------------------------

#-Young_SCs_integrated (ALL celltypes) DEGs-------------------------------------------------------------------------

library(Seurat)
library(MAST)
library(purrr)
library(dplyr)
library(openxlsx)

# Define cell types
cell_types <- unique(Young_SCs_integrated$predicted.celltype.l2)
cell_types <- cell_types[!is.na(cell_types)]

# Run MAST for each cell type
all_mast_degs <- map_dfr(cell_types, function(ct) {
  subset_obj <- subset(Young_SCs_integrated, predicted.celltype.l2 == ct)
  
  # Check if "SCs" and "Young" exist in this subset
  if (!all(c("SCs", "Young") %in% unique(subset_obj$Age_range))) {
    message(sprintf("Skipping cell type %s: missing comparison groups", ct))
    return(NULL)
  }
  
  tryCatch({
    FindMarkers(
      subset_obj,
      ident.1 = "SCs",
      ident.2 = "Young",
      group.by = "Age_range",
      test.use = "MAST",
      latent.vars = "nCount_RNA",
      min.pct = 0.1,
      logfc.threshold = 0.25
    ) %>%
      mutate(
        celltype = ct,
        gene = rownames(.)  # Add gene names as a new column
      )
  }, error = function(e) {
    message(sprintf("Error in cell type %s: %s", ct, e$message))
    return(NULL)
  })
})

# Filter significant DEGs
significant_degs <- all_mast_degs %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  select(gene, everything())  # Put gene column first

# Filter significant DEGs
significant_degs <- all_mast_degs %>%
  dplyr::as_tibble() %>%  # Explicitly use dplyr's as_tibble
  dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  dplyr::select(gene, dplyr::everything())  # Explicit namespace

# Save as Excel file
write.xlsx(
  significant_degs, 
  file = "MAST_DEGs_AgeRange_by_CellType_ALL.xlsx", 
  rowNames = FALSE,  # Don't write row names since we have gene column
  firstRow = TRUE,
  headerStyle = createStyle(textDecoration = "bold")  # Make header bold
)

#-PATHWAYS for UPreg genes in SCs-------------------------------------------------------------------------------

head(GOenrich_upregSCs)

library(clusterProfiler)
library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse
library(enrichplot)

org_db <- org.Hs.eg.db

# Initialize list to store results
pathway_results <- list()

# Loop through each cell type
for (ct in names(GOenrich_downregSCs)) {
  gene_list <- GOenrich_downregSCs[[ct]]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org_db)
  
  # Run KEGG or GO enrichment
  enrich_res <- enrichGO(gene         = entrez_ids$ENTREZID,
                         OrgDb        = org_db,
                         keyType      = "ENTREZID",
                         ont          = "MF",     # Biological Process
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  # Store result
  pathway_results[[ct]] <- enrich_res
}

# View top pathways for a cell type
head(pathway_results[["gdT"]])

# Dotplot
dotplot(pathway_results[["CD16 Mono"]], showCategory = 10) + ggtitle("CD14 Mono Pathway Enrichment")

# Barplot
barplot(enrich_res, showCategory = 10)


#-KEGG new way for GO enrich--------------------------------------------------------------------------------------

library(dplyr)

# Loop through each cell type
for (ct in names(GOenrich_downregSCs)) {
  gene_list <- GOenrich_downregSCs[[ct]]
}

library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db if you're working with mouse data

entrez_up <- bitr(gene_list, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# KEGG enrichment analysis
kegg_res <- enrichKEGG(gene = entrez_up$ENTREZID,
                       organism = 'hsa',  # 'mmu' for mouse
                       pvalueCutoff = 0.05)

# View top pathways
head(kegg_res)

library(ReactomePA)

reactome_res <- enrichPathway(gene = entrez_up$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)

# View enriched pathways
head(reactome_res)

# KEGG
barplot(kegg_res, showCategory = 10, title = "KEGG Pathways")

# Reactome
barplot(reactome_res, showCategory = 10, title = "Reactome Pathways")

# Filter descriptions containing your terms
terms_of_interest <- c("immune", "cytokine", "aging", "lifespan", "senescence", "autophagy")

# Get matching descriptions
filtered_terms <- as.data.frame(reactome_res) %>%
  dplyr::filter(grepl(paste(terms_of_interest, collapse = "|"), Description, ignore.case = TRUE))

# Plot only filtered pathways manually (custom plot)
library(ggplot2)
ggplot(filtered_terms, aes(x = reorder(Description, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Filtered Reactome Pathways",
       x = NULL,
       y = "Gene Count") +
  theme_minimal()

ggplot(filtered_terms, aes(x = reorder(Description, -Count), y = Count, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log10(adj p-value)") +
  labs(title = "Filtered Reactome Pathways",
       x = NULL,
       y = "Gene Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14))

#-For GSEA------------------------------------------------------------------------------------------
# Named vector of logFCs
gene_list <- Table_SX_KEGG$avg_log2FC
names(gene_list) <- Table_SX_KEGG$gene

# Convert gene names to Entrez
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list) <- gene_list_entrez$ENTREZID

# GSEA using KEGG
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = 'hsa',
                     nPerm = 1000,
                     minGSSize = 10,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

dotplot(gsea_kegg, showCategory = 10)








