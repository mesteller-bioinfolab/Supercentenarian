# ==========================
# Required Libraries
# ==========================
library(readxl)         # To read Excel files
library(openxlsx)       # To write and modify Excel workbooks
library(NOISeq)         # For normalization and batch correction (e.g., ARSyNseq)
library(ggplot2)        # For plotting (PCA, etc.)
library(ggfortify)      # For autoplot PCA
library(limma)          # For differential expression analysis
library(clusterProfiler) # For GO enrichment analysis
library(org.Hs.eg.db)   # Human gene annotation for GO analysis
library(GOSemSim)       # To compute semantic similarity between GO terms
library(dplyr)          # Data manipulation
library(tidyr)          # For pivot_longer and reshaping data
library(circlize)       # For circos plot
library(tibble)         # For tidy data frames


# ================================
# Input Data Format Requirements
# ================================
#
# The input file should be an excel or CSV file structured as follows:
#
# 1. **First Column**: Contains *protein abbreviations* or identifiers 
#    (e.g., UniProt IDs, gene symbols).
#
# 2. **Remaining Columns**: Contain *quantitative expression values* for each sample.
#    - Each column must correspond to a biological replicate or condition.
#    - Column names should be the *group labels* (e.g., "M116" and "Ypm").
#      In this case: 6 columns for "M116" and 24 columns for "Ypm" (each woman in sextuplicate).
#
# ⚠ Data cleaning steps (applied before analysis):
# - **Replace zero values**:
#     --> Any zero (`0`) value should be replaced by (row-wise minimum non-zero value) / 2.
#
# - **Remove non-quantified proteins**:
#     --> Rows corresponding to proteins that were detected but not quantified 
#         (i.e., all values are zero or NA) should be excluded from the analysis.
#

# =======================================================================
# Cleaning, Replicate Equalization, batch effect and TMM Normalization
# =======================================================================
# Read data from Excel
data <- read_excel("file.xlsx", sheet = "data", col_names = TRUE)
data <- as.data.frame(data)
original_data <- data  # Save copy of original data

# Prepare treatment labels and protein names
treatment <- colnames(data)
treatment <- treatment[treatment != "proteins"]
protein_names <- data$proteins

# Set rownames and remove the 'proteins' column
rownames(data) <- make.unique(as.character(protein_names))
data <- data[, -1]

# Convert all values to numeric (if needed)
data <- as.data.frame(lapply(data, function(x) as.numeric(as.character(x))))

# Clean column names and define treatment groups
treatment <- gsub("\\.\\.\\.[0-9]+", "", treatment)
myfactors <- data.frame(Treatment = treatment)
myfactors$Treatment <- factor(myfactors$Treatment, levels = c("M116", "Ypm"))

# Create NOISeq object and perform TMM normalization
# ⚠ If needed, add batch effect correction by setting `batch = TRUE`
mydata <- readData(data = data, factors = myfactors)
mydatacorr <- ARSyNseq(mydata, batch = FALSE, beta = 1.5, norm = "TMM", logtransf = FALSE)

# Extract normalized expression matrix and log transform
corrected_TMM <- assayData(mydatacorr)$exprs
corrected_TMM <- log2(corrected_TMM + 1)
rownames(corrected_TMM) <- rownames(data)
colnames(corrected_TMM) <- colnames(data)

log2_data_with_proteins <- cbind(Protein = rownames(corrected_TMM), corrected_TMM)
# Write both original and normalized data to Excel
wb <- createWorkbook()
addWorksheet(wb, "Original_data")
writeData(wb, sheet = "Original_data", original_data)

addWorksheet(wb, "TMM_normalized_data")
writeData(wb, sheet = "TMM_normalized_data", log2_data_with_proteins)

saveWorkbook(wb, "file.xlsx", overwrite = TRUE)

# =============
# PCA Analysis 
# =============
#
# This section performs a PCA to check technical and biological replicate distribution and group distribution.
# If any replicate clusters abnormally (e.g., outlier or technical issue),
# it should be removed and the normalization step repeated.

df <- read.xlsx("file.xlsx", sheet = "TMM_normalized_data")

# Define treatment groups 
treatments <- c(rep("M116", 6),
                rep("yPM1", 6),
                rep("yPM2", 6),
                rep("yPM3", 6),
                rep("yPM4", 6))  # 24 YPM replicates total

# Transpose data for PCA
a <- df[, -1]  # Remove 'Protein' column
df_transposed <- as.data.frame(t(a))
colnames(df_transposed) <- df$Protein
rownames(df_transposed) <- make.unique(colnames(df)[-1])

# Add group labels
df_transposed$condition <- as.character(treatments)


numeric_data <- as.data.frame(sapply(df_transposed[, -ncol(df_transposed)], as.numeric))
numeric_data$condition <- df_transposed$condition

# Remove proteins with zero variance 
zero_var_cols <- apply(numeric_data[, -ncol(numeric_data)], 2, var, na.rm = TRUE) == 0
filtered_data <- numeric_data[, !zero_var_cols]
filtered_data$condition <- numeric_data$condition

# Perform PCA
pca_all <- prcomp(filtered_data[, -ncol(filtered_data)], scale. = TRUE)
pca_plot_all <- autoplot(pca_all, data = filtered_data, colour = 'condition', size = 2,
                         loadings = FALSE, loadings.label = FALSE) +
  ggtitle("PCA – TMM Normalized Data") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

print(pca_plot_all)


# ============================================================
# Averaging Technical Replicates – TMM Normalized Data
# ============================================================

# Load normalized data
data <- read.xlsx("file.xlsx", sheet = "TMM_normalized_data")

# Extract protein names and convert expression values to numeric
averaged_results <- data.frame(Abbreviation = data[[1]])
data_matrix <- data[, -1]
data_matrix[] <- lapply(data_matrix, function(x) as.numeric(as.character(x)))

# Define biological group labels (5 groups × 6 replicates)
group_names <- c("M116", "Ypm1", "Ypm2", "Ypm3", "Ypm4")
replicates_per_group <- 6
expected_cols <- length(group_names) * replicates_per_group

# Check column count
if (ncol(data_matrix) != expected_cols) {
  stop("Mismatch: expected ", expected_cols, " data columns for ", length(group_names), " groups × ", replicates_per_group, " replicates.")
}

# Average every 6 columns
for (i in seq(1, ncol(data_matrix), by = replicates_per_group)) {
  group_data <- data_matrix[, i:(i + replicates_per_group - 1)]
  group_avg <- rowMeans(group_data, na.rm = TRUE)
  group_index <- (i - 1) / replicates_per_group + 1
  group_label <- group_names[group_index]
  averaged_results[[group_label]] <- group_avg
}

# Save to new sheet in Excel
wb <- loadWorkbook("file.xlsx")
addWorksheet(wb, "TMM_normalized_averaged")
writeData(wb, sheet = "TMM_normalized_averaged", averaged_results)
saveWorkbook(wb, "file.xlsx", overwrite = TRUE)



# =================================
# Differential expression analysis
# =================================

# This block performs differential expression analysis using the averaged TMM-normalized values.
# It compares the M116 group versus Ypm using limma's linear modeling.

data <- read.xlsx("file.xlsx", sheet = "TMM_normalized_averaged")

# Set the row names to the values from the 'Abbreviation' column
protein_names <- data$Abbreviation


# Remove the 'Abbreviation' column (first column) from 'data' for analysis
data <- data[,-1]

# Extract the treatment information from the column names (excluding 'Abbreviation')

treatment <-c("M116", "Ypm", "Ypm", "Ypm", "Ypm")
data <- as.data.frame(lapply(data, as.numeric))

# Create 'myfactors' data frame for treatments
myfactors <- data.frame(Treatment = treatment)
myfactors$Treatment <- factor(myfactors$Treatment, levels = c("M116", "Ypm"))

# Update the NOISeq object with TMM normalized data
# Ensure that the data passed to readData has no row names and only numeric data
mydata2 <- readData(data = as.matrix(data), factors = myfactors)

# Create the design matrix for limma
# Ensure that the number of rows in the design matrix matches the number of samples (columns) in mydata2
design <- model.matrix(~ 0 + treatment, data = myfactors)  
# Check the dimensions of both mydata2 and design
colnames(design) <- sub("^treatment", "", colnames(design))
# Fit linear model using limma
fit <- lmFit(mydata2, design)
# Define contrasts for pairwise comparisons between all Ypm groups and M116 
contrast.matrix <- makeContrasts(
  M116_vs_Ypm = M116 - Ypm,
  
  levels = design
)

# Apply contrasts to the linear model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract top table 
result <- topTable(fit2, coef = 1, adjust = "BH", number = Inf, genelist = protein_names)

# Filtered results: adjusted p-value < 0.05 and |logFC| > 1
filtered_result <- subset(result, adj.P.Val < 0.05 & abs(logFC) > 1)

# Filtered results: raw p-value < 0.05 and |logFC| > 1
filtered_notadjusted <- subset(result, P.Value < 0.05 & abs(logFC) > 1)


addWorksheet(wb, "results_M116vsYpm")
addWorksheet(wb, "filtered_M116vsYpm")
addWorksheet(wb, "filt(notAdj)_M116vsYpm")

writeData(wb, sheet = "results_M116vsYpm", result)
writeData(wb, sheet = "filtered_M116vsYpm", filtered_result)
writeData(wb, sheet = "filt(notAdj)_M116vsYpm", filtered_notadjusted)

# Save the workbook
saveWorkbook(wb, "file.xlsx", overwrite = TRUE)



# ======================================================
# Gene Ontology-Biological Processes and clustering
# ======================================================
#Clustering: Use the dendrogram as a reference to group enriched GO terms into higher-level biological function clusters.
#You may define the number of clusters visually (by cutting the tree) or based on domain knowledge.
##Manual Refinement: After auto-clustering, check whether GO term descriptions align biologically within clusters.If not, manually reassign terms by biological relevance.
#Terms Excluded from Similarity Matrix: Some GO terms may not appear in the similarity matrix due to missing semantic data. You may manually assign these to relevant clusters based on their biological meaning and enrichment context.

# Load data from the Excel file
df <- read.xlsx("file.xlsx", sheet = "filtered_M116vsYpm")#significant results obtained with the differential analysis

# Step 2: Extract the gene list
gene_list <- df$ID

# Perform GO enrichment analysis
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",   # Use Entrez Gene IDs
                ont = "BP",           # Biological Process
                pvalueCutoff = 0.05,  
                qvalueCutoff = 0.05)
ego_results<-as.data.frame(ego)
write.xlsx(ego_results, "EGO_results.xlsx") # save GO results in a file

d <- godata('org.Hs.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method="Wang", semData = d)

# Extract similarity matrix from the pairwise_termsim result
sim_matrix <-as.matrix(ego2@termsim)
hc <- hclust(as.dist(sim_matrix), method = "complete")  

# Plot the dendrogram
plot(hc, main = "Dendrogram of GO Terms", xlab = "GO Terms", sub = "", cex = 0.7)
go_ordered_ids <- hc$labels[hc$order]
go_ordered_ids<-as.data.frame(go_ordered_ids)
write.xlsx(go_ordered_ids,"sim_matrix_orderedID.xlsx") #file with the GO IDs in dendrogram order

#Use the dendrogram as a reference to group enriched GO terms into higher-level biological function clusters
#We clustered them into 8:
   ##Cluster 1: Coagulation, Hemostasis, and Wound Healing
   ##Cluster 2: Immune System and Inflammatory Responses
   ##Cluster 3: Lipid and Cholesterol Transport and Metabolism
   ##Cluster 4: Apoptotic Signaling and Cellular Response to Stress
   ##Cluster 5: Protein Processing and Translation
   ##Cluster 6: Cellular  Detoxification
   ##Cluster 7: Cell Adhesion and Structural Organization
   ##Cluster 8: Regulation of mRNA 

# ========
# FIGURES
# ========

#Create a new sheet in file.xlsx with the Genes, log fold change, and clusters-in each cluster column there should be the geneID from the GO results corresponding to each term of each cluster even if repeated
#We had the following columns: gene, logFC, Cluster1, Cluster2, Cluster3, Cluster4, Cluster5, Cluster6, Cluster7, Cluster8.

# Read the Excel file
data <- read.xlsx("file.xlsx", sheet = "FIGURES")

# ============
# CIRCOS PLOT
# ============

# Identify cluster columns
cluster_columns <- grep("^Cluster", names(data), value = TRUE)

# Reshape the data
data_long <- data %>%
  pivot_longer(
    cols = all_of(cluster_columns),
    names_to = "cluster",
    values_to = "gene_in_cluster"
  ) %>%
  filter(!is.na(gene_in_cluster))

# Map logFC values correctly
final_data <- data_long %>%
  mutate(logFC = data$logFC[match(gene_in_cluster, data$Gene)]) %>%
  rename(gene = gene_in_cluster)

# Select desired columns
final_data <- final_data[, c("gene", "logFC", "cluster")]

# Sort the data by logFC in descending order
final_data <- final_data %>%
  arrange(desc(logFC))

# Add gene counts
gene_counts <- table(final_data$gene)
final_data$gene_count <- gene_counts[final_data$gene]

# Define color ramp for logFC
library(circlize)
col_fun <- colorRamp2(
  breaks = c(min(final_data$logFC, na.rm = TRUE), 0, max(final_data$logFC, na.rm = TRUE)),
  colors = c("blue", "white", "red")
)

# Map gene colors consistently based on unique genes and their logFC values
unique_genes <- unique(final_data$gene)
logFC_values <- sapply(unique_genes, function(g) {
  unique(final_data$logFC[final_data$gene == g])
})
gene_colors <- setNames(col_fun(logFC_values), unique_genes)

# Define cluster colors
unique_clusters <- unique(final_data$cluster)
cluster_colors <- c(
  "Cluster1" = "#BAAE00",  
  "Cluster2" = "gold",  
  "Cluster3" = "#B3E0A6",  # Green
  "Cluster4" = "chocolate",  # Red
  "Cluster5" = "#E982B5",   # Purple
  "Cluster6" = "grey",  # Green
  "Cluster7" = "#17becf",  # Red
  "Cluster8" = "#B4219C"   # Purple
)

# Combine colors
grid.col <- c(gene_colors, cluster_colors)

# Create link colors
final_data$link_col <- cluster_colors[final_data$cluster]

# Ensure minimum sector width
factor_widths <- pmax(
  c(as.numeric(table(final_data$gene)) / sum(table(final_data$gene)), 
    rep(0.1, length(unique_clusters))),
  0.005
)
factor_widths <- factor_widths / sum(factor_widths)

# Calculate custom gap degrees
gene_count <- length(unique(final_data$gene))
cluster_count <- length(unique(final_data$cluster))

# Uniform gap between genes and clusters, with larger gap between the two groups
gap_degrees <- c(
  rep(0.3, gene_count - 1),  # Uniform gap between genes
  6,                       # Larger gap between genes and clusters
  rep(0.5, cluster_count - 1),  # Uniform gap between clusters
  5                        # Larger gap to close the circle
)

# Initialize circos parameters with custom gap degrees
circos.clear()
circos.par(
  "gap.degree" = gap_degrees,  # Apply custom gaps
  "track.margin" = c(0.005, 0.005),
  "start.degree" = 90
)
# Create chord diagram with genes ordered by logFC

custom_cluster_order <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", 
                          "Cluster5", "Cluster6", "Cluster7", "Cluster8")
final_data$cluster <- factor(final_data$cluster, levels = custom_cluster_order)

# Ensure that the cluster colors are ordered according to the custom order
cluster_colors <- cluster_colors[custom_cluster_order]
chordDiagram(
  x = final_data[, c("gene", "cluster")],  # Gene-cluster relationships
  grid.col = grid.col,                     # Colors for sectors
  col = final_data$link_col,               # Colors for links
  transparency = 0.7,                        # Transparency for the links
  annotationTrack = "grid",                # Include annotation track grid
  annotationTrackHeight = c(0.1, 0.2),     # Height for annotation tracks
  order = c(unique(final_data$gene), custom_cluster_order),  # Order genes and clusters
  directional = 1,                         # Add directional links
  diffHeight = mm_h(8.5),
  target.prop.height = mm_h(8),
  scale = TRUE)

# Add labels using circos.labels
circos.labels(
  sectors = unique(final_data$gene),  # Genes as sectors
  labels = unique(final_data$gene),   # Labels for each sector (genes)
  side = "outside",                   # Place labels outside the circle
  x = rep(1, length(unique(final_data$gene))),  # Position labels outside the circle
  cex = 0.38,                          # Adjust font size
  niceFacing = TRUE                    # Adjust label facing for readability
)


## LEGEND##

# Legend for Fold Change Ramp
col_fun = colorRamp2(c(min(final_data$logFC, na.rm = TRUE), 0, max(final_data$logFC, na.rm = TRUE)), c("blue", "white", "red"))
fc_lgd <- Legend(
  col_fun = col_fun, 
  title = "logFC",
  at = c(-5, -2.5, 0, 2.5, 5)  # Force legend to display full range
)

# Define cluster colors
cluster_colors <- c(
  "Cluster 1" = "#BAAE00",  
  "Cluster 2" = "gold",  
  "Cluster 3" = "#B3E0A6",  
  "Cluster 4" = "chocolate",  
  "Cluster 5" = "#E982B5",   
  "Cluster 6" = "grey",  
  "Cluster 7" = "#17becf",  
  "Cluster 8" = "#B4219C"   
)

# Cluster descriptions
cluster_labels <- c(
  "Coagulation, Hemostasis, and Wound Healing",
  "Immune System and Inflammatory Responses",
  "Lipid and Cholesterol Transport and Metabolism",
  "Apoptotic Signaling and Cellular Response to Stress",
  "Protein Processing and Translation",
  "Cellular Detoxification",
  "Cell Adhesion and Structural Organization",
  "Regulation of mRNA"
)

# Legend for Cluster Colors with Square Indicators
cluster_lgd <- Legend(
  at = names(cluster_colors),  # Cluster names
  labels = cluster_labels,     # Corresponding descriptions
  legend_gp = gpar(fill = cluster_colors),  # Use your defined cluster_colors
  labels_gp = gpar(fontsize = 10),  # Adjust font size for labels
  title = "GO Biological Process Clusters",
  type = "grid"  # Square color indicators
)

# Combine Legends and Display
draw(packLegend(fc_lgd, cluster_lgd))


# ============
# PIE PLOT
# ============

# Count all gene-cluster occurrences
cluster_counts <- final_data %>%
  group_by(cluster) %>%
  summarise(count = n(), .groups = "drop")  # n() counts rows

# Calculate percentages
cluster_counts <- cluster_counts %>%
  mutate(percentage = (count / sum(count)) * 100)

# Define order
custom_order <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", 
                  "Cluster5", "Cluster6", "Cluster7", "Cluster8")

cluster_counts <- cluster_counts %>%
  mutate(cluster = factor(cluster, levels = custom_order)) %>%
  arrange(cluster)

# Generate labels with percentages
labels <- paste0(round(cluster_counts$percentage, 1), "%")

colors <- cluster_colors[(cluster_counts$cluster)]
pie(
  x = cluster_counts$count,
  labels = labels,
  clockwise=TRUE,
  col = colors,
  main = "Cluster Distribution"
)


