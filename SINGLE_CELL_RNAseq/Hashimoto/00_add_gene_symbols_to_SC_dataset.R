library(tidyverse)
library(Seurat)
library(biomaRt)

dataset <- "Hashimoto_et_al"
step <-  "0_add_gene_symbols"

projectRoot <- "~/Documents/PROJECTS/07.2_Supercentenaria/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)


# Read data (takes a long time to load)
# --------------------------------------------------
raw_data_path <- file.path(projectRoot, "raw_data", dataset)

umi_counts <- read.table(file.path(raw_data_path, "01.UMI.txt"), header = TRUE, row.names = 1, sep = "\t")
cell_barcodes <- read.table(file.path(raw_data_path, "03.Cell.Barcodes.txt"), header = FALSE)
data <- CreateSeuratObject(counts = umi_counts, project = "supercentenarians")

saveRDS(data, file.path(saveProcessedData, "merged_obj.RDS"))
# --------------------------------------------------

data <- readRDS(file.path(saveProcessedData, "merged_obj.RDS"))



# CONVERT ENSMBL IDS TO SYMBOL IDS
#-------------------------------------------------------------------------------
#Select the Ensembl database and specify the dataset for Mus musculus (mouse)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for gene names, removing version numbers from the IDs
genes <- rownames(data)
length(unique(genes))
result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id", 
                values = genes, 
                mart = mart)

print(head(result))
length(unique(result$hgnc_symbol))
result <- result %>% mutate(hgnc_symbol_unique = make.unique(as.character(result$hgnc_symbol), sep = "_"))
length(unique(result$hgnc_symbol))
saveRDS(result, file.path(saveProcessedData, "result.RDS"))
print(head(result))

matrix <- data@assays$RNA



raw_counts_gene_symbol <- raw_counts %>% inner_join(result, by = c("ensembl_gene_id_v2"="ensembl_gene_id"))
raw_counts_gene_symbol <- raw_counts_gene_symbol %>% as_tibble() %>% 
  dplyr::select(-c(ensembl_gene_id_v2, mgi_symbol, ensembl_gene_id)) %>% as.data.frame()

rownames(raw_counts_gene_symbol) <- raw_counts_gene_symbol$mgi_symbol_unique
raw_counts_gene_symbol$mgi_symbol_unique <-  NULL

write.csv(raw_counts_gene_symbol, file.path(projectRoot, "raw_data", dataset, "raw_counts_gene_symbol.csv"), row.names = TRUE)

print(dim(raw_counts))
print(dim(raw_counts_gene_symbol))

#-------------------------------------------------------------------------------



# Define Ensembl connection
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve annotations
ensembl_gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                              filters = 'ensembl_gene_id', 
                              values = rownames(data), 
                              mart = ensembl.con)

# Filter and ensure unique mappings
annotated_data <- ensembl_gene_symbols %>%
  filter(external_gene_name != "" & !is.na(external_gene_name)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Subset data for valid Ensembl IDs
valid_ensembl_ids <- annotated_data$ensembl_gene_id
data <- subset(data, features = valid_ensembl_ids)

# Update rownames to gene symbols
gene_symbols_ordered <- annotated_data$external_gene_name[
  match(rownames(data), annotated_data$ensembl_gene_id)
]
rownames(data) <- gene_symbols_ordered

# Validate barcode integration
unmatched_barcodes <- setdiff(data$barcodes, cell_barcodes$barcodes)
cat("Number of unmatched barcodes:", length(unmatched_barcodes), "\n")

# Integrate metadata
data@meta.data <- data@meta.data %>%
  inner_join(cell_barcodes, by = "barcodes")











#--------











ensembl_ids <- rownames(data)
# convert ensemble ids to gene symbols
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
ensembl_gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                              filters = 'ensembl_gene_id', 
                              values = ensembl_ids,
                              mart = ensembl.con)

matched_indices <- match(ensembl_ids, ensembl_gene_symbols$ensembl_gene_id)
is.na(matched_indices) %>% table()
# FALSE  TRUE 
# 22979   405 
matched_indices <- na.omit(matched_indices)
ensembl_gene_symbols <- ensembl_gene_symbols %>%
  dplyr::slice(matched_indices)

ensembl_gene_symbols <- ensembl_gene_symbols %>%
  filter(external_gene_name != "" & !is.na(external_gene_name))
ensembl_gene_symbols <- ensembl_gene_symbols %>% mutate(gene_symbol_duplicated = duplicated(external_gene_name))
ensembl_gene_symbols <- ensembl_gene_symbols %>%
  filter(gene_symbol_duplicated == FALSE)
valid_ensembl_ids <- ensembl_gene_symbols$ensembl_gene_id

data <- subset(data, features = valid_ensembl_ids)
gene_symbols_ordered <- ensembl_gene_symbols$external_gene_name[match(rownames(data), ensembl_gene_symbols$ensembl_gene_id)]

rownames(data) <- gene_symbols_ordered


data$barcodes <- rownames(data@meta.data)
cell_barcodes <- cell_barcodes %>% dplyr::rename(Donor_id = V2) %>% dplyr::rename(barcodes = V1) %>% dplyr::rename(condition = V3)
cell_barcodes <- cell_barcodes %>% mutate(barcodes = str_replace_all(barcodes, "-", "."))
data@meta.data <- data@meta.data %>% inner_join(cell_barcodes, by = "barcodes")
rownames(data@meta.data) <- data$barcodes 

# Save object
saveRDS(data, "~/Documents/PROJECTS/07_Supercentenaria/public_studies/processed-data/supercentenarians/SC_counts_with_gene_symbols.RDS")
