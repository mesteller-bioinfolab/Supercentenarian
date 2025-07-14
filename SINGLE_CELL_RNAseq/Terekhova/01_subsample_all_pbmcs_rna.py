import scanpy as sc
import pandas as pd
import numpy as np


# Load the data
adata = sc.read_h5ad("/mnt/beegfs/mcasado/PROJECTS/07.2_Supercentenarian/raw_data/Terekhova_et_al/all_pbmcs/all_pbmcs_rna.h5ad")
print("DATA LOADED")
print(adata.var)

print(adata.obs.head(10))
print(adata.obs.columns)
print("adata.var_names[:5]:")
print(adata.var_names[:5])

print("adata.var")
print(adata.var)

# Print the first 5 cell barcodes
print("First 5 cell barcodes from adata.obs:")
print(adata.obs.index[:5])

# Print the first 5 gene names
print("First 5 gene names from adata.var:")
print(adata.var.index[:5])


# Load the metadata, ensuring the index is read as a string
meta = pd.read_csv('/mnt/beegfs/mcasado/PROJECTS/07.2_Supercentenarian/raw_data/Terekhova_et_al/all_pbmcs/all_pbmcs_metadata.csv', index_col=0, dtype=str)
print("METADATA LOADED")
print(meta.head())

# Ensure indices are strings for both adata.obs and meta
adata.obs.index = adata.obs.index.astype(str)
meta.index = meta.index.astype(str)

# Filter adata to keep only cells present in meta
cells_to_keep = adata.obs.index.intersection(meta.index)
adata = adata[cells_to_keep].copy()

# Merge the metadata into adata.obs
adata.obs = adata.obs.merge(meta, left_index=True, right_index=True, how='left')

print("METADATA ADDED")
print(adata.obs.head(10))
print(adata.obs.columns)
print(f"Number of cells in adata: {adata.n_obs}")

# Subsample
subsampled_groups = []

for age_group in adata.obs['Age_group'].unique():
    # Filter for current age group
    current_group_data = adata[adata.obs['Age_group'] == age_group]
    
    # Check if the current group has at least 10,000 cells
    if current_group_data.n_obs > 10000:
        # Subsample 50,000 cells from the current group
        sc.pp.subsample(current_group_data, n_obs=10000)
    # If the group has fewer than 10,000 cells, just take them all without subsampling
    subsampled_groups.append(current_group_data)

# Concatenate all the subsampled groups back into a single AnnData object
adata_subsampled = sc.concat(subsampled_groups, join='outer')

print(f"Number of cells in adata: {adata_subsampled.n_obs}")

counts = adata_subsampled.obs.groupby(['Donor_id', 'Age_group']).size().unstack(fill_value=0)
print(counts.head())

counts = adata_subsampled.obs.groupby(['Age_group', 'Donor_id']).size().unstack(fill_value=0)
print(counts.head())



adata_subsampled.write_h5ad('/mnt/beegfs/mcasado/PROJECTS/07.2_Supercentenarian/processed_data/Terekhova_et_al/01_subsampling/all_pbmcs_rna_subsampled_10000_by_ag.h5ad')

# Print the first 5 cell barcodes
print("First 5 cell barcodes from adata_subsampled.obs:")
print(adata_subsampled.obs.index[:5])

# Print the first 5 gene names
print("First 5 gene names from adata_subsampled.var:")
print(adata_subsampled.var.index[:5])


