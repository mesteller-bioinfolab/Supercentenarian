#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models

saveResults = f"/mnt/beegfs/mcasado/PROJECTS/07.2_Supercentenarian/results/Terekhova_et_al/02_celltypist"

adata = sc.read_h5ad("/mnt/beegfs/mcasado/PROJECTS/07.2_Supercentenarian/processed_data/Terekhova_et_al/01_subsampling/all_pbmcs_rna_subsampled_10000_by_ag.h5ad")
check = adata.X.expm1().sum(axis = 1)
print(check[:5])

predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
predictions.predicted_labels.to_csv(f"{saveResults}/celltypist_pred_Immune_All_Low_majority_voting.csv")

predictions = celltypist.annotate(adata, model = 'Immune_All_High.pkl', majority_voting =  True)
predictions.predicted_labels.to_csv(f"{saveResults}/celltypist_pred_Immune_All_High_majority_voting.csv")

