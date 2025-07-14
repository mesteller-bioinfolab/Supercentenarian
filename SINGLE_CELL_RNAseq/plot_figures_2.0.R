library(Seurat)
library(tidyverse)  

dataset <- "figures"

projectRoot <- "~/Eva/IJC/Supercentenarian/07.2_Supercentenaria"
source(file.path(projectRoot, "misc", "colors.R"))

savePlots <- file.path(projectRoot, "plots", dataset); dir.create(savePlots, recursive = T)

# M116
data_path_MM16 <- file.path(projectRoot, "processed_data", "M116", "04.3.add_CellTyspist_pred", "data.RDS")
data_M116 <- readRDS(data_path_MM16)

read.csv(file.path("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/07.2_Supercentenarian/results/Terekhova_et_al/02_celltypist/celltypist_pred_Immune_All_High_majority_voting.csv"))
read.csv(file.path("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/07.2_Supercentenarian/results/Terekhova_et_al/02_celltypist/celltypist_pred_Immune_All_Low_majority_voting.csv"))

# Terekhova_et_al
celltypist_high_Terekhova_et_al <- read.csv(file.path(projectRoot, "Terekhova_et_al_ANALYSIS_HPC", "07.2_Supercentenarian", "results", "Terekhova_et_al", "02_celltypist", "celltypist_pred_Immune_All_High_majority_voting.csv")) 
celltypist_low_Terekhova_et_al <- read.csv(file.path(projectRoot, "Terekhova_et_al_ANALYSIS_HPC", "07.2_Supercentenarian", "results", "Terekhova_et_al", "02_celltypist", "celltypist_pred_Immune_All_Low_majority_voting.csv")) 
colnames(celltypist_low_Terekhova_et_al) <- paste0(colnames(celltypist_low_Terekhova_et_al), "_Immune_All_Low")
colnames(celltypist_high_Terekhova_et_al) <- paste0(colnames(celltypist_high_Terekhova_et_al), "_Immune_All_High")
celltypist_low_Terekhova_et_al$majority_voting_Immune_All_Low <- gsub("/", " or ", celltypist_low_Terekhova_et_al$majority_voting_Immune_All_Low)
celltypist_high_Terekhova_et_al$majority_voting_Immune_All_High  <- gsub("/", " or ", celltypist_high_Terekhova_et_al$majority_voting_Immune_All_High)

# Hashimoto_et_al
#data_path_Hashimoto <- file.path(projectRoot, "processed_data", "Hashimoto_et_al", "04.3.add_CellTyspist_pred", "data.RDS")
data_Hashimoto <- readRDS("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/07.2_Supercentenarian/processed_data/Hashimoto_et_al/04.3.add_CellTyspist_pred/data.RDS")


# Barplot proportions
df_M116 <- data_M116@meta.data %>% dplyr::select(majority_voting_Immune_All_Low) %>% mutate(sample = "M116") 
df_Hashimoto <- data_Hashimoto@meta.data %>% dplyr::select(majority_voting_Immune_All_Low, sample) %>% 
  filter(!is.na(majority_voting_Immune_All_Low))

read_csv("~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/Terekhova_et_al_ANALYSIS_HPC/07.2_Supercentenarian/raw_data/Terekhova_et_al/all_pbmcs/all_pbmcs/all_pbmcs_metadata.csv")
celltypist_low_Terekhova_et_al <- celltypist_low_Terekhova_et_al %>%
  dplyr::rename(barcodes = `X_Immune_All_Low`)
metadata_Terekhova_et_al <- read_csv(file.path(projectRoot, "Terekhova_et_al_ANALYSIS_HPC", "07.2_Supercentenarian", "raw_data", "Terekhova_et_al", "all_pbmcs", "all_pbmcs", "all_pbmcs_metadata.csv")) 
names(metadata_Terekhova_et_al)[1] <- "barcodes"

celltypist_low_Terekhova_et_al <- celltypist_low_Terekhova_et_al %>% left_join(metadata_Terekhova_et_al %>% dplyr::select(barcodes, Age_group), by = "barcodes")
celltypist_low_Terekhova_et_al <- celltypist_low_Terekhova_et_al %>%
  dplyr::select(majority_voting_Immune_All_Low, Age_group) %>%
  dplyr::rename(sample = Age_group)

# Merge

df <- rbind(df_M116, celltypist_low_Terekhova_et_al)
df <- rbind(df, df_Hashimoto)

# Remove HSC or MPP
table(df$majority_voting_Immune_All_Low)
# HSC or MPP 
# 44
df <- df %>% filter(majority_voting_Immune_All_Low != "HSC or MPP")
df$sample <- factor(df$sample, levels = c("A","B","C","D","E","SC1","SC2","SC3","SC4","SC5","SC6","SC7","M116"))

plot <- df %>% 
  group_by(sample) %>% 
  count(majority_voting_Immune_All_Low) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = sample, y = prop, fill = majority_voting_Immune_All_Low)) +
  geom_col() +
  scale_fill_manual(values = cols_majority_voting_Immune_All_Low) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.7), # Add rectangle border
    panel.grid = element_blank(), # Optional: Remove grid for a cleaner look
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )  +
  labs(x = "", y = "Fraction", fill = "Cell type") 
ggsave(plot, filename = file.path(savePlots, "barplot_celltype.pdf"), width = 6.5, height = 5.5)  


geom_col(color = "black", size = 0.2)

# 1. Get your original order from table() output
celltype_order <- names(table(df$majority_voting_Immune_All_Low))
celltype_order <- celltype_order[celltype_order != "HSC or MPP"]

# 2. Create the plot with perfectly matched orders
plot <- df %>%
  # Convert to factor with original order
  mutate(majority_voting_Immune_All_Low = factor(majority_voting_Immune_All_Low, 
                                                 levels = celltype_order)) %>%
  
  # Calculate proportions
  group_by(sample, majority_voting_Immune_All_Low) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  
  # Visualization
  ggplot(aes(x = sample, y = prop, 
             fill = majority_voting_Immune_All_Low,
             color = "black")) +
  geom_col(position = position_stack(),  # Regular stacking
           linewidth = 0.3) +           # Border thickness
  scale_fill_manual(
    values = cols_majority_voting_Immune_All_Low,
    breaks = celltype_order,            # Original table order
    guide = guide_legend(ncol = 1)      # No reverse here
  ) +
  scale_color_manual(values = "black", guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 0.7),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(x = "", y = "Fraction", fill = "Cell Type")
ggsave(plot, filename = "~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/plots/figures/barplot_celltype.pdf", width = 6.5, height = 5.5)


plot <- DimPlot(data_M116, group.by = "majority_voting_Immune_All_Low", label = F, cols = cols_majority_voting_Immune_All_Low, pt.size = 0.005) +
  labs(title = "Cell type") 
ggsave(plot, filename = "~/Eva/IJC/Supercentenarian/07.2_Supercentenaria/plots/figures/UMAP_M116.pdf", width = 7.5, height = 5.5)

