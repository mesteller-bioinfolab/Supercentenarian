####################################
# DIFFERENTIAL METHYLATION ANALYSIS#
####################################

# Load libraries
library (minfi)
library (limma)
library (IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library (PCAtools)
library (ComplexHeatmap)
library (Rtsne)
library (dplyr)
library (stringr)
library (enrichR)
library (tidyr)
library (data.table)
library (REMP)
library(ggthemes)

WDir<-"./"

# Load annotation
ann850k<-as.data.frame(fread(paste0(WDir,"Annotation.csv")))
ann850k<-ann850k[,c(1:4:7,9)]
head(as.data.frame(ann850k))

# Load Sample Sheet
targets<-readxl::read_xlsx(paste0(WDir,"targets.xlsx"))

# Load processed bVals
bVals<-data.frame(fread("C:/Users/anoguera/PROJECTS/SUPERGRANNY/RESULTS/bVals_M116_SALIVA_LAL.csv"))
rownames(bVals)<-bVals$V1
bVals<-bVals[-1]
bVals

# Find methylation differneces between M116 and CTRL

M116<-as.data.frame(bVals)[1]
colnames(M116)<-"M116"
M116$CpG_Name<-rownames(M116)
dim(M116)
CTRL<-as.data.frame(bVals)[-1]
dim(CTRL)
CTRL<-as.data.frame(rowMeans(CTRL))
colnames(CTRL)<-"CTRL"
CTRL$CpG_Name<-rownames(CTRL)

Mean<-merge(M116,CTRL,"CpG_Name")
Mean$AB_Meth<-Mean$M116-Mean$CTRL
Mean_filt<-Mean[abs(Mean$AB_Meth) >=0.50,]
Mean_filt

DMPs<-merge(as.data.frame(ann850k),Mean_filt, "CpG_Name")
DMPs<-DMPs[order(-abs(DMPs$AB_Meth)),]

DMPs$Methylation_Status<-ifelse(DMPs$AB_Meth > 0, "Hypermethylated","Hypomethylated")

writexl::write_xlsx(DMPs,paste0(WDir,"DMPs_M116vsCTRL.xlsx"))

DMPs

# t-SNE with DMPs

bVals_DMPs<-bVals[rownames(bVals) %in% DMPs$CpG_Name,] 
targets<-targets[match(colnames(bVals_DMPs), targets$Sample_Name),]

# Perform the t-SNE analysis (It is important to transpose the dataframe): CpGs on the column and samples on the row
set.seed(524)
tsne <- Rtsne(t(bVals_DMPs), perplexity=5, check_duplicates = FALSE)

colnames(tsne$Y)<-c("t_SNE1","t_SNE2")

# Keep only the V1 and V2 that are the ones that we are going to represent
tSNE <- tsne$Y %>% 
  as.data.frame() %>%
  mutate(ID=row_number())

# Add a column with ID
targets_PCA<-targets
targets_PCA$ID<-1:nrow(targets_PCA)

# Considering that the row numbers are 1:52, we will add the same numbers to the PCA_metadata
tsne <- tSNE %>% inner_join(targets_PCA , by="ID")

# Represent the analysis into graphics with colours and shapes 
palPA <-colorRamp2::colorRamp2(seq(from=20,to=120,length.out=10),
                               rev(viridis::inferno(10)))

Plot<-tsne %>%ggplot(aes(x = t_SNE1, 
                         y = t_SNE2,
                         fill = Age))+
  geom_point(pch=21, size = 6)+
  scale_fill_gradient(low="darkblue", high = "orange")+
  
  theme(legend.position="bottom")+
  theme_classic()+
  ggthemes::theme_few()

Plot

# Supervised HeatMap

ColAnn<-HeatmapAnnotation(Sample_Type =targets$Project,
                          col = list("Sample_Type" = c("Aging"="#5ADBFF",
                                                       "LookALike"="#FE9000")),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

distance <- dist(t(bVals_DMPs), method="canberra") # calcular las distancias
clust.model <- hclust(distance, method="ward.D2") # clustering

distance <- dist(bVals_DMPs, method="canberra") # calcular las distancias
clust.model_r <- hclust(distance, method="ward.D2") # clustering

Heatmap<-draw(Heatmap(as.matrix(bVals_DMPs), name = "Methylation", 
                      col = gplots::greenred(75), 
                      bottom_annotation = ColAnn,
                      cluster_columns = clust.model,
                      cluster_rows = clust.model_r,
                      show_row_names = FALSE, 
                      show_column_names = FALSE, 
                      show_row_dend = FALSE))
