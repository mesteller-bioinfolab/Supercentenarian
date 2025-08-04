################################
# Repetitive Elements Analysis #
################################

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
library (ggthemes)

# Functions
EnricheR_Table<-function(Genes_Vector){
  Enriched_Results<-enrichr(as.vector(Genes_Vector),dbs_selection)
  
  Enriched_Table_Results<-data.frame()
  
  for (i in 1:length(dbs_selection)) {
    Database<-dbs_selection[i]
    Enriched_Table<-Enriched_Results[[Database]]
    Enriched_Table$Database<-Database
    Enriched_Table<-Enriched_Table[,-c(5:8)]
    Enriched_Table$Obs_Genes<-as.numeric(gsub("/..*","",Enriched_Table$Overlap))
    Enriched_Table$Set_Genes<-as.numeric(gsub("..*/","",Enriched_Table$Overlap))
    Enriched_Table<-Enriched_Table[Enriched_Table$Adjusted.P.value<=0.05,]
    if(!nrow(Enriched_Table)==0)
      Enriched_Table_Results<-rbind(Enriched_Table_Results,Enriched_Table)
    else{
      next 
    }
  }
  
  if(nrow(Enriched_Table_Results)>0){
    Enriched_Table_Results$Gene_Ratio<-as.numeric(Enriched_Table_Results$Obs_Genes/Enriched_Table_Results$Set_Genes)
    Enriched_Table_Results<-Enriched_Table_Results[,c(6,1:5,9)]
    Enriched_Table_Results<-Enriched_Table_Results[order(-Enriched_Table_Results$Gene_Ratio),]
    return(Enriched_Table_Results)
  } else {
    return("No differential pathways found")
  }
}

# Function to analyze Repetitive Elements with methylation microarray data

Analysis_RE<-function(remp.res,targets){
  
  list_Report<-list()
  
  # TRIM
  print("Trimming ...")
  remp.res_Trim <- rempTrim(remp.res, threshold = 1.7, missingRate = 0.2)
  
  # Aggregate
  print("Aggregating ...")
  remp.res_Agg <- rempAggregate(remp.res_Trim, NCpG = 2)
  
  # Extract B-values
  print("Extracting B-values ...")
  RE_Meth<-as.data.frame(rempB(remp.res_Agg))
  list_Report[["Bvals"]]<-RE_Meth
  
  # Extracting the Annotation of RE
  print("Extracting Annotation ...")
  Annotation_RE<-decodeAnnot(remp.res_Agg)
  Annotation_RE<-as.data.frame(rempAnnot(Annotation_RE))
  list_Report[["Annotation"]]<-Annotation_RE
  
  
  
  return(list_Report)
  
  print("ALL DONE ENJOY YOUR RESULTS!!!")
  
}

# Load annotation
ann850k<-as.data.frame(fread("./EPICv1_UCSC_Annotation_Full.csv"))

# Open Sample Sheet
targets<-readxl::read_xlsx("./Sample_Sheet.xlsx")

# Open bVals data
bVals_RE <- as.data.frame(data.table::fread("./bVals_RE.csv"))
rownames(bVals_RE)<-bVals_RE$V1
bVals_RE<-bVals_RE[-1]

remparcel <- initREMP(arrayType = "EPIC",
                      REtype = "Alu", # Change the RE type for each RE of interest
                      annotation.source = "UCSC",
                      genome = "hg38")

remp.res <- remp(bVals_RE,
                 REtype = 'Alu',
                 parcel = remparcel,
                 seed = 564)

save(remp.res, file = "./RESULTS/RE_Prediction_ALL.RData")


load("./CLUSTER_FILES/RE_PREDICTION_ARRAY/RE_ALU_Prediction_ALL.RData")
remp.res_ALU <- remp.res
RE_ALU_ANALYSIS<-Analysis_RE(remp.res_ALU,targets)

RE_ALU_ANALYSIS$Annotation

RE_Plot_Individual<-function(Obj){
  
  bVals_RE<-Obj$Bvals
  bVals_RE<-na.omit(bVals_RE)
  
  df<-data.frame()
  
  for (i in 1:ncol(bVals_RE)) {
    
    Pat<-colnames(bVals_RE)[i]
    temp<-bVals_RE[i]
    colnames(temp)<-"Meth_Val"
    temp$ER_ID<-rownames(temp)
    temp$Sample_Name<-Pat
    rownames(temp)<-NULL
    
    df<-rbind(df,temp)
  }
  
  df$Condition<-ifelse(df$Sample_Name == "M116_Saliva","M116","Control")
  df<-merge(df,targets[c(1,12)],"Sample_Name")
  
  
  df2 <- df %>%
    mutate(Sample_Name = reorder(Sample_Name, Meth_Val)) %>%
    arrange(Sample_Name)  # Arrange the dataframe to match the x-axis order
  
  temp<-targets[c(1,12)]
  temp<-temp[match(unique(df2$Sample_Name),temp$Sample_Name),]
  
  Dic<-setNames(temp$Age,temp$Sample_Name)
  
  ggplot_GLobal <- ggplot(df2, aes(x=reorder(Sample_Name,Meth_Val), y=Meth_Val, fill=Condition)) + 
    scale_fill_manual(values = c("M116" = "#FC8D59",
                                 "Control" = "#91BFDB")) + 
    ylim(0,1)+
    geom_boxplot() +
    theme_few() +
    theme(#axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+ 
    scale_x_discrete(labels= Dic)
  
  
  return(ggplot_GLobal)
}

ALU_Plot<-RE_Plot_Individual(RE_ALU_ANALYSIS)
Graphic<-ggplot_GLobal

load("./CLUSTER_FILES/RE_PREDICTION_ARRAY/RE_L1_Prediction_ALL.RData")
RE_L1_ANALYSIS<-Analysis_RE(remp.res_L1,targets)

L1_Plot<-RE_Plot_Individual(RE_L1_ANALYSIS)
L1_Plot

load("./CLUSTER_FILES/RE_PREDICTION_ARRAY/RE_ERV_Prediction_ALL.RData")
RE_ERV_ANALYSIS<-Analysis_RE(remp.res_ERV,targets)

ERV_Plot<-RE_Plot_Individual(RE_ERV_ANALYSIS)
