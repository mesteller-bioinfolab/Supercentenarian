#################
# METHYLATION QC#
#################

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

## HUMAN ANNOTATION

# Load annotation
ann850k<-as.data.frame(fread(paste0(WDir,"Annotation.csv")))
ann850k<-ann850k[,c(1:4:7,9)]
head(as.data.frame(ann850k))

# Load Sample Sheet
targets<-readxl::read_xlsx(paste0(WDir,"targets.xlsx"))

rgSet <- read.metharray.exp(targets=targets, force = TRUE)
sampleNames(rgSet) <- targets[["Sample_Name"]]

detP <- detectionP(rgSet)
keep <- colMeans(detP) < 0.05
table(keep)

keep<-colnames(rgSet)
rgSet <- rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]
mSetSq <- preprocessNoob(rgSet)
mSetSq

detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.01) >=0.9 * ncol(mSetSq)
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt <- mapToGenome(mSetSqFlt)
keep <- !(featureNames(mSetSqFlt) %in% ann850k$CpG_Name[ann850k$CpG_chrm %in% c("chrX","chrY")])
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt <- maxprobes::dropXreactiveLoci(mSetSqFlt)
mSetSqFlt

bVals <- getBeta(mSetSqFlt)
detP <- detP[match(rownames(bVals),rownames(detP)),]
failed_cpgs<- NULL
failed_cpgs <- which(detP > 0.01, arr.ind = TRUE)
bVals[failed_cpgs] <- NA
bVals<-t(apply(bVals,1,function(x) replace(x, is.na(x), median(x,na.rm=TRUE))))
dim(bVals)

head(as.data.frame(bVals))

write.csv(bVals,paste0(WDir,"bVals.csv"),row.names = TRUE)
