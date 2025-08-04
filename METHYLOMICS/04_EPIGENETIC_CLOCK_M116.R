#############################
# EPIGENETIC CLOCK ANALYSIS #
#############################

# Load libraries
library (tibble)
library (methylclock)
library(preprocessCore)

# Load targets
targets_M116 <- read.metharray.sheet("~/M116/idats_torun", pattern="30.09.24_M116.csv")

# Load bVals 
bVals_M116_clock<-as.data.frame(data.table::fread("./bVals.csv"))

bVals_M116_clock <- bVals_M116
bVals_M116_clock <- as.data.frame(bVals_M116_clock)

bVals_M116_clock <- tibble::rownames_to_column(bVals_M116_clock, "CpGNAme")

cpgs.missing <- checkClocks(bVals_M116_clock)
cpgs.missing.GA <- checkClocksGA(bVals_M116_clock)

DNAmAge(bVals_M116_clock)
age.M116 <- DNAmAge(bVals_M116_clock)