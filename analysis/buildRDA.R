############################################
### Creating the .rda file from raw data ###
############################################

library(vcfR)

getwd()
allLoci <- read.vcfR("./extData/allLoci.vcf")
sumStats.noOutliers <- read.delim("./extData/batch_1.sumstats.tsv", sep="\t", header = TRUE)
sumStats.pops <- read.delim("./extData/batch_1.sumstats.pops.csv", sep="\t", header = FALSE)
envData <- read.csv("./extData/envData.csv", header = TRUE)
outlierLoci <- read.csv("./extData/outlierLoci.csv", header = TRUE)

save(allLoci, sumStats.noOutliers, sumStats.pops, envData, outlierLoci, file = "./data/WFUrbanAdaptation.rda")
