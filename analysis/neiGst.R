###################################################
####       Population pairwise (Nei 1973)      ####
###################################################

##### Load required libraries
library(hierfstat)
library(vcfR)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

##### Convert vcfR object to genind object and add @pop vector
dat.genind <- vcfR2genind(allLoci)
pops <- c(rep("URB-1",12),rep("URB-2",12),rep("URB-4",12),rep("RUR-3",12),rep("RUR-1",12),rep("RUR-4",12),rep("URB-3",12),rep("RUR-2",12))
dat.genind@pop <- as.factor(pops)

##### Calculate pairwise population Fst
pairwise.fsts <- pairwise.fst(dat.genind, res.type = "matrix")

