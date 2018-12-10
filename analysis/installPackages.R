#########################################################
#####         Install packages needed for           #####
#####          WFUrbanAdaptation analyses           #####
#########################################################

list.of.packages <- c("vegan", "tidyverse",
                      "adegenet", "vcfR",
                      "lme4", "stringr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

#### Install qvalue from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")


