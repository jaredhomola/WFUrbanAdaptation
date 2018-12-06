###########################################################################
#####   Calculate population-specific measures of genetic diversity   #####
###########################################################################

#### Load required libraries
library(tidyverse)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

#### Set up tibble
drops <- c("Batch.ID", "Chr", "Col", "N", "Smoothed.Fis", "Smoothed.Fis.P.value", "Smoothed.Pi",
           "Smoothed.Pi.P.value", "Obs.Hom", "P.Nuc", "Q.Nuc", "BP", "Fis", "Exp.Hom", "Exp.Het")
dat <- sumStats.noOutliers[ , !(names(sumStats.noOutliers) %in% drops)]

#### Measures for diversity table #####
## Subset data
dat.5 <- subset(dat, Pop.ID =="5")
dat.8 <- subset(dat, Pop.ID =="8")
dat.4 <- subset(dat, Pop.ID =="4")
dat.6 <- subset(dat, Pop.ID =="6")
dat.1 <- subset(dat, Pop.ID =="1")
dat.2 <- subset(dat, Pop.ID =="2")
dat.7 <- subset(dat, Pop.ID =="7")
dat.3 <- subset(dat, Pop.ID =="3")

## Number of sites
sum(dat$Pop.ID == "5")
sum(dat$Pop.ID == "8")
sum(dat$Pop.ID == "4")
sum(dat$Pop.ID == "6")
sum(dat$Pop.ID == "1")
sum(dat$Pop.ID == "2")
sum(dat$Pop.ID == "7")
sum(dat$Pop.ID == "3")

## Percent polymorphic sites
nrow(subset(dat, Pop.ID =="5" & Obs.Het > 0)) / sum(dat$Pop.ID == "5")
nrow(subset(dat, Pop.ID =="8" & Obs.Het > 0)) / sum(dat$Pop.ID == "8")
nrow(subset(dat, Pop.ID =="4" & Obs.Het > 0)) / sum(dat$Pop.ID == "4")
nrow(subset(dat, Pop.ID =="6" & Obs.Het > 0)) / sum(dat$Pop.ID == "6")
nrow(subset(dat, Pop.ID =="1" & Obs.Het > 0)) / sum(dat$Pop.ID == "1")
nrow(subset(dat, Pop.ID =="2" & Obs.Het > 0)) / sum(dat$Pop.ID == "2")
nrow(subset(dat, Pop.ID =="7" & Obs.Het > 0)) / sum(dat$Pop.ID == "7")
nrow(subset(dat, Pop.ID =="3" & Obs.Het > 0)) / sum(dat$Pop.ID == "3")

## Percent of private alleles
nrow(subset(dat, Pop.ID =="5" & Private > 0)) / sum(dat$Pop.ID == "5")
nrow(subset(dat, Pop.ID =="8" & Private > 0)) / sum(dat$Pop.ID == "8")
nrow(subset(dat, Pop.ID =="4" & Private > 0)) / sum(dat$Pop.ID == "4")
nrow(subset(dat, Pop.ID =="6" & Private > 0)) / sum(dat$Pop.ID == "6")
nrow(subset(dat, Pop.ID =="1" & Private > 0)) / sum(dat$Pop.ID == "1")
nrow(subset(dat, Pop.ID =="2" & Private > 0)) / sum(dat$Pop.ID == "2")
nrow(subset(dat, Pop.ID =="7" & Private > 0)) / sum(dat$Pop.ID == "7")
nrow(subset(dat, Pop.ID =="3" & Private > 0)) / sum(dat$Pop.ID == "3")

## Average frequency of the major allele (P)
mean(dat.5$P)
mean(dat.8$P)
mean(dat.4$P)
mean(dat.6$P)
mean(dat.1$P)
mean(dat.2$P)
mean(dat.7$P)
mean(dat.3$P)

## Average observed heterozygosity (Obs.Het)
mean(dat.5$Obs.Het)
mean(dat.8$Obs.Het)
mean(dat.4$Obs.Het)
mean(dat.6$Obs.Het)
mean(dat.1$Obs.Het)
mean(dat.2$Obs.Het)
mean(dat.7$Obs.Het)
mean(dat.3$Obs.Het)

## Average nucleotide diversity (Pi)
mean(dat.5$Pi)
mean(dat.8$Pi)
mean(dat.4$Pi)
mean(dat.6$Pi)
mean(dat.1$Pi)
mean(dat.2$Pi)
mean(dat.7$Pi)
mean(dat.3$Pi)
