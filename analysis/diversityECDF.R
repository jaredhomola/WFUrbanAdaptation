######################################################################
#####  Create empirical cumulative density function to compare   #####
#####    genetic diversity statistics among environment types    #####
######################################################################
##########       Decoding sumstats file (Pop.ID = Pop)      ##########
##########               5 = Rep 1 - Rural                  ##########
##########               8 = Rep 2 - Rural                  ##########
##########               4 = Rep 3 - Rural                  ##########
##########               6 = Rep 4 - Rural                  ##########
##########               1 = Rep 1 - Urban                  ##########
##########               2 = Rep 2 - Urban                  ##########
##########               7 = Rep 3 - Urban                  ##########
##########               3 = Rep 4 - Urban                  ##########
######################################################################

#### Load required libraries
library(tidyverse)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

#### Set up tibble
drops <- c("Batch.ID", "Chr", "Col", "N", "Smoothed.Fis", "Smoothed.Fis.P.value", "Smoothed.Pi",
           "Smoothed.Pi.P.value", "Obs.Hom", "P.Nuc", "Q.Nuc", "BP", "Fis", "Exp.Hom", "Exp.Het")
dat <- sumStats.noOutliers[ , !(names(sumStats.noOutliers) %in% drops)]

dat$envType <- NA
dat$envType[dat$Pop.ID == 1] <- "Urban"
dat$envType[dat$Pop.ID == 2] <- "Urban"
dat$envType[dat$Pop.ID == 3] <- "Urban"
dat$envType[dat$Pop.ID == 4] <- "Rural"
dat$envType[dat$Pop.ID == 5] <- "Rural"
dat$envType[dat$Pop.ID == 6] <- "Rural"
dat$envType[dat$Pop.ID == 7] <- "Urban"
dat$envType[dat$Pop.ID == 8] <- "Rural"

#### Subset data and compare ####
dat.urban <- subset(dat, envType == "Urban")
dat.rural <- subset(dat, envType == "Rural")

## Observed heterozygosity
P.urban = ecdf(dat.urban$Obs.Het)
plot(P.urban, ylim=c(0.4, 1.0))
P.rural = ecdf(dat.rural$Obs.Het)
lines(P.rural, col = "gray70")

ks.test(dat.rural$Obs.Het, dat.urban$Obs.Het)
mean(dat.urban$Obs.Het)
mean(dat.rural$Obs.Het)

## Nucleotide diversity
P.urban = ecdf(dat.urban$Pi)
plot(P.urban, ylim=c(0.4, 1.0))
P.rural = ecdf(dat.rural$Pi)
lines(P.rural, col = "gray70")

ks.test(dat.rural$Pi, dat.urban$Pi)
mean(dat.urban$Pi)
mean(dat.rural$Pi)
