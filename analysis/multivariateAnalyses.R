################################################
####        Multivariate analyses including ####
####            1. DAPC of all loci         ####
####            2. DAPC of outliers         ####
####            3. PCoA of outliers         ####
################################################
###   List of outliers provided in .rda,     ###
### but created using outlierAnalysis script ###
################################################

#### Load required libraries
library(tidyverse)
library(vcfR)
library(adegenet)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

#### Load data, subset outliers, establish site categories (population, environment type, and replicate)
genind.dat <- vcfR2genind(allLoci)

outliers.dat <- genind.dat[loc=outlierLoci$locusName]

pops <- c(rep("URB-1",12),rep("URB-2",12),rep("URB-4",12),rep("RUR-3",12),rep("RUR-1",12),rep("RUR-4",12),rep("URB-3",12),rep("RUR-2",12))
envType <- c(rep("Urban",12),rep("Urban",12),rep("Urban",12),rep("Rural",12),rep("Rural",12),rep("Rural",12),rep("Urban",12),rep("Rural",12))
rep <- c(rep("Rep 1",12),rep("Rep 2",12),rep("Rep 4",12),rep("Rep 3",12),rep("Rep 1",12),rep("Rep 4",12),rep("Rep 3",12),rep("Rep 2",12))



##### DAPC with all loci
## Factor = Replicate
genind.dat@pop <- as.factor(rep)
dapc.all <- dapc(genind.dat, var.contrib = TRUE, scale = FALSE, n.pca = 25, n.da = nPop(genind.dat) - 1)
scatter(dapc.all, cell = 0, pch = 19, cstar = 0, lwd = 2, cex = 2, lty = 2)
## Reassignment probability
summary(dapc.all)$assign.per.pop*100

## Factor = Population
genind.dat@pop <- as.factor(pops)
dapc.all <- dapc(genind.dat, var.contrib = TRUE, scale = FALSE, n.pca = 25, n.da = nPop(genind.dat) - 1)
scatter(dapc.all, cell = 0, pch = 19, cstar = 0, lwd = 2, cex = 2, lty = 2)
## Reassignment probability
summary(dapc.all)$assign.per.pop*100



##### DAPC with outliers
## Factor = Environment Type
outliers.dat@pop <- as.factor(envType)
dapc.ol <- dapc(outliers.dat, var.contrib = TRUE, scale = FALSE, n.pca = 25, n.da = nPop(outliers.dat) - 1)
scatter(dapc.ol, cell = 0, pch = 19, cstar = 0, lwd = 2, cex = 2, lty = 2)

## Loading plot
loadingplot(dapc.ol$var.contr, axis = 1, thres = 0.04,
          lab.jitter = 1, srt=90, adj = 0)

## Reassignment probability
summary(dapc.ol)$assign.per.pop*100



##### PCoA with outliers
X <- tab(outliers.dat, freq=TRUE, NA.method="mean")
pcoa.ol <- dudi.pco(dist(X), scannf=FALSE, nf=3)
col <- c("red", "blue")
s.class(pcoa.ol$li, pop(outliers.dat), xax=1, yax=2, col=transp(col, 0.6),
        axesell=FALSE, cellipse = 0, cstar=0, cpoint=3, pch=16:17, grid=FALSE)

