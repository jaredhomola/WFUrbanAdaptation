#####################################################
####  Genotype-environment association analysis  ####
####           using logistic regression         ####
#####################################################

#### Load required libraries
library(tidyverse)
library(vcfR)
library(lme4)
library(qvalue)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

##### Set up locus tibble
all.data <- vcfR2tidy(allLoci,
                      single_frame = TRUE,
                      format_fields = c("GT"))

drop.cols <- c("CHROM", "QUAL", "FILTER")

df <- all.data$dat %>%
  separate(Indiv, c("Site1", "Site2", "Ind")) %>%
  unite("Site", c("Site1", "Site2")) %>%
  select(-one_of(drop.cols))

temp <- recode_factor(df$Site,
                      "BAN_7" = "Urban",
                      "ORO_1" = "Rural",
                      "BRU_3" = "Urban",
                      "BRU_2" = "Rural",
                      "YAR_4" = "Urban",
                      "FRE_3" = "Rural",
                      "WEL_4" = "Urban",
                      "TAT_1" = "Rural")

temp2 <- recode_factor(df$Site,
                       "BAN_7" = "1",
                       "ORO_1" = "1",
                       "BRU_3" = "2",
                       "BRU_2" = "2",
                       "YAR_4" = "3",
                       "FRE_3" = "3",
                       "WEL_4" = "4",
                       "TAT_1" = "4")

temp3 <- recode_factor(df$gt_GT,
                       "0/0" = "0",
                       "0/1" = "1",
                       "1/1" = "2")

df <- df %>% mutate(envType = temp,
                    rep = temp2,
                    recode.GT = temp3) %>%
  filter(!is.na(recode.GT)) %>%
  filter(ID != "40030_41")

##### Describe the regression function
locus.logReg<-function(df) {
  glm(envType~recode.GT + rep, data=df, family=binomial)
}

##### Loop though all loci with the regression function
result.p <- list()
sample <- list()
j = 0

for (i in unique(df$ID)) {
  test.df <- subset(df, ID == i)
  fit <- locus.logReg(test.df)
  result.p[[length(result.p)+1]] = coef(summary(fit))[2,4]
  sample[[length(sample)+1]] = i
  j = j + 1
  print(j)
}


###### Organize results and do FDR
results.p.df <- do.call("rbind", lapply(result.p, as.data.frame))
sample.df <- do.call("rbind", lapply(sample, as.data.frame))
all.df <- cbind(results.p.df, sample.df)
colnames(all.df)[1] <- "p"
colnames(all.df)[2] <- "loc"

## Histograph of p values
p.list <- as.numeric(all.df$p)
hist(all.df$p, xlab = "p value")

## Identify outlier loci from outlierAnalysis.R with p < 0.05
t <- subset(all.df, loc %in% outlierLoci$locusName )
pSig <- subset(t, p < 0.05)

## Assess using false discovery rate of 0.05
alpha <- 0.05
qvals <- qvalue(p = p.list)$qvalues
outliers <- which(qvals<alpha) ## No loci significant at q = 0.05

