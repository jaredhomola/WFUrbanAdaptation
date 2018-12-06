###################################################
####  Outlier analysis to determine systematic ####
####    preference for certain genotypes in    ####
####        urban or rural environments        ####
###################################################

##### Load required libraries
library(tidyverse)
library(vcfR)
library(lme4)
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

##### Cacluate Nei's Fst for each locus within each replicate ####
rep1 <- allLoci[, c(1:13, 49:60)]
rep1.pop <- as.factor(rep(c('BAN7', 'ORO1'), each = 12))
rep1.diff <- genetic_diff(rep1, rep1.pop, method = "nei")

rep2 <- allLoci[, c(1, 14:25, 86:97)]
rep2.pop <- as.factor(rep(c('BRU3', 'BRU2'), each = 12))
rep2.diff <- genetic_diff(rep2, rep2.pop, method = "nei")

rep3 <- allLoci[, c(1, 38:49, 74:85)]
rep3.pop <- as.factor(rep(c('FRE3', 'YAR4'), each = 12))
rep3.diff <- genetic_diff(rep3, rep3.pop, method = "nei")

rep4 <- allLoci[, c(1, 26:37, 62:73)]
rep4.pop <- as.factor(rep(c('WEL4', 'TAT1'), each = 12))
rep4.diff <- genetic_diff(rep4, rep4.pop, method = "nei")

loci <- all.data$dat %>% distinct(ID)
loci <- loci$ID

##### Place all Gst values into a new df
all.gst <- as.data.frame(cbind(as.factor(1:8344), loci, as.numeric(rep1.diff$Gst), rep2.diff$Gst, rep3.diff$Gst, rep4.diff$Gst))
all.gst <- setNames(all.gst, c("locNumber", "locusName", "Rep1","Rep2","Rep3", "Rep4"))
all.gst[all.gst =="NaN"] <- 0 ## NaNs can't be calculated because they're all fixed. So, Fst = 0

##### Move data into tibble
fst.df <- as_tibble(all.gst) %>%
          mutate(Rep1 = as.numeric(levels(Rep1))[Rep1],
                 Rep2 = as.numeric(levels(Rep2))[Rep2],
                 Rep3 = as.numeric(levels(Rep3))[Rep3],
                 Rep4 = as.numeric(levels(Rep4))[Rep4])

###### Assign each locus a percentile based on the observed range of Gst within each replicate ####
fst.df <- fst.df %>%
  gather('Rep1', 'Rep2', 'Rep3', 'Rep4', key = "Replicate", value = "Gst") %>%
  mutate(Replicate = factor(Replicate), locNumber = factor(locNumber)) %>%
  group_by(Replicate) %>%
  mutate(Percentile.Gst = Gst / max(Gst))

##### Create new tibble that includes the mean GST percentile and its min
fst.df.percent <- fst.df %>%
  group_by(locusName, locNumber) %>%
  summarize(mean.Gst = mean(Gst), mean.Percentile.Gst = (mean(Percentile.Gst)) * 100,
            min.Gst = min(Gst), se.Gst = (sd(Gst) / sqrt(4))*100) %>%
  mutate(outlier = ifelse(mean.Percentile.Gst-se.Gst > 20, "Outlier" , "Not Outlier" ))

fst.df.percent$rank <- rank(fst.df.percent$mean.Percentile.Gst, ties.method = "random")

##### Isolate list of outlier loci
outlierLoci <- subset(fst.df.percent, outlier =="Outlier")
write.csv(outlierLoci, "./extData/outlierLoci.csv")

##### Plot percentiles across loci including error bars
ggplot(fst.df.percent, aes(locNumber, mean.Percentile.Gst, color=factor(outlier))) +
  geom_errorbar(aes(ymax=mean.Percentile.Gst+se.Gst,ymin=mean.Percentile.Gst-se.Gst), width=0.5) +
  geom_point(size=2, shape=16) +
  scale_colour_manual(values=c("gray85", "black")) +
  ylab(expression(paste('Mean F'['ST'], ' Percentile'))) +
  expand_limits(x = c(0, 8400), y = c(0,50)) +
  scale_x_discrete(expression("Locus Number"), breaks=c(0,2000,4000,6000,8000)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position="none", axis.line = element_line(colour = "black"))

##### Plot in order of mean Gst
ggplot(fst.df.percent, aes(rank, mean.Percentile.Gst, color=factor(outlier))) +
  geom_point(size=2, shape=16) +
  scale_colour_manual(values=c("gray85", "black")) +
  ylab(expression(paste('Mean F'['ST'], ' Percentile'))) +
  scale_x_continuous(expression(paste('F'['ST'], ' Percentile Rank')), breaks=c(0,2000,4000,6000,8000)) +
  expand_limits(x = c(0, 8400), y = c(0,50)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position="none", axis.line = element_line(colour = "black"), legend.title=element_blank())


##### Generate genotype frequency plots
## Create tibble
drop.cols <- c("REF", "ALT", "POS", "NS", "AF", "gt_GT", "recode.GT")
df.plot <- df %>%
  select(-one_of(drop.cols)) %>%
  group_by(ID, envType, rep, gt_GT_alleles) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n))

df.plot2 <- within(df.plot,
                   envType <- factor(envType,
                                     levels = names(sort(
                                       table(envType),
                                       decreasing = TRUE
                                     ))))

## Plots created locus-by-locus.
locus <- "71089_22" ## Specify locus of interest
df.plot2 %>%
  filter(ID == locus) %>%
  ggplot(aes(x = envType, y = freq, fill = gt_GT_alleles)) +
  geom_boxplot() +
  ylim(0.0,1.0) +
  ylab("Genotype frequency") + xlab("") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title=element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 22))
