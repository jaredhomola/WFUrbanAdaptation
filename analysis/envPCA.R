######################################################################
#####  Principal components analysis of site environmental data  #####
######################################################################

#### Load required libraries
library(vegan)
library(ggplot2)
library(WFUrbanAdaptation)
data(WFUrbanAdaptation)

##### Log transform environmental data
envData.log <- log(envData[4:7]) ## Log transform numerical data
envData.log <- cbind(envData$Urban.rural, envData.log) ## Add environment type (Urban.rural) column
names(envData.log) <- c("envType", "Imperv", "canopyCover", "distNearestRd", "percentDevelopment")

##### Perform PCA
pca <- prcomp(envData.log[2:5], center = TRUE, scale = TRUE)

## Perform ADONIS
adonis(envData.log[2:5] ~ envData.log$envType, method='eu', permutations = 99999)

## Review scree plot and summary
plot(pca, type = "l")
summary(pca)

## Plot PCA
pca.df<-data.frame(pca$x, envType = envData.log$envType)

pca.plot <- ggplot(pca.df, aes(x = PC1, y = PC2, group = envType)) +
  geom_point(size = 8, aes(shape = envType)) +
  labs(shape = "Site type") +
  xlab("PC1 (88.89%)") +
  ylab("PC2 (7.08%)") +
  scale_color_manual(values = c("#0033cc", "#ff0000")) +
  theme(
    text = element_text(size = 22),
    legend.position = c(0.85, 0.15),
    legend.key = element_blank(),
    legend.background = element_rect(
      color = "black",
      fill = "white",
      size = 1,
      linetype = "solid"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
