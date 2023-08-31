library(corrplot)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(patchwork)
library(ggplot2)

####oxobased
Correlation_data_oxobased <- read.csv("Files/MS1/Oxobased_PseudoMALDI.csv", header = TRUE, sep = ";") %>% select( - mass)
Correlation_data_oxobased_reduced <- Correlation_data_oxobased  %>% select( - BAT1, - BAT2, - bladder1, - bladder2, - stomach1, - stomach2, -mammaryGland1, - mammaryGland2)
Correlation_oxobased_reduced <- cor(Correlation_data_oxobased_reduced,method = "pearson")

Correlationplot_oxobased <- pheatmap(Correlation_oxobased_reduced, fontsize = 10, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
         clustering_method = "ward.D")

png("Correlation_oxobased.png", width     = 8,
    height    = 7,
    units     = "in",
    res       = 1200,
    pointsize = 4)
Correlationplot_oxobased
dev.off()
####targeted extended
load(file ="./RDA/PseudoMALDI_targeted_extended.rda") 
Correlation_data_targeted_reduced <- PseudoMALDI_targeted_extended %>% select( - mass, - BAT1, - BAT2, - Bladder1,- Bladder2, - Stomach1, - Stomach2, - MammaryGland1, - MammarygGland2)
Correlation_targeted_reduced <- cor(Correlation_data_targeted_reduced,method = "pearson")
corrplot(Correlation_targeted, method = "color",order = "hclust",is.corr = TRUE,  tl.cex = 0.6,col.lim = c(0,1),tl.col = "black")
Correlationplot_targeted <- pheatmap(Correlation_targeted_reduced, fontsize = 10, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
         clustering_method = "ward.D2")

png("Correlation_targeted.png", width     = 7,
    height    = ,
    units     = "in",
    res       = 1200,
    pointsize = 4)
Correlationplot_targeted
dev.off()

#####SNOG
Correlation_data_SNOG8 <- read.csv("Files/MS1/SNOG8_SetforClustering.csv", header = TRUE, sep = ";") %>% select( - mass,
                                                                                                                 - colon1,
                                                                                                                 - esophagus1,
                                                                                                                 - esophagus2,
                                                                                                                 - muscle1,
                                                                                                                 - muscle2,
                                                                                                                 - muscle3re,
                                                                                                                 - bladder1,
                                                                                                                 - bladder2,
                                                                                                                 - BAT1,
                                                                                                                 - BAT2,
                                                                                                                 - stomach1,
                                                                                                                 - stomach2,
                                                                                                                 - mammaryGland1,
                                                                                                                 - mammaryGland2)
                                                                                                                 
                                                                                                                

Cor_SNOG8 <- cor(Correlation_data_SNOG8, method = "pearson")
Correlation_SNOG6 <- pheatmap(Cor_SNOG8, fontsize = 10, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
         clustering_method = "ward.D")

png("Correlation_SNOG6.png", width     = 8,
    height    = 7,
    units     = "in",
    res       = 1200,
    pointsize = 4)
Correlation_SNOG6
dev.off()
