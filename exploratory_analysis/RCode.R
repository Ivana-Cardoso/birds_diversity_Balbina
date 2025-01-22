# EXPLORATORY ANALYSIS
#
# Diversity and structuring process
# of understory bird assemblages in Balbina
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: January 21, 2025

rm(list = ls())

# Load packages
library(ggplot2)
library(ggpubr)
library(iNEXT)

# Set work directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina")

# Import data
comm <- read.csv("balbina_comm.csv", row.names = 1)
traits <- read.csv("balbina_understory_bird_traits.csv", row.names = 2)
traits <- traits[,-1]

# Heatmap of the community data to assess:
# 1. Islands with low number of species and individuals
# 2. Which species are rare
heatmap(as.matrix(t(comm)), scale = "row", col= viridis(10))

# Check species sampling coverage
iNEXT_result <- iNEXT::iNEXT(t(comm), q = 0, datatype = "abundance", size = NULL)
iNEXT_result$DataInfo
# Adeus, Bacaba, Panema, Torem and Xibe had Sample Coverage below 0.5. This is because species that occurred only a few times, or even just once, were recorded at these locations. Examples include Geotrygon montana (Torem), Ceratopipra erythrocephala (Bacaba and Xibe), Glaucis hirsutus (Xibe), Schiffornis turdina (Adeus), Momotus momota (Bacaba), and Selenidera piperivora (Panema).

coverage <- iNEXT_result$iNextEst$coverage_based
coverage <- subset(coverage, Method == "Rarefaction")

coverage_plot <- 
  ggplot(coverage, aes(x = m, y = SC, color = Assemblage)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Sample Coverage Curves by Assemblage",
    x = "Number of individuals",
    y = "Sample Coverage (SC)") +
  facet_wrap(~Assemblage) +  
  theme_minimal() +
  theme(legend.position = "none")


