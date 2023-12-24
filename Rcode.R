## Diversity and structuring process
## of understory bird assemblages in Balbina
##
## Ivana Cardoso - ivanawaters@gmail.com
##
## Last modification: December 24, 2023

# Set working directory
setwd("C:/Users/Ivana/OneDrive/PhD_INPA/1.Diversity_question/Analises/birds_diversity_Balbina")

# Load packages
library(ape)

# Import data
tree = ape::read.nexus(unzip("1000_trees_vertnet.zip", "output.nex"))
