# Diversity and structuring process
# of understory bird assemblages in Balbina
#
# EXPLORATORY ANALYSES
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: January 28, 2025

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

# Plot the distribution of collected continuous traits to check for outliers
hist(traits$Peso_real)
hist(traits$CT)
hist(traits$Asa)
hist(traits$Remiges_primarias)
hist(traits$Remiges_secundarias)
hist(traits$Cauda)
hist(traits$Bico_comprimento)
hist(traits$Bico_largura)
hist(traits$Bico_altura)
hist(traits$Tarso_comprimento)
hist(traits$Tarso_largura)

# Plot categorical traits to check the balance of levels
barplot(table(traits$Habitat), xlab = "Type of Habitat")
barplot(table(traits$Trophic.Niche), xlab = "Trophic Niche")
# The majority of species are forest invertivores



### Species richness X area complete data (2015-16 and 2023-24)
# Load required packages
library(ggplot2)
library(ggpubr)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/Data")

# Import data
comm_2023_24 <- read.csv("balbina_comm_islands_2023-24.csv", row.names = 1)
env_2023_24 <- read.csv("balbina_environmental_islands_2022.csv", row.names = 1)

comm_2015_16 <- read.csv("balbina_comm_islands_2015-16.csv", row.names = 1)
env_2015 <- read.csv("balbina_environmental_islands_2015.csv", row.names = 1)


#### TAXONOMIC DIVERSITY - Species Richness ####
#### Data of 2023 and 2024

# Calculate species richness
env_2023_24$richness = rowSums(ifelse(comm_2023_24 > 0, 1, 0))

# Log-transform the variables (log10) and
# Calculate simple linear regression model (sr ~ area)
mod1 = lm(log10(env_2023_24$richness) ~ log10(env_2023_24$area.ha.2022))
summary(mod1)

# Test the normality of the residuals from the linear regression model 
# to ensure the validity of the statistical inference
hist(mod1$residuals)
shapiro.test(mod1$residuals) # It is normal (p=0.1503)

# Plot the graph showing the relationship between number of bird 
# species and island area
SR_plot = 
  ggplot() +
  
  labs( x = "Island area (ha)",
        y = "Number of bird species (n)") +
  
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_log10(limits = c(1,50),
                breaks = c(1, 3, 10, 30 , 50),
                labels = c("1", "3", "10", "30", "50")) +
  annotation_logticks() +
  
  geom_smooth(data = env_2023_24,
              mapping = aes(x = area.ha.2022, y = richness),
              method = lm,
              size = 1.5,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = env_2023_24, aes(x = area.ha.2022,
                                     y = richness)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = env_2023_24, aes(x = area.ha.2022,
                                     y = richness)) +
  
  theme_bw(base_size = 20) +
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm")) +
  
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = c(0.95, 0.05),
        legend.justification = c(0.95, 0.05),
        legend.background = element_rect(colour = "black", size = 0.5),
        legend.key = element_rect(fill = NA))+
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.647"))) +
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = 0.278"))) +
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value < 0.001")))

SR_plot


#### Data of 2015 and 2016

# Calculate species richness
env_2015$richness = rowSums(ifelse(comm_2015_16 > 0, 1, 0))

# Log-transform the variables (log10) and
# Calculate simple linear regression model (sr ~ area)
mod2 = lm(log10(env_2015$richness) ~ log10(env_2015$area.ha))
summary(mod2)

# Test the normality of the residuals from the linear regression model 
# to ensure the validity of the statistical inference
hist(mod2$residuals)
shapiro.test(mod2$residuals) # It is normal (p=0.767)

# Plot the graph showing the relationship between number of bird 
# species and island area
SR_plot2 = 
  ggplot() +
  
  labs( x = "Island area (ha)",
        y = "Number of bird species (n)") +
  
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_log10(limits = c(1,50),
                breaks = c(1, 3, 10, 30 , 50),
                labels = c("1", "3", "10", "30", "50")) +
  annotation_logticks() +
  
  geom_smooth(data = env_2015,
              mapping = aes(x = area.ha, y = richness),
              method = lm,
              size = 1.5,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = env_2015, aes(x = area.ha,
                                  y = richness)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = env_2015, aes(x = area.ha,
                                  y = richness)) +
  
  theme_bw(base_size = 20) +
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm")) +
  
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = c(0.95, 0.05),
        legend.justification = c(0.95, 0.05),
        legend.background = element_rect(colour = "black", size = 0.5),
        legend.key = element_rect(fill = NA))+
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.697"))) +
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = 0.246"))) +
  
  annotate("text", x = 200, y = 3,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value < 0.001")))

SR_plot2


ggpubr::ggarrange(SR_plot2, SR_plot, ncol = 2, nrow = 1,
                  labels = c("2015-16", "2023-24"))




