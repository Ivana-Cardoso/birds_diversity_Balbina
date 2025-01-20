# Organizing bird captures data
# Balbina 2023-2024
#
# Ivana Cardoso
# ivanawaters@gmail.com
#
# Last modified: January 20, 2025

rm(list = ls())

# Load packages
library(dplyr)
library(patchwork)
library(data.table)

# Set working directory and importing data
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina")
data = read.csv("balbina_bird_captures_2023_2024.csv")
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Div.funcional/data")
avonet = read.csv("AVONET_BirdLife.csv")

# Organizing data
data$Especie = gsub(" ", "_", data$Especie)
avonet$Species1 = gsub(" ", "_", avonet$Species1)
colnames(avonet)[2] = "Especie"

data$Ocorrencia = 1
data_subset <- subset(data, Status != "Recaptura")


# Creating a community data matrix
comm = as.data.frame(tapply(data_subset$Ocorrencia,
                            list(data_subset$Sitio, data_subset$Especie),
                            FUN = sum))

comm[is.na(comm)] = 0

# Creating a trait data matrix
comm_t = as.data.frame(t(comm))
comm_t$Especie = row.names(comm_t)
row.names(comm_t) = NULL
comm_2 = comm_t[,c(38, 1:37)]
comm_2$Especie[comm_2$Especie == "Chrysuronia_versicolor"] = "Amazilia_versicolor"
comm_2$Especie[comm_2$Especie == "Maschalethraupis_surinamus"] = "Maschalethraupis_surinama"
comm_2$Especie[comm_2$Especie == "Troglodytes_musculus"] = "Troglodytes_aedon"
comm_2$Especie[comm_2$Especie == "Tamatia_tamatia"] = "Nystactes_tamatia"
comm_2$Especie[comm_2$Especie == "Thraupis_episcopus"] = "Tangara_episcopus"

traits = inner_join(comm_2, avonet, by = "Especie")
traits = traits[,-c(2:39)]

comm_2$Especie == traits$Especie

colnames(traits)[1] = "BIRDLIFE_2016"
traits$CBRO_2021 = comm_t$Especie
traits$Especie = comm_t$Especie
traits = traits[,c(38,37, 1:36)]

remove(comm_2)
remove(comm_t)
remove(avonet)

# Calculating mean values of collected traits
collected_traits = as.data.table(data_subset[,c(4, 10:21)]) # Only numerical data
collected_traits$Peso_real = collected_traits$Peso - collected_traits$Saco # Calculating the bird's actual weight (bird's weight in the bag - weight of the bag)
collected_traits = collected_traits[,c(1,14, 4:13)]

mean_traits = collected_traits[, lapply(.SD, mean, na.rm = T), by = Especie]

traits = merge(traits, mean_traits, by = "Especie")
traits = traits[,-c(6:11, 23:27,29:30,33:38)]

remove(collected_traits)
remove(mean_traits)
remove(data)
remove(data_subset)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina")
write.csv(comm, file = "balbina_comm.csv", row.names = TRUE, col.names = TRUE)
write.csv(traits, file = "balbina_understory_bird_traits.csv", row.names = TRUE)
