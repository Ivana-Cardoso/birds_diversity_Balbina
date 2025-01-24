# Organizing bird captures data
# Balbina 2023-2024
#
# Ivana Cardoso
# ivanawaters@gmail.com
#
# Last modified: January 24, 2025

# Clean environment
rm(list = ls()) # Clear all objects in memory
gc() # Garbage collection to free memory

# Load necessary packages
library(dplyr)
library(patchwork)
library(data.table)
library(terra)
library(landscapemetrics)
library(sf)

#### STEP 1 - ORGANIZE COMMUNITY AND TRAIT DATA ####
# Set working directory and importing data
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina")
data = read.csv("balbina_bird_captures_2023_2024.csv")
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Div.funcional/data")
avonet = read.csv("AVONET_BirdLife.csv")

# Removing blank spaces
data$Especie = gsub(" ", "_", data$Especie)
avonet$Species1 = gsub(" ", "_", avonet$Species1)
colnames(avonet)[2] = "Especie"

# Removing recaptured birds
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

# Renaming birds following BirdLife 2016 nomenclature
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

# Calculating mean values of collected traits
collected_traits = as.data.table(data_subset[,c(4, 10:21)]) # Only numerical data
collected_traits$Peso_real = collected_traits$Peso - collected_traits$Saco # Calculating the bird's actual weight (bird's weight in the bag - weight of the bag)
collected_traits = collected_traits[,c(1,14, 4:13)]

mean_traits = collected_traits[, lapply(.SD, mean, na.rm = T), by = Especie]

traits = merge(traits, mean_traits, by = "Especie")
traits = traits[,-c(6:11, 23:27,29:30,33:38)]

setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina")
write.csv(comm, file = "balbina_comm.csv", row.names = TRUE)
write.csv(traits, file = "balbina_understory_bird_traits.csv", row.names = TRUE)


#### STEP 2 - CALCULATE ENVIRONMENTAL VARIABLE: AREA ####
# Load necessary packages
library(terra)
library(landscapemetrics)
library(sf)

# Import raster data from MAPBIOMAS 2022 Sentinel (10m)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/MAPBIOMAS_raster")
Balbina <- rast("balbina_sentinel_2022.tif")

# Check the CRS of the imported raster
crs(Balbina)

# Reproject the raster to UTM Zone 21S with nearest neighbor resampling
Balbina <- project(Balbina, "EPSG:32721", method = "near")

# Round the raster values to avoid fractional classes
Balbina <- round(Balbina)

# Check unique values in the reprojected raster
unique_values <- unique(values(Balbina))
print(unique_values)

# Check if the raster has the correct properties
check_landscape(Balbina)

# Select only forest pixels (class 3)
forest_formation <- Balbina == 3

# Check properties
check_landscape(forest_formation)

# If no longer needed, remove the 'Balbina' object from memory to free up space (optional)
remove(Balbina)

# Import environmental variables for the sites
setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/Data")
sites <- read.csv("balbina_environmental_variables_2023_2024.csv")

# Prepare coordinates data frame from the imported sites
coordinates <- data.frame(
  X = sites$longitude.WGS84, 
  Y = sites$latitude.WGS84, 
  id = sites$site
)

# Remove rows corresponding to Continuous Forest sites (7 to 11)
coordinates <- coordinates[-c(7:11),] 

# Transform points into spatial coordinates
coordinates_wgs84 <- st_as_sf(coordinates, coords = c("X", "Y"), crs = 4326)
coordinates_utm <- st_transform(coordinates_wgs84, crs = 32721)

# Calculate the area of the island where each coordinate is located
area <- extract_lsm(forest_formation, y = coordinates_utm, 
                    what = "lsm_p_area", extract_id = coordinates_utm$id)

# Rename columns to merge area and sites data frames
colnames(sites)[1] = "extract_id"
sites <- merge(sites, area, by = "extract_id")
sites <- sites[,c(1:3,9)]
colnames(sites)[1] <- "sites"
colnames(sites)[4] <- "area_2022"

setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/Data")
write.csv(sites, "balbina_environmental_variables_2023_2024.csv")

