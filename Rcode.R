# Diversity and structuring process
# of understory bird assemblages in Balbina
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 03, 2025

# Clean environment
rm(list = ls()) # Clear all objects in memory
gc() # Garbage collection to free memory

# Disable scientific notation
options(scipen = 999)

# Load required packages
library(ape)
library(ggplot2)
library(picante)
library(ggpubr)
#library(reshape2)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/Data")

# Import data
comm <- read.csv("balbina_comm_islands_2023-24.csv", row.names = 1)
env <- read.csv("balbina_environmental_islands_2022.csv", row.names = 1)
trees <- ape::read.nexus(unzip("trees_vertlife_hackett.zip", "output.nex"))


#### TAXONOMIC DIVERSITY - Species Richness ####

# Calculate species richness
env$richness = rowSums(ifelse(comm > 0, 1, 0))

# Log-transform the variables (log10) and
# Calculate simple linear regression model (sr ~ area)
mod1 = lm(log10(env$richness) ~ log10(env$area.ha.2022))
summary(mod1) #adjr=0.6466; p < 0.001; z = 0.27812

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
  scale_y_log10() +
  annotation_logticks() +
  
  geom_smooth(data = env,
              mapping = aes(x = area.ha.2022, y = richness),
              method = lm,
              size = 1.5,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = env, aes(x = area.ha.2022,
                             y = richness)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = env, aes(x = area.ha.2022,
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



#### PHYLOGENETIC DIVERSITY - Faith's PD ####

# PD reflects the total phylogenetic richness of the community, calculated as the sum of tree branches length linking the species of the community.
#
# A high PD value indicates that the community includes a large amount of unique 
# evolutionary history.
#
# A low PD value indicates that the community includes less unique 
# evolutionary history.


# Fix the species names that are different between the community data (comm) 
# and the trees. Trees follow VertLife nomenclature, but community data 
# follows the Comite Brasileiro de Registros Ornitologicos (CBRO), 2021.
setdiff(colnames(comm), trees[[13]]$tip.label)

names(comm)[names(comm) == "Ceratopipra_erythrocephala"] <- "Pipra_erythrocephala"
names(comm)[names(comm) == "Cercomacroides_tyrannina"] <- "Cercomacra_tyrannina"
names(comm)[names(comm) == "Certhiasomus_stictolaemus"] <- "Deconychura_stictolaema"
names(comm)[names(comm) == "Chrysuronia_versicolor"] <- "Amazilia_versicolor"
names(comm)[names(comm) == "Cyanoloxia_rothschildii"] <- "Cyanocompsa_cyanoides"
names(comm)[names(comm) == "Epinecrophylla_gutturalis"] <- "Myrmotherula_gutturalis"
names(comm)[names(comm) == "Isleria_guttata"] <- "Myrmotherula_guttata"
names(comm)[names(comm) == "Maschalethraupis_surinamus"] <- "Tachyphonus_surinamus"
names(comm)[names(comm) == "Myrmoderus_ferrugineus"] <- "Myrmeciza_ferruginea"
names(comm)[names(comm) == "Pheugopedius_coraya"] <- "Thryothorus_coraya"
names(comm)[names(comm) == "Pseudopipra_pipra"] <- "Pipra_pipra"
names(comm)[names(comm) == "Sporophila_angolensis"] <- "Oryzoborus_angolensis"
names(comm)[names(comm) == "Tamatia_tamatia"] <- "Bucco_tamatia"
names(comm)[names(comm) == "Troglodytes_musculus"] <- "Troglodytes_aedon"

setdiff(colnames(comm), trees[[13]]$tip.label)

# Calculating Faith's PD for each site using 1000 trees
# For each tree, there are 32 PD and SR values
pds <- lapply(trees, function(x) pd(comm, x))  

# Select only PD values
pd <- lapply(pds, function(x) x[, 1])

# Transform it into a table with sites in the rows and PD 
# values from each tree in the columns
pd_df <- do.call(cbind, pd)

# Calculate mean PD for each site
pd_mean <- rowMeans(pd_df)

# Transform into a data frame
pd_sr <- data.frame(pd = pd_mean, sr = pds[[1]][, 2])

# Calculating confidence interval (1000 PD values) for each island
t <- list()
conf.inter <- list()

for (i in 1:32) {
  t[[i]] <- t.test(pd_df[i, ])
  conf.inter <- lapply(t, function(x) x[["conf.int"]])
}

conf.inter <- do.call(cbind, conf.inter)
conf.inter <- as.data.frame(t(conf.inter))
pd_sr$ymin <- conf.inter$V1
pd_sr$ymax <- conf.inter$V2

# How many values deviate from the mean for each island? 
# (Coefficient of variation around the mean)
coef_var <- sapply(1:32, function(i) (sd(pd_df[i,])/mean(pd_df[i, ]))*100)
coef_var <- as.data.frame(coef_var)
row.names(coef_var) <- row.names(comm)

# Testing if the values are normal
normal <- apply(pd_df, 1, shapiro.test)
shapiro_p <- lapply(normal, function(x) x[[2]])
shapiro_p_values <- do.call(c, shapiro_p)

subset(shapiro_p_values, shapiro_p_values < 0.05)  # All p-values were < 0.05, so not normal

# Calculating the Mode of the 1000 PD values
pd_mode <- lapply(1:32, function(i) {
  hist_info <- hist(pd_df[i, ], breaks = seq(min(pd_df[i, ]), max(pd_df[i, ]), length.out = 11))
  hist_info$mids[which.max(hist_info$counts)]
})
pd_mode <- do.call(cbind, pd_mode)
pd_mode = as.vector(pd_mode)

# Testing if mean and mode PD are correlated
plot(pd_mode ~ pd_mean)
cor.test(pd_mean, pd_mode)  # PD mean and mode are extremely correlated (0.999), so I will use PD mean.

# Plotting PD x Area, PD X SR
pd_sr$area <- env$area.ha.2022

row.names(pd_sr) <- row.names(env)

mod2 <- lm(log10(pd_sr$pd) ~ log10(pd_sr$area))
summary(mod2)  # adjr = 0.5802; p < 0.001; z = 0.18473


PD_plot = 
  ggplot() +
  
  labs(x = "Island area (ha)",
       y = "Mean Phylogenetic Diversity (Faith's PD)") +
  
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_log10() +
  annotation_logticks() +
  
  geom_smooth(data = pd_sr,
              mapping = aes(x = area, y = pd),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = pd_sr, aes(x = area,
                               y = pd)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = pd_sr, aes(x = area,
                               y = pd)) +
  
# geom_errorbar(data = pd_sr, aes(x = area, ymin = ymin, ymax = ymax),
#              col = "red") + 
# The variation of PD values around the mean for each island is minimal. 
# Therefore, I disregarded this part and will only use the mean.
  
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
  
  annotate("text", x = 200, y = 200,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.580"))) +
  
  annotate("text", x = 200, y = 200,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = 0.185"))) +
  
  annotate("text", x = 200, y = 200,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value < 0.001")))

PD_plot

# Visualize the Pearson correlation (R) between PD and SR
PDxSR <- ggpubr::ggscatter(pd_sr, x = "pd", y = "sr",
                                xlab = "Phylogenetic Diversity",
                                ylab = "Number of species",
                                add = "reg.line",
                                conf.int = TRUE,
                                color = "#636363",
                                shape = 16) +
  stat_cor(p.accuracy = 0.001) +
  theme_pubr()
PDxSR # Highly correlated (R = 0.97)


#### MEAN PAIRWISE DISTANCE (MPD) - Webb et al. 2002 ####

# MPD reflects the mean phylogenetic diversity among the species in the community
#
# A high MPD indicates that the species are phylogenetically distant from each 
# other, suggesting a more diverse community in evolutionary terms.
#
# A low MPD indicates that the species are phylogenetically close, suggesting 
# a less diverse community in evolutionary terms.

# Calculating MPD for each site using 1000 trees
mpds <- lapply(trees, function(x) mpd(comm, cophenetic(x)))

# Transform it into a table with sites in the rows and MPD 
# values from each tree in the columns
mpd_df <- do.call(cbind, mpds)
mpd_df = mpd_df[-9,] # One site, Formiga, has 0 MPD, because only one species were captured there. So, there is no way to calculate MPD. I will have to exclude

# Calculate mean MPD for each island
mpd_mean <- as.data.frame(rowMeans(mpd_df))

# Calculating confidence interval (1000 MPD values) for each island
t2 <- list()
conf.inter2 <- list()

for (i in 1:31) {
  t2[[i]] <- t.test(mpd_df[i, ])
  conf.inter2 <- lapply(t2, function(x) x[["conf.int"]])
}

conf.inter2 <- do.call(cbind, conf.inter2)
conf.inter2 <- as.data.frame(t(conf.inter2))
mpd_mean$ymin <- conf.inter2$V1
mpd_mean$ymax <- conf.inter2$V2

# How many values deviate from the mean for each island? 
# (Coefficient of variation around the mean)
coef_var2 <- sapply(1:31, function(i) (sd(mpd_df[i,])/mean(mpd_df[i, ]))*100)
coef_var2 <- as.data.frame(coef_var2)
row.names(coef_var2) <- row.names(comm)

# Testing if the values are normal
normal <- apply(mpd_df, 1, shapiro.test) 
shapiro_p <- lapply(normal, function(x) x[[2]])
shapiro_p_values <- do.call(c, shapiro_p)

subset(shapiro_p_values, shapiro_p_values < 0.05)  # All p-values were < 0.05, so not normal

# Calculating the Mode of the 1000 MPD values
mpd_mode <- lapply(1:31, function(i) {
  hist_info <- hist(mpd_df[i, ], breaks = seq(min(mpd_df[i, ]), max(mpd_df[i, ]), length.out = 11))
  hist_info$mids[which.max(hist_info$counts)]
})
mpd_mode <- do.call(cbind, mpd_mode)
mpd_mode = as.vector(mpd_mode)

# Testing if mean and mode MPD are correlated
plot(mpd_mode ~ mpd_mean[,1])
cor.test(mpd_mode, mpd_mean[,1])  # MPD mean and mode are extremely correlated (0.996), so I will use MPD mean.


mpd_mean$area <- env$area.ha.2022[-9]
colnames(mpd_mean)[1] = "mpd"
row.names(mpd_mean) <- row.names(env[-9,])

mod3 <- lm(log10(mpd_mean$mpd) ~ log10(mpd_mean$area))
summary(mod3)  # adjr=0.096; p = 0.0497; z = -0.014852


MPD_plot = 
  ggplot() +
  
  labs(x = "Island area (ha)",
       y = "Mean Pairwise Distance (MPD mean)") +
  
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_log10() +
  annotation_logticks() +
  
  geom_smooth(data = mpd_mean,
              mapping = aes(x = area, y = mpd),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = mpd_mean, aes(x = area,
                               y = mpd)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = mpd_mean, aes(x = area,
                               y = mpd)) +
  
# geom_errorbar(data = mpd_mean, aes(x = area, ymin = ymin, ymax = ymax),
#               col = "red") + 
# The variation of MPD values around the mean for each island is minimal. 
# Therefore, I disregarded this part and will only use the mean.
  
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
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.096"))) +
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = -0.015"))) +
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value = 0.05")))

MPD_plot


#### MEAN NEAREST TAXON DISTANCE (MNTD) - Webb et al. 2002 ####

# MNTD reflects the average distance between an individual and the most closely related (non-conspecific) individual
#
# A high MNTD indicates that species in the community are distant from their closest relatives (phylogenetically dispersed community).
#
# A low MNTD indicates that species are close to their closest relatives (phylogenetically clustered community).

# Calculating MNTD for each site using 1000 trees
mntds <- lapply(trees, function(x) mntd(comm, cophenetic(x)))

# Transform it into a table with sites in the rows and MNTD 
# values from each tree in the columns
mntd_df <- do.call(cbind, mntds)
mntd_df = mntd_df[-9,] # One site, Formiga, has 0 MNTD, because only one species were captured there. So, there is no way to calculate it. I will have to exclude

# Calculate mean MNTD for each island
mntd_mean <- as.data.frame(rowMeans(mntd_df))

# Calculating confidence interval (1000 MNTD values) for each island
t3 <- list()
conf.inter3 <- list()

for (i in 1:31) {
  t3[[i]] <- t.test(mntd_df[i, ])
  conf.inter3 <- lapply(t3, function(x) x[["conf.int"]])
}

conf.inter3 <- do.call(cbind, conf.inter3)
conf.inter3 <- as.data.frame(t(conf.inter3))
mntd_mean$ymin <- conf.inter3$V1
mntd_mean$ymax <- conf.inter3$V2

# How many values deviate from the mean for each island? 
# (Coefficient of variation around the mean)
coef_var3 <- sapply(1:31, function(i) (sd(mntd_df[i,])/mean(mntd_df[i, ]))*100)
coef_var3 <- as.data.frame(coef_var3)
row.names(coef_var3) <- row.names(comm[-9,])

# Testing if the values are normal
normal <- apply(mntd_df, 1, shapiro.test) 
shapiro_p <- lapply(normal, function(x) x[[2]])
shapiro_p_values <- do.call(c, shapiro_p)

subset(shapiro_p_values, shapiro_p_values < 0.05)  # All p-values, but one, were < 0.05, so not normal

# Calculating the Mode of the 1000 MNTD values
mntd_mode <- lapply(1:31, function(i) {
  hist_info <- hist(mntd_df[i, ], breaks = seq(min(mntd_df[i, ]), max(mntd_df[i, ]), length.out = 11))
  hist_info$mids[which.max(hist_info$counts)]
})
mntd_mode <- do.call(cbind, mntd_mode)
mntd_mode = as.vector(mntd_mode)

# Testing if mean and mode MNTD are correlated
plot(mntd_mode ~ mntd_mean[,1])
cor.test(mntd_mode, mntd_mean[,1])  # MNTD mean and mode are extremely correlated (0.999), so I will use MNTD mean.


mntd_mean$area <- env$area.ha.2022[-9]
colnames(mntd_mean)[1] = "mntd"
row.names(mntd_mean) <- row.names(env[-9,])

mod4 <- lm(log10(mntd_mean$mntd) ~ log10(mntd_mean$area))
summary(mod4)  # adjr = 0.4398; p < 0.001; z = -0.11171


MNTD_plot = 
  ggplot() +
  
  labs(x = "Island area (ha)",
       y = "Mean Nearest Taxon Distance (MNTD mean)") +
  
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_log10() +
  annotation_logticks() +
  
  geom_smooth(data = mntd_mean,
              mapping = aes(x = area, y = mntd),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  
  geom_point(shape = 16, colour = "black", 
             size = 5, alpha = 0.3,
             data = mntd_mean, aes(x = area,
                                   y = mntd)) +
  geom_point(shape = 21, colour = "black", 
             size = 5, 
             data = mntd_mean, aes(x = area,
                                   y = mntd)) +
  
# geom_errorbar(data = mntd_mean, aes(x = area, ymin = ymin, ymax = ymax),
#                col = "red") + 
# The variation of MNTD values around the mean for each island is minimal. 
# Therefore, I disregarded this part and will only use the mean.
  
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
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.439"))) +
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = -0.112"))) +
  
  annotate("text", x = 200, y = 160,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value < 0.001")))

MNTD_plot


ggpubr::ggarrange(SR_plot, NULL, NULL,
                  PD_plot, MPD_plot, MNTD_plot,
                  ncol = 3, nrow = 2,
                  heights = c(1,1))



#### PHYLOGENETIC DIVERSITY - Regional null model ####
# 
# I will use the Regional null model as proposed by Elliot Miller et al at 
# http://onlinelibrary.wiley.com/doi/10.1111/ecog.02070/abstract
#
# The regional null model simulates neutral dispersal of individuals into plots 
# by sampling from a regional pool that is unaffected by local dynamics
#
# To install the metricTester package, which is no longer available on CRAN, 
# I had to download the spacodiR package, which is also not on CRAN. 
# Additionally, I had to use an older version of geiger that includes the 
# rescale function, which metricTester calls during installation.
#
# More information about the functions and arguments can be found in:
# http://cran.nexr.com/web/packages/metricTester/metricTester.pdf
# http://www.ecography.org/sites/ecography.org/files/appendix/ecog-02070.pdf

library(devtools)
install_github("gdauby/spacodiR")
library(spacodiR)

remotes::install_version("geiger", version = "2.0.10")  # Replace with the version that has rescale function

install_github("eliotmiller/metricTester", dependencies = TRUE)
library(metricTester)

# Calculate regional null model
abundance_vector = metricTester::abundanceVector(comm)

null_model = list()
set.seed(13)

for (i in 1:999) {
  null_model[[i]] <- metricTester::oldRegionalNull(picante.cdm = comm,
                                                   tree = trees[[i]],
                                                   regional.abundance = abundance_vector)}




expPDs <- lapply(trees, function(x) ses.pd(comm, x, null.model = "independentswap",
                                           runs = 999, iterations = 1000,
                                           include.root = FALSE))

expPD <- lapply(expPDs, function(x) x[, 3])  # Select only expected PD values
expPD <- do.call(cbind, expPD)
expPD_mean <- rowMeans(expPD)
pd_sr$expPD <- expPD_mean

# SESPD = Obs-Exp/SDexp, where Obs = the observed PD, Exp = the mean of the 999 simulated values, and SDexp = the standard deviation of this mean.
SESPD <- (pd_sr$pd - pd_sr$expPD) / sd(pd_sr$expPD)
pd_sr$SESPD <- SESPD

mod.sespd <- lm(pd_sr$SESPD ~ log(pd_sr$area))
summary(mod.sespd)  # p=0.133, adjr=0.06158

# Plotting sesPD x Area, PD X sesPD
t.test(pd_sr$SESPD)


pd_sr$color[which(pd_sr$SESPD > 0.2)] = "blue" #overdispersion
pd_sr$color[which((pd_sr$SESPD > -0.2) & (pd_sr$SESPD < 0.2))] = "black" #null model
pd_sr$color[which(pd_sr$SESPD < -0.2)] = "red" #clustering
pd_sr$color = as.factor(pd_sr$color)



plot.sespd.log <- ggplot(data = pd_sr,
                         mapping = aes(x = area, y = SESPD)) +
  labs(x = "Island area (ha)",
       y = "Mean Standardized Effect Size of PD (sesPD)") +
  scale_x_log10(limits = c(NA, NA),
                breaks = c(1, 3, 10, 30, 100, 300, 1000, 16988.40),
                labels = c("1", "3", "10", "30", "100", "300", "1000", "CF")) +
  annotation_logticks() +
  geom_point(data = pd_sr,
             mapping = aes(x = area, y = SESPD),
             color = "#000000", fill = pd_sr$color,
             shape = 21, size = 3) +
  geom_smooth(data = pd_sr[c(1:7, 13:38), ],
              mapping = aes(x = area, y = SESPD),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  annotate(geom = "text", x = 6, y = 0.5,
           label = "Adj R² = 0.06158, p = 0.133") +
  theme_pubr(base_size = 20) +
  theme(axis.ticks = element_line(size = 0.25))

plot.pd.sespd <- ggpubr::ggscatter(pd_sr, x = "pd", y = "SESPD",
                                   xlab = "Phylogenetic Diversity",
                                   ylab = "Mean Standardized Effect Size of PD (SESPD)",
                                   add = "reg.line",
                                   conf.int = TRUE,
                                   color = "#000000",
                                   fill = "#636363",
                                   shape = 16) +
  stat_cor(p.accuracy = 0.001) +
  theme_pubr()

plot.sr.sespd <- ggpubr::ggscatter(pd_sr, x = "sr", y = "SESPD",
                                   xlab = "Species richness",
                                   ylab = "Mean Standardized Effect Size of PD (SESPD)",
                                   add = "reg.line",
                                   conf.int = TRUE,
                                   color = "#000000",
                                   fill = "#636363",
                                   shape = 16) +
  stat_cor(p.accuracy = 0.001) +
  theme_pubr()

plot(sort(pd_sr$SESPD))

write.csv(pd_sr, "pd_sr.csv")
