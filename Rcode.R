# Diversity and structuring process
# of understory bird assemblages in Balbina
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: January 28, 2025

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
trees <- ape::read.nexus(unzip("trees_vertlife.zip", "output.nex"))


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

# Calculating Faith's PD for each site using 1000 trees.
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
summary(mod2)  # adjr=0.5797; p < 0.001; z = 0.19702


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
# The variation of PD values around the mean for each island is minimal. Therefore, I disregarded this part and will only use the mean.
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
  
  annotate("text", x = 500, y = 200,
           hjust = 0, vjust = 0, size = 6,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.579"))) +
  
  annotate("text", x = 500, y = 200,
           hjust = 0, vjust = 2, size = 6,
           parse = T, label = as.character(expression(italic(z)*"-value = 0.197"))) +
  
  annotate("text", x = 500, y = 200,
           hjust = 0, vjust = 4, size = 6,
           parse = T, label = as.character(expression(italic(p)*"-value < 0.001")))

PD_plot



########## CONTINUE FROM HERE (JANUARY 28, 2025)



plot.pd.sr <- ggpubr::ggscatter(pd_sr, x = "pd", y = "sr",
                                xlab = "Phylogenetic Diversity",
                                ylab = "Number of species",
                                add = "reg.line",
                                conf.int = TRUE,
                                color = "#636363",
                                shape = 16) +
  stat_cor(p.accuracy = 0.001) +
  theme_pubr()

ggarrange(plot.sr.log, plot.pd.log, ncol = 2, nrow = 1)

# Null model - Matrix-swap model
swapmodel <- list()
set.seed(13)
for (i in 1:999) {
  swapmodel[[i]] <- picante::randomizeMatrix(t(comm), null.model = "independentswap")
}

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
