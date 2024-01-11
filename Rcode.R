# Diversity and structuring process
# of understory bird assemblages in Balbina
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: January 11, 2024

# Set working directory
setwd("C:/Users/Ivana/OneDrive/PhD_INPA/1.Diversity_question/Analises/birds_diversity_Balbina")

# Load packages
library(ape)
library(reshape2)
library(picante)
library(ggpubr)

# Import data
trees <- ape::read.nexus(unzip("1000_trees_vertlife.zip", "output.nex"))
comm <- read.csv("https://ndownloader.figshare.com/files/15158531")
env <- read.csv("https://ndownloader.figshare.com/files/15158528", row.names = 1)

# Handling data
comm <- subset(comm, comm$status != "Recapture")  # Removing recaptured birds
comm <- reshape2::dcast(comm, site ~ species, value.var = "species")
row.names(comm) <- comm$site
row.names(comm) == row.names(env)

# Calculating Faith's PD
pds <- lapply(trees, function(x) pd(comm, x))  # Calculate the PD for each site using the 1000 trees. For each tree, there are 38 PD and SR values.
pd <- lapply(pds, function(x) x[, 1])  # Select only PD values
pd_df <- do.call(cbind, pd)  # Transform it into a table with sites in the rows and PD values from each tree in the columns.

pd_mean <- rowMeans(pd_df)  # Calculate mean PD for each site

pd_sr <- data.frame(pd = pd_mean, sr = pds[[1]][, 2])

# Calculating confidence interval (1000 PD values) for each island
t <- list()
conf.inter <- list()

for (i in 1:38) {
  t[[i]] <- t.test(pd_df[i, ])
  conf.inter[[i]] <- lapply(t, function(x) x[["conf.int"]])
}

conf.inter <- do.call(cbind, conf.inter)
conf.inter <- as.data.frame(t(conf.inter))
pd_sr$ymin <- conf.inter$V1
pd_sr$ymax <- conf.inter$V2

# How many values deviate from the mean for each island? (Coefficient of variation around the mean)
coef_var <- sapply(1:38, function(i) (sd(pd_df[i,])/mean(pd_df[i, ]))*100)
coef_var <- as.data.frame(coef_var)
row.names(coef_var) <- row.names(comm)

# Testing if the values are normal
normal <- apply(pd_df, 1, shapiro.test)
shapiro_p <- lapply(normal, function(x) x[[2]])
shapiro_p_values <- do.call(c, shapiro_p)

subset(shapiro_p_values, shapiro_p_values < 0.05)  # All p-values were < 0.05, so not normal

# Calculating the Mode of the 1000 PD values
pd_mode <- lapply(1:38, function(i) {
  hist_info <- hist(pd_df[i, ], breaks = seq(min(pd_df[i, ]), max(pd_df[i, ]), length.out = 11))
  hist_info$mids[which.max(hist_info$counts)]
})

# Testing if mean and mode PD are correlated
plot(pd_mode ~ pd_mean)
cor.test(pd_mean, pd_mode)  # PD mean and mode are extremely correlated (0.999), so I will use PD mean.

# Plotting PD x Area, SR X Area, PD X SR
area <- env$area.ha
pd_sr$area <- area
pd_sr[8:12, 5] <- 16988.4  # Consider CF area 10 times greater than the largest island (to construct the graphs)

row.names(pd_sr) <- row.names(env)

mod.pd <- lm(log(pd_sr$pd) ~ log(pd_sr$area))
summary(mod.pd)  # p=8.344e-11, adjr=0.6862

mod.sr <- lm(log(pd_sr$sr) ~ log(pd_sr$area))
summary(mod.sr)  # p=2.616e-13, adjr=0.7715

plot.pd.log <-
  ggplot(data = pd_sr,
         mapping = aes(x = area, y = pd)) +
  labs(x = "Island area (ha)",
       y = "Mean Phylogenetic Diversity (PD)") +
  scale_x_log10(limits = c(NA, NA),
                breaks = c(1, 3, 10, 30, 100, 300, 1000, 16988.40),
                labels = c("1", "3", "10", "30", "100", "300", "1000", "CF")) +
  scale_y_log10(breaks = c(100, 300, 500, 700, 1100),
                labels = c("100", "300", "500", "700", "1100")) +
  annotation_logticks() +
  geom_point(data = pd_sr,
             mapping = aes(x = area, y = pd),
             color = "#000000", fill = "#000000",
             shape = 21, size = 3) +
  geom_smooth(data = pd_sr[c(1:7, 13:38), ],
              mapping = aes(x = area, y = pd),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  # geom_errorbar(aes(x = area, ymin = ymin, ymax = ymax),
  #               col = "red") + 
  # The variation of PD values around the mean for each island is minimal. Therefore, I disregarded this part and will only use the mean.
  annotate(geom = "text", x = 4, y = 1200,
           label = "Adj R² = 0.69, p < 0.001") +
  theme_pubr(base_size = 20) +
  theme(axis.ticks = element_line(size = 0.25))

plot.sr.log <- ggplot(data = pd_sr,
                      mapping = aes(x = area, y = sr)) +
  labs(x = "Island area (ha)",
       y = "Number of species") +
  scale_x_log10(limits = c(NA, NA),
                breaks = c(1, 3, 10, 30, 100, 300, 1000, 16988.40),
                labels = c("1", "3", "10", "30", "100", "300", "1000", "CF")) +
  scale_y_log10(breaks = c(10, 30, 50),
                labels = c("10", "30", "50")) +
  annotation_logticks() +
  geom_point(data = pd_sr,
             mapping = aes(x = area, y = sr),
             color = "#000000", fill = "#000000",
             shape = 21, size = 3) +
  geom_smooth(data = pd_sr[c(1:7, 13:38), ],
              mapping = aes(x = area, y = sr),
              method = lm,
              color = "#000000",
              fill = "#636363") +
  annotate(geom = "text", x = 4, y = 48,
           label = "Adj R² = 0.77, p < 0.001") +
  theme_pubr(base_size = 20) +
  theme(axis.ticks = element_line(size = 0.25))

plot.pd.sr <- ggpubr::ggscatter(pd_sr, x = "pd", y = "sr",
                                xlab = "Phylogenetic Diversity",
                                ylab = "Number of species",
                                add = "reg.line",
                                conf.int = TRUE,
                                color = "#636363",
                                shape = 16) +
  stat_cor(p.accuracy = 0.001) +
  theme_pubr()

ggarrange(plot.pd.log, plot.sr.log, ncol = 1, nrow = 2)

