#############################
# Code for water mites data
############################# 
# 
# Objective:
# 
# Analyze the prevalence and intensity of water mites on Enallagma civile across landcover context and host sex. 
# 
# For landcover context, host sex, and years sampled, subset files were generated. 
#
# Intensity:
# 
# We used bias-corrected and accelerated bootstrap for confidence intervals regarding water mites.
# For statistical significance testing, we used a Wilcox test and a bootstrapped t-test for water mite intensity.
# 
# Prevalence:
# 
# We used a Clopper-Pearson confidence interval, which takes into account the binary nature of prevalence data
# A Fischer's Exact test was used to determine any statistically significant difference between water mite prevalence and predictors.
# 

library(GenBinomApps)
library(coxed)
library(MKinfer)
library(ape)

getwd()
water.mites.df <- read.csv("C:/yourdirectory/e_civile_mites.csv")

# Subset example:

# grass <- subset(water.mites.df, Landuse == "Grassland")
# View(grass)

# Water mite intensity

mean(water.mites.df$Mites.) # Mean intensity of water mites on E. civile 
sd(water.mites.df$Mites.) # Standard deviation of water mites on E. civile

# Use mean and sd to determine the bias-corrected and accelerated bootstrap
theta <- rnorm(1000, mean, sd)
bca(theta, conf.level = 0.95) 

# Wilcox test to quantify the intensity of water mites by a predictor (e.g., landuse, host sex)

wilcox.test(Mites. ~ as.factor(Landuse), data = water.mites.df, exact = TRUE)

wilcox.test(Mites. ~ as.factor(Sex), data = water.mites.df, exact = TRUE)

# Bootstrapped t-test for water mite intensity by a predictor (e.g., landuse, host sex)

boot.t.test(Mites. ~ as.factor(Landuse), data = intensit)

boot.t.test(Mites. ~ as.factor(Sex), data = mites)

# Water mite prevalence

length(which(water.mites.df$Mites == 1))/nrow(water.mites.df) * 100 # Prevalence of water mites on E. civile


clopper.pearson.ci(k = 39, # k is the # of failures and/or successes (e.g., parasitized or not)
                   n = 100000, alpha = 0.05,
                   CI = "two.sided")

# Fisher's exact test to quantify the prevalence of water mites by a predictor (e.g., landuse, host sex)

fisher.test(water.mites.df$Landuse, water.mites.df$Mites)

fisher.test(water.mites.df$Sex, water.mites.df$Mites)

# Spatial autocorrelation

# Generate distance matrix, then calculate the inverse of the matrix values,
# replacing infinite values (from having multiple samples at some sites) and diagonals with zeros

mites.dists <- as.matrix(dist(cbind(water.mites.df$Lon, water.mites.df$Lat)))
water.mites.inv <- 1/mites.dists
water.mites.inv[!is.finite(water.mites.inv)] <- 0
diag(water.mites.inv) <- 0

Moran.I(water.mites.df$Mites., water.mites.inv)
