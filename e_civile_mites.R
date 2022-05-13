###############################
# Code for water mites analyses
############################### 
# 
# Objective:
# 
# Analyze the prevalence and intensity of water mites on Enallagma civile across land-use context and host sex. 
# 
# For land use, host sex, and years sampled, subset files were generated for analysis. 
#
# Intensity:
# 
# Non-infected individuals were not included in calculation of mean intensity.
# 
# We used bias-corrected and accelerated bootstrap for confidence intervals regarding water mites (package coxed).
# 
# For statistical significance testing, we used a Wilcox test and a bootstrapped t-test (package MKinfer) 
# for water mite intensity.
# 
# Prevalence:
# 
# We used a Clopper-Pearson confidence interval (package GenBinomApps), which takes into account 
# the binary nature of prevalence data
# 
# A Fisher's exact test (package stats in base R) was used to determine any statistically significant difference  
# between water mite prevalence and predictors.
# 
# Spatial autocorrelation:
#
# We calculated Moran's I (package ape) to test the null hypothesis of no spatial autocorrelation in the number 
# of mites present based on site proximity.
#

library(ape)
library(coxed)
library(GenBinomApps)
library(MKinfer)

water.mites.df <- read.csv("C:/yourdirectory/e_civile_mites.csv")

# To ignore 0's in MitesNum column (so as to calculate true intensity rather than mite abundance):

is.na(water.mites.df[[10]]) <- water.mites.df[[10]] < 1

# To test for differences between years to determine if pooling is possible:





# Subsets:

grass <- subset(water.mites.df, Landuse == "Grassland")
crop <- subset(water.mites.df, Landuse == "Cropland")


# Water mite intensity over all sampling dates (both years):

mean(water.mites.df$MitesNum, na.rm=TRUE) # Mean intensity of water mites on E. civile 
sd(water.mites.df$MitesNum, na.rm=TRUE) # Standard deviation of water mites on E. civile

# Use mean and sd to determine the bias-corrected and accelerated bootstrap:

theta <- rnorm(1000, mean, sd)
bca(theta, conf.level = 0.95) 

# Wilcox test to quantify the intensity of water mites by a predictor (e.g., landuse, host sex):

wilcox.test(MitesNum ~ as.factor(Landuse), data = water.mites.df, exact = TRUE)

wilcox.test(MitesNum ~ as.factor(Sex), data = water.mites.df, exact = TRUE)

# Bootstrapped t-test for water mite intensity by a predictor (e.g., landuse, host sex):

boot.t.test(MitesNum ~ as.factor(Landuse), data = intensit)

boot.t.test(MitesNum ~ as.factor(Sex), data = mites)

# Water mite prevalence over all sampling dates (both years):

length(which(water.mites.df$Mites == 1))/nrow(water.mites.df) * 100 # Prevalence of water mites on E. civile

clopper.pearson.ci(k = 39, # k is the # of failures and/or successes (e.g., parasitized or not)
                   n = 100000, alpha = 0.05,
                   CI = "two.sided")

# Fisher's exact test to quantify the prevalence of water mites by a predictor (e.g., landuse, host sex):

fisher.test(water.mites.df$Landuse, water.mites.df$Mites)

fisher.test(water.mites.df$Sex, water.mites.df$Mites)

# Spatial autocorrelation

# Generate distance matrix, then calculate the inverse of the matrix values,
# replacing infinite values (from having multiple samples at some sites) and diagonals with zeros:

mites.dists <- as.matrix(dist(cbind(water.mites.df$Lon, water.mites.df$Lat)))
water.mites.inv <- 1/mites.dists
water.mites.inv[!is.finite(water.mites.inv)] <- 0
diag(water.mites.inv) <- 0
Moran.I(water.mites.df$MitesNum, water.mites.inv, na.rm=TRUE) #intensity
Moran.I(water.mites.df$MitesPres, water.mites.inv) #prevalence
