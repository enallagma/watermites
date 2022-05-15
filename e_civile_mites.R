###############################
# Code for water mites analyses
############################### 
# 
# Objective:
# 
# Analyze the prevalence and intensity of water mite parasitism on Enallagma civile across land-use context and host sex. 
# 
# Intensity:
# 
# Non-infected individuals were not included in calculation of mean intensity.
# 
# We used bias-corrected and accelerated bootstrap for confidence intervals regarding water mites (package coxed).
# 
# For statistical significance testing of water mite intensity, we used a bootstrapped t-test (package MKinfer). 
# 
# Prevalence:
# 
# We used a Clopper-Pearson confidence interval (package GenBinomApps), which takes into account 
# the binary nature of prevalence data
# 
# A Fisher's exact test (package stats in base R) was used to determine any statistically significant difference  
# in water mite prevalence with landuse, host sex.
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

# To ignore 0's in MitesNum column (so as to calculate true intensity of infection rather than mite abundance):

is.na(water.mites.df[[10]]) <- water.mites.df[[10]] < 1

# To determine if pooling is justified:

# Fisher's exact probability test for equal proportions between years (prevalences):
res <- prop.test(x=c(16,33), n=c(48,82))
res
# Conclusion: no significant difference between years

# Fligner-Killeen test for homogeneity of variances between years:
fligner.test(water.mites.df$MitesNum ~ water.mites.df$Year)
# Conclusion: no significant difference between years

# Therefore, pooling data for both years is acceptable.

# Subsets:

grass <- subset(water.mites.df, Landuse == "Grassland")
crop <- subset(water.mites.df, Landuse == "Cropland")
males <- subset(water.mites.df, Sex == "M")
females <- subset(water.mites.df, Sex == "F")

# Water mite intensity over all sampling dates (both years):

mean(water.mites.df$MitesNum, na.rm=TRUE) # Mean intensity of water mites on E. civile 
sd(water.mites.df$MitesNum, na.rm=TRUE) # Standard deviation of water mites on E. civile

# Bootstrapped t-test for water mite intensity by landuse, host sex:

boot.t.test(water.mites.df$MitesNum ~ as.factor(water.mites.df$Landuse))

boot.t.test(water.mites.df$MitesNum ~ as.factor(water.mites.df$Sex))

# Use mean and sd to determine the bias-corrected and accelerated bootstrap:

theta <- rnorm(1000, mean, sd)
bca(theta, conf.level = 0.95) 

# Water mite prevalence over all sampling dates (both years):

length(which(water.mites.df$MitesPres == 1))/nrow(water.mites.df) * 100 # Average prevalence of water mites on E. civile

# Fisher's exact test to quantify the prevalence of water mites by landuse, host sex:

fisher.test(water.mites.df$Landuse, water.mites.df$Mites)

fisher.test(water.mites.df$Sex, water.mites.df$Mites)

binomCI(x = 49, n = 130, method = "clopper-pearson")

# Spatial autocorrelation

# Generate distance matrix, then calculate the inverse of the matrix values,
# replacing infinite values (from having multiple samples at some sites) and diagonals with zeros:

mites.dists <- as.matrix(dist(cbind(water.mites.df$Lon, water.mites.df$Lat)))
water.mites.inv <- 1/mites.dists
water.mites.inv[!is.finite(water.mites.inv)] <- 0
diag(water.mites.inv) <- 0
Moran.I(water.mites.df$MitesNum, water.mites.inv, na.rm=TRUE) #intensity
Moran.I(water.mites.df$MitesPres, water.mites.inv, na.rm=TRUE) #prevalence
