###############################
# Code for water mites analyses
############################### 
# 
# Objective:
# 
# Analyze the prevalence and intensity of water mite parasitism on Enallagma civile across land-use context and host sex. 
# 
# All p-values were Bonferroni-adjusted.
#
# Intensity:
# 
# Non-infected individuals were not included in calculation of mean intensity.
# 
# For statistical significance testing of water mite intensity, we used a bootstrapped t-test (package MKinfer). 
# Because this is a bootstrapped test, the p-value will vary every time the test is run.
# 
# We used bias-corrected and accelerated bootstrap for confidence intervals of mean intensity (package coxed).
# 
# Prevalence:
# 
# A Fisher's exact test (package stats in base R) was used to determine any statistically significant difference  
# in water mite prevalence with landuse, host sex.
# 
# We calculated 95% Clopper-Pearson confidence intervals (package GenBinomApps), taking into account 
# the binary nature of prevalence data.
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

# Subsets:
Yr2006 <- subset(water.mites.df, Year == "2006")
Yr2007 <- subset(water.mites.df, Year == "2007")

# 2006:

# Water mite intensity:
  
# Bootstrapped t-test for water mite intensity by landuse, host sex:
boot.t.test(Yr2006$MitesNum ~ as.factor(Yr2006$Landuse), na.rm=TRUE)
boot.t.test(Yr2006$MitesNum ~ as.factor(Yr2006$Sex), na.rm=TRUE)

mean06 <- mean(Yr2006$MitesNum, na.rm=TRUE) # Mean intensity of water mites on E. civile 
sd06 <- sd(Yr2006$MitesNum, na.rm=TRUE) # Standard deviation of water mites on E. civile

# Use mean and sd to determine the bias-corrected and accelerated bootstrap 95% CI:
theta06 <- rnorm(1000, mean06, sd06)
bca(theta06, conf.level = 0.95) 

# Water mite prevalence:

# Average prevalence of water mites on E. civile:
mean(Yr2006$MitesPres, na.rm=TRUE)

# Fisher's exact test to quantify the prevalence of water mites by landuse, host sex:
fisher.test(Yr2006$Landuse, Yr2006$MitesPres)
fisher.test(Yr2006$Sex, Yr2006$MitesPres)
clopper.pearson.ci(16, 48, alpha=0.05)

# 2007:

# Water mite intensity:

# Bootstrapped t-test for water mite intensity by landuse, host sex:
boot.t.test(Yr2007$MitesNum ~ as.factor(Yr2007$Landuse), na.rm=TRUE)
boot.t.test(Yr2007$MitesNum ~ as.factor(Yr2007$Sex), na.rm=TRUE)

mean07 <- mean(Yr2007$MitesNum, na.rm=TRUE) # Mean intensity of water mites on E. civile 
sd07 <- sd(Yr2007$MitesNum, na.rm=TRUE) # Standard deviation of water mites on E. civile

# Use mean and sd to determine the bias-corrected and accelerated bootstrap 95% CI:
theta07 <- rnorm(1000, mean07, sd07)
bca(theta07, conf.level = 0.95) 

# Water mite prevalence:

# Average prevalence of water mites on E. civile:
mean(Yr2006$MitesPres, na.rm=TRUE)

# Fisher's exact test to quantify the prevalence of water mites by landuse, host sex:
fisher.test(Yr2007$Landuse, Yr2007$MitesPres)
fisher.test(Yr2007$Sex, Yr2007$MitesPres)
clopper.pearson.ci(16, 82, alpha=0.05)

# Years pooled:

# Water mite intensity over all sampling dates (both years):

# Bootstrapped t-test for water mite intensity by landuse, host sex:
boot.t.test(water.mites.df$MitesNum ~ as.factor(water.mites.df$Landuse), na.rm=TRUE)
boot.t.test(water.mites.df$MitesNum ~ as.factor(water.mites.df$Sex), na.rm=TRUE)

mean <- mean(water.mites.df$MitesNum, na.rm=TRUE) # Mean intensity of water mites on E. civile 
sd <- sd(water.mites.df$MitesNum, na.rm=TRUE) # Standard deviation of water mites on E. civile

# Use mean and sd to determine the bias-corrected and accelerated bootstrap 95% CI:
theta <- rnorm(1000, mean, sd)
bca(theta, conf.level = 0.95) 

# Water mite prevalence over all sampling dates (both years):

# Average prevalence of water mites on E. civile:
mean(water.mites.df$MitesPres, na.rm=TRUE)

# Fisher's exact test to quantify the prevalence of water mites by landuse, host sex:
fisher.test(water.mites.df$Landuse, water.mites.df$Mites)
fisher.test(water.mites.df$Sex, water.mites.df$Mites)
clopper.pearson.ci(49, 130, alpha=0.05)

# Spatial autocorrelation

# Generate distance matrix, then calculate the inverse of the matrix values,
# replacing infinite values (from having multiple samples at some sites) and diagonals with zeros:
#
# 2006:
mites.dists.06 <- as.matrix(dist(cbind(Yr2006$Lon, Yr2006$Lat)))
water.mites.inv.06 <- 1/mites.dists.06
water.mites.inv.06[!is.finite(water.mites.inv.06)] <- 0
diag(water.mites.inv.06) <- 0
Moran.I(Yr2006$MitesNum, water.mites.inv.06, na.rm=TRUE) #intensity
Moran.I(Yr2006$MitesPres, water.mites.inv.06, na.rm=TRUE) #prevalence
#
# 2007:
mites.dists.07 <- as.matrix(dist(cbind(Yr2007$Lon, Yr2007$Lat)))
water.mites.inv.07 <- 1/mites.dists.07
water.mites.inv.07[!is.finite(water.mites.inv.07)] <- 0
diag(water.mites.inv.07) <- 0
Moran.I(Yr2007$MitesNum, water.mites.inv.07, na.rm=TRUE) #intensity
Moran.I(Yr2007$MitesPres, water.mites.inv.07, na.rm=TRUE) #prevalence
#
