##############################
# Code for water mites data
##############################
#
# Analyzing the prevalence and intensity of water mites on E. civile across landcover 
# context and host sex. 
# For landcover context, separate files were generated per type. 
#
# We used a chi-square and Kruskal-Wallis test. 
#

getwd()
water.mites.df <- read.csv("C:/yourdirectory")


# total water mites intensity, standard error, and prevalence.
mean(water.mites.df$Mites.)
se <- function(x) sqrt(var(x)/length(x))
se(water.mites.df$Mites.) 
range(water.mites.df$Mites.)
prevalence.total <- length(which(water.mites.df$Mites == 1))/nrow(water.mites.df) * 100
prevalence.total

## Kruskal-Wallis test with a post-hoc Dunn's test.

mites.kr <- kruskal.test(Mites. ~ as.factor(Landuse), data = water.mites.df)
summary(mites.kr)

hostsex.kr <- kruskal.test(Mites. ~ as.factor(Sex), data = water.mites.df)

# chi-square test for binary response and binary predictor variable.
chisq.test(water.mites.df$Landuse, water.mites.df$Mites, correct = FALSE)
chisq.test(water.mites.df$Sex, water.mites.df$Mites, correct = FALSE)
