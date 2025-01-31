#Script to run Klaza glm models
## EARLY WINTER
#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

#Reviewed Nov 20, 2024

library(easystats)
library(see)
library(effectsize)
library(parameters)
library(plyr)
library(tidyverse)
library(terra)
library(sf)
library(dplyr)
library(sp)
library(corrplot)
library(report)
library(terra)
library(predicts)
library(car)
library(AICcmodavg)

set.seed(123)

#Import used and available points with covariates previously extracted to them
ew_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_ew_data3.csv')

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(ew_data$landcover2015)
ew_data$landcover2015 <- factor(ew_data$landcover2015)

# View the current levels
levels(ew_data$landcover2015) #already set during raster prep

str(ew_data$firedec)
ew_data$firedec <- factor(ew_data$firedec)
levels(ew_data$firedec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(ew_data), 0.8*nrow(ew_data))
ew_train = ew_data[n,]
ew_test = ew_data[-n,] 

##########################################################
#Univariate Model Selection
#Univariate environmental glms - unscaled versions - only necessary for covariates where effect of different scales needs to be checked.
#Allows me to define a unique set of covariates for each season.

#Landcover
ew_lc30_1 <- glm(presence~landcover2015, family=binomial, data=ew_train)
ew_lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=ew_train)
ew_lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=ew_train)

summary(ew_lc30_1)
plot(ew_lc30_1)
summary(ew_lc100_1)
plot(ew_lc100_1)
summary(ew_lc250_1)
plot(ew_lc250_1) #LC 30 is the best based on AIC model selection. Low, negative z scores.
#AIC - Land cover
ew_lc_glms <- list(ew_lc30_1, ew_lc100_1, ew_lc250_1)
ew_lc_glms_names <- c('ew_lc30_1', 'ew_lc100_1', 'ew_lc250_1')
aictab(cand.set = ew_lc_glms, modnames = ew_lc_glms_names)

#Lichen
ew_lichbi30 <- glm(presence~lichen_bi, family=binomial, data=ew_train)
ew_lichenbi100 <- glm(presence~lichen_bi100, family=binomial, data=ew_train)
ew_lichenbi250 <- glm(presence~lichen_bi250, family=binomial, data=ew_train)
summary(ew_lich30_1)
summary(ew_lichbi30)# Lichen binomial 100 m is the best based on AIC model selection. High z scores. 
summary(ew_lich100_1)
summary(ew_lich250_1) 
#AIC - lichen PFT
ew_lich_glms <- list(ew_lichbi30, ew_lichenbi100, ew_lichenbi250)
ew_lich_glms_names <- c('ew_lichbi30', 'ew_lichenbi100', 'ew_lichenbi250')
aictab(cand.set = ew_lich_glms, modnames = ew_lich_glms_names)

#Conifer
ew_conbi30 <- glm(presence~conifer_bi, family=binomial, data=ew_train)
ew_conbi100 <- glm(presence~conifer_bi100, family=binomial, data=ew_train)
ew_conbi250 <- glm(presence~conifer_bi250, family=binomial, data=ew_train)
summary(ew_con30_1)
summary(ew_conbi100)# Conifer bivariate 250 m is the best based on AIC model selection. Medium negative z scores. 
summary(ew_conbi250)
#AIC - conifer PFT
ew_con_glms <- list(ew_conbi30, ew_conbi100, ew_conbi250)
ew_con_glms_names <- c('ew_conbi30', 'ew_conbi100', 'ew_conbi250')
aictab(cand.set = ew_con_glms, modnames = ew_con_glms_names)

#Deciduous shrub
ew_decidbi30 <- glm(presence~decid_bi, family=binomial, data=ew_train)
ew_decidbi100 <- glm(presence~decid_bi100, family=binomial, data=ew_train)
ew_decidbi250 <- glm(presence~decid_bi250, family=binomial, data=ew_train)
summary(ew_decidbi30) #deciduous shrub bivariate 250m is the best based on AIC model selection, 30m is close second. 
summary(ew_decid100_1)
summary(ew_decid250_1) 
#AIC - deciduous shrub PFT
ew_decid_glms <- list(ew_decidbi30, ew_decidbi100, ew_decidbi250)
ew_decid_glms_names <- c('ew_decidbi30', 'ew_decidbi100', 'ew_decidbi250')
aictab(cand.set = ew_decid_glms, modnames = ew_decid_glms_names)

#graminoid
ew_grambi30 <- glm(presence~gram_bi, family=binomial, data=ew_train)
ew_grambi100 <- glm(presence~gram_bi100, family=binomial, data=ew_train)
ew_grambi250 <- glm(presence~gram_bi250, family=binomial, data=ew_train)
summary(ew_grambi100) #graminoid bivariate 30 m is the best based on AIC model selection. 
ew_gram_glms <- list(ew_grambi30, ew_grambi100, ew_grambi250)
ew_gram_glm_names <- c('ew_grambi30', 'ew_grambi100', 'ew_grambi250')
aictab(cand.set = ew_gram_glms, modnames = ew_gram_glm_names)

#Fire
ew_fire_1 <- glm(presence ~ fire2015, family = binomial, data=ew_train)
ew_firedec_1 <- glm(presence ~ firedec, family = binomial, data=ew_train)
ew_firesum_1 <- glm(presence ~ firesum, family = binomial, data=ew_train)

summary(ew_fire_1) #fire variables all have very similar delta AIC scores, will test with other environmental variables to find best fit. 
summary(ew_firedec_1) 
summary(ew_firesum_1) 
ew_fire_glms <- list(ew_fire_1, ew_firedec_1, ew_firesum_1)
ew_fire_glms_names <- c('ew_fire_1', 'ew_firedec_1', 'ew_firesum_1')
aictab(cand.set = ew_fire_glms, modnames = ew_fire_glms_names)


#Include: elev (q), slope (q), landcover30, lichenbi100, conbi250, grambi30, fire
#excluding deciduous shrub because not biologically important in winter

##############################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Early Winter:

#Only include variables identified from above process
# Add the transformed terms to the dataset to check collinearity
ew_data$slope2 <- ew_data$slope^2 #quadratic transformation
ew_data$kdem2 <- ew_data$kdem^2

ew_dat.cor = ew_data[,c("kdem", "kdem2", "slope", "slope2", "landcover2015", "lichen_bi100", "conifer_bi250", "gram_bi", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015d", "road2015_250", "road2015_500", "road2015_1000")]
names(ew_dat.cor) = c("kdem", "kdem2", "slope", "slope2", "landcover2015", "lichen_bi100", "conifer_bi250", "gram_bi", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015d", "road2015_250", "road2015_500", "road2015_1000")
cor_matrix <- cor(ew_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(ew_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/ew_correlation_matrix.csv", row.names = TRUE)

########################################################################
## EARLY WINTER GLM MODEL SELECTION ##

#Initial exploration
#all variables
ew_test1 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#no slope or landcover
ew_test2 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#no slope
ew_test3 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#no binomial veg rasters
ew_test4 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015 + scale(fire2015), family=binomial,data=ew_train) 
#no fire
ew_test5 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi), family=binomial,data=ew_train) 
#no landcover
ew_test6 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 

summary(ew_test6)
#elevation is primary driver of model. Lichen also important. Neither slope term significant and only the grasslands and non-vegetated landcover classes are significant. 
#general agreement among landcover classes and binomial vegetation predictors
#slope consistently not significant, in some models both slope terms are negative. Proceed without. 
#fire has a significant negative effect


#Starting with GLM models first
#Unique model for each season
# Phase 1: Find the best environmental base model
#Starting with glm (no mixed effect), standardized predictors
#exlude distance to water and deciduous shrub
#include: elev (q), slope (q), landcover30, lichenbi100, conbi250, grambi30, fire
 
#All top performing predictors from evaluations and fire by year
ew_glm_1e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#All top performing predictors from evaluations, no fire
ew_glm_2e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi), family=binomial,data=ew_train) 
#All top performing variables from evaluations, no graminoid
ew_glm_3e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(fire2015), family=binomial,data=ew_train) 
#All top performing variables from evaluations, no conifer
ew_glm_4e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#All top performing variables from evaluations, no lichen 
ew_glm_5e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train)                 
#All top performing variables from evaluations, no landcover
ew_glm_6e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#all top performing variables from evaluations, no elevation
ew_glm_7e <- glm(presence~landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015), family=binomial,data=ew_train) 
#null model
ew_glm_8e <- glm(presence~1,family=binomial,data=ew_train) 
#no slope, no graminoid
#ew_glm_12e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(fire2015), family=binomial,data=ew_train)


ew_e_glms <- list(ew_glm_1e, ew_glm_2e, ew_glm_3e, ew_glm_4e, ew_glm_5e, ew_glm_6e, ew_glm_7e, ew_glm_8e)
ew_e_glms_names <- c('ew_glm_1e', 'ew_glm_2e', 'ew_glm_3e', 'ew_glm_4e', 'ew_glm_5e', 'ew_glm_6e', 'ew_glm_7e', 'ew_glm_8e')
aictab(cand.set = ew_e_glms, modnames = ew_e_glms_names)

summary(ew_glm_1e) #best performing model - all variables

#graminoid is least important variable but still contributes to model fit. Fire is second least important variable, then landcover. 
#lichen, elevation, conifer are all important variables
report(ew_glm_1e)
#Some firedec and landcover categories are not significant. 

##################################################
## Model Calibration - environmental base model
#plot predicted vs observed probabilities of occurrence:
# Step 1: Fit the GLM model on the training data (already done in your case with lw_train)
# Assuming the model is already fit and called lw_glm_1e

# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(ew_glm_1e, newdata = ew_test, type = "response")

# Step 3: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = ew_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()
# Plot a histogram of predicted probabilities
ggplot(results_test, aes(x = predicted)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(x = "Predicted Probabilities", y = "Frequency", 
       title = "Distribution of Predicted Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Bin the predicted probabilities using deciles.

results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()
#Good calibration, the model preforms well across range of probabilities. 
#Very Slight under prediction throughout model.

#Evaluation with AUC / ROC
#unscaled model:
ew_glm_1e_us <- glm(presence~kdem + I(kdem^2) + landcover2015 + lichen_bi100 + conifer_bi250 + gram_bi + fire2015, family=binomial,data=ew_train)
#evaluate
ew_glm_1e_eval <- pa_evaluate(predict(ew_glm_1e_us, ew_test[ew_test$presence==1, ]), predict(ew_glm_1e_us, ew_test[ew_test$presence==0, ]))
print(ew_glm_1e_eval)
plot(ew_glm_1e_eval, "ROC")
#AUC=0.952 #Very high AUC score, excellent performance
##################
#Model performance with easystats
r2(ew_glm_1e_us) #Tjur's R2 = 0.648



#Predicted occurrence map
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast2 <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
#envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kfirerast2015_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif')%>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast2$kdem) #dem to numeric - should have been done at raster prep script
landcover2015 <- as.factor(envrast2$landcover2015)
slope <- envrast2$slope
lichen_bi100 <- rast100m$lichen_bi100 |>
  resample(mask, method='near')
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
gram_bi <- envrast2$gram_bi
fire2015 <- firerast$fire2015

rasters <- c(kdem, slope, landcover2015, lichen_bi100, conifer_bi250, gram_bi, fire2015) #combine rasters. Make sure they are separated by commas.
#run predict function on disturbance model
p1 <- predict(rasters, ew_glm_1e_us, type="response")
plot(p1) #plots!

#Add gps points to map
#filtered used points w/o covariates extracted to them:
ew_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_ew_daily.shp')

# Plot the raster and add points on top
plot(st_geometry(ew_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

######################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation. 
#Only going to buffered disturbance 500 m because used/available distributions switched at 1000 m buffer

#disturbance 30 m buffer
ew_glm_1d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(dist2015_30), family=binomial,data=ew_train) 
#disturbance 250 m buffer
ew_glm_2d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(dist2015_250), family=binomial,data=ew_train) 
#disturbance 500 m buffer
ew_glm_3d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(dist2015_500), family=binomial,data=ew_train) 
#mining 30 m buffer
ew_glm_4d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_30), family=binomial,data=ew_train) 
#mining 250 m buffer
ew_glm_5d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_250), family=binomial,data=ew_train) 
#mining 500 m buffer
ew_glm_6d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_500), family=binomial,data=ew_train) 
#mining 1000 m buffer
ew_glm_7d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_1000), family=binomial,data=ew_train) 
#mining 2000 m buffer
ew_glm_8d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_2000), family=binomial,data=ew_train) 
#mining 3000 m buffer
ew_glm_9d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_3000), family=binomial,data=ew_train) 
#road 30 m buffer
ew_glm_10d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(road2015_30), family=binomial,data=ew_train) 
#road 250 m buffer
ew_glm_11d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(road2015_250), family=binomial,data=ew_train) 
#road 500 m buffer
ew_glm_12d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(road2015_500), family=binomial,data=ew_train) 
#road 1000 m buffer
ew_glm_13d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(road2015_1000), family=binomial,data=ew_train) 
#disturbance density 250 m
ew_glm_14d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(distden2015_250), family=binomial,data=ew_train) 
#disturbance density 500 m
ew_glm_15d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(distden2015_500), family=binomial,data=ew_train) 
#disturbance density 1000 m
ew_glm_16d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(distden2015_1000), family=binomial,data=ew_train) 
#disturbance density 2000 m
ew_glm_17d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(distden2015_2000), family=binomial,data=ew_train) 
#disturbance density 3000 m
ew_glm_18d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(distden2015_3000), family=binomial,data=ew_train) 
#line density 250 m
ew_glm_19d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(lineden2015_250), family=binomial,data=ew_train) 
#line density 500 m
ew_glm_20d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(lineden2015_500), family=binomial,data=ew_train) 
#line density 1000 m
ew_glm_21d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(lineden2015_1000), family=binomial,data=ew_train) 
#line density 2000 m
ew_glm_22d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(lineden2015_2000), family=binomial,data=ew_train) 
#line density 3000 m
ew_glm_23d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(lineden2015_3000), family=binomial,data=ew_train) 
#null
ew_glm_24d <- glm(presence~1,family=binomial,data=ew_train) 
#mine 3000 m buff and road 500 m buff
ew_glm_25d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_3000) + scale(road2015_500), family=binomial,data=ew_train)
#mine 3000 m buff and distden 2000 
ew_glm_26d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_3000) + scale(distden2015_2000), family=binomial,data=ew_train)
#mine 3000 m buff and road 500 m buff and distden 2000
ew_glm_27d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015 + scale(lichen_bi100) + scale(conifer_bi250) + scale(gram_bi) + scale(fire2015) + scale(mine2015_3000) + scale(road2015_500) + scale(distden2015_2000), family=binomial,data=ew_train)


ew_d_glms <- list(ew_glm_1d, ew_glm_2d, ew_glm_3d, ew_glm_4d, ew_glm_5d, ew_glm_6d, ew_glm_8d, ew_glm_9d, ew_glm_10d, ew_glm_12d, ew_glm_14d, ew_glm_15d, ew_glm_16d, ew_glm_17d, ew_glm_18d, ew_glm_19d, ew_glm_20d, ew_glm_21d, ew_glm_22d, ew_glm_23d, ew_glm_24d, ew_glm_25d, ew_glm_26d, ew_glm_27d, ew_glm_1e)
ew_d_glms_names <- c('ew_glm_1d', 'ew_glm_2d', 'ew_glm_3d', 'ew_glm_4d', 'ew_glm_5d', 'ew_glm_6d', 'ew_glm_8d', 'ew_glm_9d', 'ew_glm_10d', 'ew_glm_12d', 'ew_glm_14d', 'ew_glm_15d', 'ew_glm_16d', 'ew_glm_17d', 'ew_glm_18d', 'ew_glm_19d', 'ew_glm_20d', 'ew_glm_21d', 'ew_glm_22d', 'ew_glm_23d', 'ew_glm_24d', 'ew_glm_25d', 'ew_glm_26d', 'ew_glm_27d', 'ew_glm_1e')
aictab(cand.set = ew_d_glms, modnames = ew_d_glms_names)

#mining 3000 m buffer is top disturbance variable, followed by mining 2000 m buff and mining 500 m buff
#road 500 m buff best buffered road variable 
#lineden 3000 m best density variable
#distden 2000 best disturbance density variable

#mine 3000 m and road 500 m can go in same model
#mine 3000 m and lineden 3000 are correlated
#mine 3000 m and distden 2000 not correlated, can go in same model
#road 500 and lineden 3000 not correlated, can go in same model
#road 500 and distden 2000 not correlated, can go in same model
#lineden 3000 and distden 2000 are correlated
summary(ew_glm_9d) 
summary(ew_glm_25d)


# Model Calibration
## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(ew_glm_25d, newdata = ew_test, type = "response")

# Step 3: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = ew_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Bin the predicted probabilities using deciles.
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()
#good calibration, the model performs well. Over and under prediction seems slightly increased from base environmental model but AUC score (see next code) is higher.

#Test the model / model evaluation
#top model with unscaled predictors
ew_glm_25d_us <- glm(presence~kdem + I(kdem^2) + landcover2015 + lichen_bi100 + conifer_bi250 + gram_bi + fire2015 + mine2015_3000 + road2015_500, family=binomial,data=ew_train)
summary(ew_glm_25d_us)
#evaluate
ew_glm25d_eval <- pa_evaluate(predict(ew_glm_25d_us, ew_test[ew_test$presence==1, ]), predict(ew_glm_25d_us, ew_test[ew_test$presence==0, ]))
print(ew_glm25d_eval)
plot(ew_glm25d_eval, "ROC")
#AUC=0.95 excellent discrimination.
#################
#Model performance with easystats
r2(ew_glm_25d_us) #Tjur's R2 = 0.525

windows()
check_model(ew_glm_25d_us)
#effect size plot with see pkg from easystats
plot(effectsize(ew_glm_25d_us)) +
  ggplot2::labs(title = "Klaza - Early Winter 2015") +
  ggplot2::scale_y_discrete(labels = c(
    "kdem" = "Elevation",
    "I(kdem^2)" = "Elevation^2",
    "landcover20155" = "Deciduous/mixed 30 m",
    "landcover20158" = "Shrublands 30 m",
    "landcover201510" = "Grasslands 30 m",
    "landcover201514" = "Wetlands 30 m",
    "landcover201516" = "Non-vegetated 30 m",
    "lichen_bi100" = "Lichen 100 m",
    "conifer_bi250" = "Conifer 250 m",
    "gram_bi" = "Graminoid 30 m",
    "fire2015" = "Burns ≤ 50 yrs",
    "mine2015_3000" = "Mining 3000 m",
    "road2015_500" = "Roads Trails 500 m"
  )) 

plot(parameters(ew_glm_25d_us)) +
  ggplot2::labs(title = "Early Winter")

#Save models
saveRDS(ew_glm_25d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_ew_scaled.rds")
saveRDS(ew_glm_25d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_ew_unscaled.rds")

#Try out predict function 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast2 <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kfirerast2015_1.tif')
minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kminebuffrast2015_1.tif')
roadbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kroadbuffrast2015_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif')%>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast2$kdem) #dem to numeric - should have been done at raster prep script
landcover2015 <- as.factor(envrast2$landcover2015) 
lichen_bi100 <- rast100m$lichen_bi100 |>
  resample(mask, method='near')
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
gram_bi <- envrast$gram_bi
fire2015 <- firerast$fire2015
mine2015_3000 <- minebuffrast$mine2015_3000
road2015_500 <- roadbuffrast$road2015_500

rasters <- c(kdem, landcover2015, lichen_bi100, conifer_bi250, gram_bi, fire2015, mine2015_3000, road2015_500) #combine rasters. Make sure they are separated by commas.
#run predict function on disturbance model
p1 <- predict(rasters, ew_glm_25d_us, type="response")
plot(p1) #plots!

writeRaster(p1,'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_ew_mod25d_prediction.tif' , overwrite=TRUE)

#Add gps points to map
#filtered used points w/o covariates extracted to them:
ew_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_ew_daily.shp')
str(ew_points)
st_crs(ew_points)

#Plot the raster and add points on top
plot(st_geometry(ew_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

####Calculate peak of quadratic relationship for elevation

summary(ew_glm_25d_us) #get unscaled coefficients for kdem & kdem2
#kdem = 0.05136, kdem2 = -1.871*10^-5
ew_dem_peak = -(0.05136 / (2 * (-1.871*10^-5)))
ew_dem_peak # peak occurs at 1372.528 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -36.50          # Intercept
beta_1 <- 0.05136         # Coefficient for kdem
beta_2 <- -1.871e-05      # Coefficient for I(kdem^2)

# Create a sequence of kdem values around the observed range
kdem_values <- seq(500, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * kdem_values + beta_2 * kdem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(kdem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Early Winter")
abline(v = 1372.53, col = "red", lty = 2)  # Add vertical line at the peak
