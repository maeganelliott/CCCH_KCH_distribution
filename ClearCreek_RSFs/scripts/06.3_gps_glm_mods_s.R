#script to create seasonal RSF models for clear creek gps data
## SUMMER
#DeMars et al. (2020) used the glmmTMB function from the glmmTMB package
#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

library(easystats)
library(see)
library(effectsize)
library(parameters)
library(DHARMa)
library(patchwork)
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
library(ggplot2)
library(car)
library(AICcmodavg)

set.seed(123)

#Import used and available points with covariates previously extracted to them
s_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_s_data_core.csv')

## SUMMER

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(s_data$landcover2015_100m)
s_data$landcover2015_100m <- factor(s_data$landcover2015_100m)
levels(s_data$landcover2015_100m)
# Relabel the levels by first creating a vector of the current levels. Make sure to use correct resolution of LC for the season.
landcover_levels <- levels(s_data$landcover2015_100m)
# Modify specific levels by direct replacement
landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
# Assign the modified levels back to the factor
levels(s_data$landcover2015_100m) <- landcover_levels
# Verify the new levels
levels(s_data$landcover2015_100m)

str(s_data$ccfiredec)
s_data$ccfiredec <- factor(s_data$ccfiredec)
levels(s_data$ccfiredec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(s_data), 0.8*nrow(s_data))
s_train = s_data[n,]
s_test = s_data[-n,] 

#############################################################################

#Run univariate model selection on covariates that were unclear from above evaluations
#Landcover
lc30_1 <- glm(presence~cclandcover2015, family=binomial, data=s_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=s_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=s_train)

summary(lc100_1) #LC 100 is the best based on AIC model selection. large, positive z scores.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Lichen
lichbi30 <- glm(presence~cclichen_bi, family=binomial, data=s_train)
lichenbi100 <- glm(presence~cclichen_bi100, family=binomial, data=s_train)
lichenbi250 <- glm(presence~cclichen_bi250, family=binomial, data=s_train)
summary(lichbi30) # Lichen binomial 30 m is the best based on AIC model selection. High, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~ccconifer_bi, family=binomial, data=s_train)
conbi100 <- glm(presence~ccconifer_bi100, family=binomial, data=s_train)
conbi250 <- glm(presence~ccconifer_bi250, family=binomial, data=s_train)
summary(conbi100) #conbi100 is best based on AIC scores. High, negative z score.
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~ccdecid_bi, family=binomial, data=s_train)
decidbi100 <- glm(presence~ccdecid_bi100, family=binomial, data=s_train)
decidbi250 <- glm(presence~ccdecid_bi250, family=binomial, data=s_train)
summary(decid250_1) #deciduous bi 30 m is the best based on AIC model selection. High, negative z score. 
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~ccgram_bi, family=binomial, data=s_train)
grambi100 <- glm(presence~ccgram_bi100, family=binomial, data=s_train)
grambi250 <- glm(presence~ccgram_bi250, family=binomial, data=s_train)
summary(grambi30) #graminoid binomial 30 m is the best based on AIC scores. High, positive z score. Gram 30% is best percent cover covariate.
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)


#Include: dem (q), slope (q), LC100, Lichenbi30, coniferbi100, decidbi30, grambi30

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Summer:
#add quadratic terms
s_data$ccslope2 <- s_data$ccslope^2 #quadratic transformation
s_data$ccdem2 <- s_data$ccdem^2

#missing_cols <- setdiff(c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m", "cclichen_bi", "ccconifer_bi100", "ccdecid_bi", "ccgram_bi", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_100", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000"), colnames(s_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.


s_dat.cor = s_data[,c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m", "cclichen_bi", "ccconifer_bi100", "ccdecid_bi", "ccgram_bi", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")]
names(s_dat.cor) = c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m", "cclichen_bi", "ccconifer_bi100", "ccdecid_bi", "ccgram_bi", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")
cor_matrix <- cor(s_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(s_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/s_correlation_matrix_core.csv", row.names = TRUE)


############################################################################
#Summer GLMs
#Some experimentation  ****************when re-running: elevation and conifer are negatively correlated (-0.7). Proceed with elevation only.
#exclude distance to water
#exclude fire (summer)
#core data: dem (q), slope (q), LC100, Lichenbi30, coniferbi100, decidbi30, grambi30
#all variables 
s_test1 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train)
#all variables, no landcover
s_test2 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train)
#removed binomial variables
s_test3 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m, family=binomial,data=s_train)
#all variables, no decid
s_test4 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)


s_test_glms <- list(s_test1, s_test2, s_test3, s_test4)
s_test_glms_names <- c('s_test1', 's_test2', 's_test3', 's_test4')
aictab(cand.set = s_test_glms, modnames = s_test_glms_names)

summary(s_test4)
#elevation is primary driver of model
#most landcover classes are not significant. 
#in full model landcover deciduous forest (sig) and shrublands have negative effects, grassland has sig positive effect, wetlands and non-vegeetated negative effects
#in model without binomial vars landcover classes are all significant and shrublands changed to a positive effect. 
#deciduous shrub weakest predictor
#some correlation between conifer cover and landcover (-0.49) but below cutoff.

#Unique model for each season
# Phase 1: Find the best environmental base model
#INCLUDE: dem (q), slope (q), LC100, Lichenbi30, coniferbi100, decidbi30, grambi30
#Starting with glm (no mixed effect), standardized predictors:

#All top performing predictors from evaluations
s_glm_1e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing predictors from evaluations
#All top performing variables from evaluation, no graminoid
s_glm_2e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi), family=binomial,data=s_train) # all top performing variable,s no graminoid
#All top performing variables from evaluations, no deciduous shrub
s_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no deciduous shrub
#All top performing variables from evaluations, no conifer
s_glm_4e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no conifer
#All top performing variables from evaluations, no lichen
s_glm_5e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no lichen                 
#All top performing variables from evaluations, no landcover
s_glm_6e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no landcover
#All top performing variables from evaluations, no slope
s_glm_7e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no slope
#all top performing variables from evaluations, no elevation
s_glm_8e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi) + scale(ccgram_bi), family=binomial,data=s_train) #all top performing variables from evaluations, no elevation
#no landcover, no decid
s_glm_9e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#null model
s_glm_10e <- glm(presence~1,family=binomial,data=s_train) #null model

s_e_glms <- list(s_glm_1e, s_glm_2e, s_glm_3e, s_glm_4e, s_glm_5e, s_glm_6e, s_glm_7e, s_glm_8e, s_glm_9e, s_glm_10e)
s_e_glms_names <- c('s_glm_1e', 's_glm_2e', 's_glm_3e', 's_glm_4e', 's_glm_5e', 's_glm_6e', 's_glm_7e', 's_glm_8e', 's_glm_9e', 's_glm_10e')
aictab(cand.set = s_e_glms, modnames = s_e_glms_names)
#Models without fire performed the best. Fire was a close second but was non-significant (previous version)
#distance to water second least important variable but still improves model fit (previous version)
#elevation, conifer, deciduous shrub, lichen all important variables. 
#model 1e performed the best - all variables 
summary(s_glm_3e)
summary(s_glm_1e)
summary()

report(lw_glm_4e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(s_glm_3e, newdata = s_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = s_test$presence,       # Binary response variable in the test data
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
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(s_glm_3e, newdata = s_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = s_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#Good calibration. model performs reasonably well. Predicted probabilities generally align with the observed rates across the 10 bins. 
#the model slightly under-predicts in the lower range (0.05 - 0.25) and slightly over-predicts for higher probability bin (>0.5)

#Evaluation with AUC / ROC
#unscaled model:
s_glm_3e_us <- glm(presence~ccdem + I(ccdem^2) + ccslope + I(ccslope^2) + landcover2015_100m + cclichen_bi + ccconifer_bi100 + ccgram_bi, family=binomial,data=s_train)
#evaluate
s_glm_3e_eval <- pa_evaluate(predict(s_glm_3e_us, s_test[s_test$presence==1, ]), predict(s_glm_3e_us, s_test[s_test$presence==0, ]))
print(s_glm_3e_eval)
plot(s_glm_3e_eval, "ROC")
#AUC=0.944 #high AUC score, excellent model performance.
#########################
#Model performance with easystats
r2(s_glm_3e_us) #Tjur's R2 = 0.533

#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m2.tif') %>%
  #resample(ccmask, method="near")
landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 250m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16) 
cclichen_bi <- envrast$cclichen_bi
ccconifer_bi100 <- rast100m$ccconifer_bi100 |>
  resample(ccmask, method='near')
ccgram_bi <- envrast$ccgram_bi

rasters <- c(ccdem, ccslope, landcover2015_100m, cclichen_bi, ccconifer_bi100, ccgram_bi) 
#run predict function
p1 <- predict(rasters, s_glm_3e_us, type="response")
plot(p1) 

#Add gps points to map
#filtered used points w/o covariates extracted to them:
s_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/GPS_collar_data/cc_gps_2018_summer_core.csv')
st_crs(s_points)
# convert from lat/long to UTM to match rasters
s_points <- st_as_sf(s_points, coords = c("location_long", "location_lat"), crs = 4326)
s_points <- st_transform(s_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(s_points)
plot(st_geometry(s_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
## Summer Disturbance Model Selection with Interaction Terms
#only test the variables that make sense based on previous evaluation. Excluding the thresholds where distributions switched.
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 

#model 5e + disturbance 30m buff
s_glm_1d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_30), family=binomial,data=s_train)
#model 5e + disturbance 250m buff
s_glm_2d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_250), family=binomial,data=s_train)
#model 5e + disturbance 500m buff
s_glm_3d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_500), family=binomial,data=s_train)
#model 5e + disturbance 1000m buff
s_glm_4d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_1000), family=binomial,data=s_train)
#dist 1000 m buff + interaction with elevation
s_glm_4di <- glm(presence~scale(ccdem) * scale(dist2018_1000) + scale(I(ccdem^2)) * scale(dist2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#dist 2000 m buff
s_glm_5d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_2000), family=binomial,data=s_train)
#model 5e + disturbance 2000m buff + interaction with elevation
s_glm_5di <- glm(presence~scale(ccdem) * scale(dist2018_2000) + scale(I(ccdem^2)) * scale(dist2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#dist 3000 m buff
s_glm_6d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(dist2018_3000), family=binomial,data=s_train)
#model 5e + disturbance 3000m buff + interaction with elevation
s_glm_6di <- glm(presence~scale(ccdem) * scale(dist2018_3000) + scale(I(ccdem^2)) * scale(dist2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#model 5e + mining 30 m buff - exclude all mining because no relationship apparent in distributions. 
#s_glm_7d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_30), family=binomial,data=s_train)
#model 5e + mining 250 m buff
#s_glm_8d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_250), family=binomial,data=s_train)
#model 5e + mining 500 m buff
#s_glm_9d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_500), family=binomial,data=s_train)
#mining 1000 m buff
#s_glm_10d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_1000), family=binomial,data=s_train)
#mining 2000 m buff
#s_glm_11d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_2000), family=binomial,data=s_train)
#mining 3000 m buff
#s_glm_12d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(mine2018_3000), family=binomial,data=s_train)
#model 5e + roads and trails 30 m buff
s_glm_13d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(road2018_30), family=binomial,data=s_train)
#model 5e + roads and trails 250 m buff
s_glm_14d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(road2018_250), family=binomial,data=s_train)
#model 5e + roads and trails 500 m buff
s_glm_15d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(road2018_500), family=binomial,data=s_train)
#model 5e + roads and trails 1000 m buff
s_glm_16d <- glm(presence~scale(ccdem) + scale(I(ccdem^2))+ scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(road2018_1000), family=binomial,data=s_train)
#model 5e + roads and trails 1000 m buff + interaction with dem
s_glm_16di <- glm(presence~scale(ccdem) * scale(road2018_1000) + scale(I(ccdem^2)) * scale(road2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#model 5e + disturbance density 250 m radii
s_glm_17d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(distden2018_250), family=binomial,data=s_train)
#model 5e + disturbance density 500 m radii
s_glm_18d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(distden2018_500), family=binomial,data=s_train)
#model 5e + disturbance density 1000 m radii
s_glm_19d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(distden2018_1000), family=binomial,data=s_train)
#model 5e + disturbance density 2000 m radii
s_glm_20d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(distden2018_2000), family=binomial,data=s_train)
#model 5e + disturbance density 3000 m radii
s_glm_21d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(distden2018_3000), family=binomial,data=s_train)
#model 5e + LF density 250 m radii
s_glm_22d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(lineden2018_250), family=binomial,data=s_train)
#model 5e + LF density 500 m radii
s_glm_23d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(lineden2018_500), family=binomial,data=s_train)
#model 5e + LF density 1000 m radii
s_glm_24d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(lineden2018_1000), family=binomial,data=s_train)
#lineden 2000 m radii
s_glm_25d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(lineden2018_2000), family=binomial,data=s_train)
#model 5e + LF density 2000 m radii + interaction with elevation
s_glm_25di <- glm(presence~scale(ccdem) * scale(lineden2018_2000) + scale(I(ccdem^2)) * scale(lineden2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#lineden 3000 m radii
s_glm_26d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi) + scale(lineden2018_3000), family=binomial,data=s_train)
#model 5e + LF density 3000 m radii + interaction with elevation
s_glm_26di <- glm(presence~scale(ccdem) * scale(lineden2018_3000) + scale(I(ccdem^2)) * scale(lineden2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#model 3e
s_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccgram_bi), family=binomial,data=s_train)
#null model
s_glm_27d <- glm(presence~1,family=binomial,data=s_train)
#dist 3000m buff + distden 3000
#s_glm_28d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi250) + scale(ccgram_bi) + scale(dist2018_3000) + scale(distden2018_3000), family=binomial,data=s_train)
#distden 3000 + mine 3000m buff
#s_glm_29d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi100) + scale(ccdecid_bi250) + scale(ccgram_bi) + scale(distden2018_3000) + scale(mine2018_3000), family=binomial,data=s_train)

#AIC
s_d_glms <- list(s_glm_1d, s_glm_2d, s_glm_3d, s_glm_4d, s_glm_4di, s_glm_5d, s_glm_5di, s_glm_6d, s_glm_6di, s_glm_13d, s_glm_14d, s_glm_15d, s_glm_16d, s_glm_16di, s_glm_17d, s_glm_18d, s_glm_19d, s_glm_20d, s_glm_21d, s_glm_22d, s_glm_23d, s_glm_24d, s_glm_25d, s_glm_25di, s_glm_26d, s_glm_26di, s_glm_27d, s_glm_3e)
s_d_glm_names <- c('s_glm_1d', 's_glm_2d', 's_glm_3d', 's_glm_4d', 's_glm_4di', 's_glm_5d', 's_glm_5di', 's_glm_6d', 's_glm_6di', 's_glm_13d', 's_glm_14d', 's_glm_15d', 's_glm_16d', 's_glm_16di', 's_glm_17d', 's_glm_18d', 's_glm_19d', 's_glm_20d', 's_glm_21d', 's_glm_22d', 's_glm_23d', 's_glm_24d', 's_glm_25d', 's_glm_25di', 's_glm_26d', 's_glm_26di', 's_glm_27d', 's_glm_3e')
aictab(cand.set = s_d_glms, modnames = s_d_glm_names)

summary(s_glm_16di)
summary(s_glm_26di)
summary(s_glm_21d)

report()

#some correlation between elevation and some disturbance variables(-0.3 to -0.4) that caused effect directions and significance to change. Tested an interaction term for these models.
#Proceed with 16di, with roads and trails with 1000 m buffer

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(s_glm_16di, newdata = s_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = s_test$presence,       # Binary response variable in the test data
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
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(s_glm_16di, newdata = s_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = s_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#Good calibration. model performs well. Predicted probabilities generally align with the observed rates across the 10 bins. 
#the model slightly under-predicts in the lower range (0.05 - 0.25) and slightly over-predicts for higher probability bin (>0.5)

#Evaluation with AUC / ROC
#unscaled model:
s_glm_16di_us <- glm(presence~ccdem * road2018_1000 + I(ccdem^2) * road2018_1000 + ccslope + I(ccslope^2) + landcover2015_100m + cclichen_bi + ccconifer_bi100 + ccgram_bi, family=binomial,data=s_train)
#evaluate
s_glm_16di_eval <- pa_evaluate(predict(s_glm_16di_us, s_test[s_test$presence==1, ]), predict(s_glm_16di_us, s_test[s_test$presence==0, ]))
print(s_glm_16di_eval)
plot(s_glm_16di_eval, "ROC")
#AUC=0.944 #high AUC score. Slight improvement from base environmental model. 
###################################
#Model performance with easystats
r2(s_glm_16di_us) #Tjur's R2 = 0.536

windows()
check_model(s_glm_16di_us)
#effect size plot with see pkg from easystats
plot(effectsize(s_glm_16di_us)) +
  ggplot2::labs(title = "Clear Creek - Summer 2018") +
  ggplot2::scale_y_discrete(labels = c(
    "ccdem" = "Elevation",
    "road2018_1000" = "Roads/Trails 1000 m",
    "I(ccdem^2)" = "Elevation^2",
    "ccslope" = "Slope",
    "I(ccslope^2)" = "Slope^2",
    "landcover2015_100m5" = "Deciduous/mixed 100 m",
    "landcover2015_100m8" = "Shrublands 100 m",
    "landcover2015_100m10" = "Grasslands 100 m",
    "landcover2015_100m14" = "Wetlands 100 m",
    "landcover2015_100m16" = "Non-vegetated 100 m",
    "cclichen_bi" = "Lichen 30 m",
    "ccconifer_bi100" = "Conifer 100 m",
    "ccgram_bi" = "Graminoid 30 m",
    "ccdem:road2018_1000" = "Elevation:Roads/Trails 1000 m",
    "road2018_1000:I(ccdem^2)" = "Elevation^2:Roads/Trails 1000 m"
  )) 

plot(parameters(s_glm_16di_us)) +
  ggplot2::labs(title = "Summer")

#Save models
saveRDS(s_glm_16di, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_s_scaled.rds")
saveRDS(s_glm_16di_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_s_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
#ndvirast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccndvi.tif')
#firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccfirerast2018_1.tif')
#distbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistbuffrast2018_1.tif')
#distdenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistdenrast2018_core.tif')
#minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccminebuffrast2018_1.tif')
roadbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccroadbuffrast2018_core.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m2.tif') %>%
# resample(ccmask, method="near")

#levels(s_data$landcover2015_100m)
#landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
#landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
#landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
#levels(landcover2015_100m)

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 2100m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16)
ccslope <- envrast$ccslope
cclichen_bi <- envrast$cclichen_bi
ccconifer_bi100 <- rast100m$ccconifer_bi100 |>
  resample(ccmask, method='near')
ccgram_bi <- envrast$ccgram_bi
#distden2018_3000 <- distdenrast$distden2018_3000
road2018_1000 <- roadbuffrast$road2018_1000

rasters <- c(ccdem, ccslope, landcover2015_100m, cclichen_bi, ccconifer_bi100, ccgram_bi, road2018_1000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, s_glm_16di_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_s_mod16di_prediction.tif')
plot(st_geometry(s_points), add = TRUE, col = "red", pch = 19, cex = 0.6) #add all used points

####################################################
####Calculate peak of quadratic relationship for elevation
windows()
summary(s_glm_16di_us) #get unscaled coefficients for elevation
#dem = 0.03064, dem2 = -8.733*10^-6
s_dem_peak = -(0.03064 / (2 * (-8.733*10^-6)))
s_dem_peak # peak occurs at 1754.265 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -27.87          # Intercept
beta_1 <- 0.03064         # Coefficient for dem
beta_2 <- -8.733*10^-6      # Coefficient for I(dem^2)

# Create a sequence of dem values around the observed range
dem_values <- seq(1000, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * dem_values + beta_2 * dem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(dem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Summer")
abline(v = 1754.265, col = "red", lty = 2)  # Add vertical line at the peak
