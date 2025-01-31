#script to create seasonal RSF models for clear creek gps data
## FALL RUT
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
f_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_f_data_core.csv')
str(f_data)
## FALL RUT

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(f_data$landcover2015_100m)
f_data$landcover2015_100m <- factor(f_data$landcover2015_100m)
levels(f_data$landcover2015_100m)
# Relabel the levels by first creating a vector of the current levels. Make sure to use correct resolution of LC for the season.
landcover_levels <- levels(f_data$landcover2015_100m)
# Modify specific levels by direct replacement
landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
# Assign the modified levels back to the factor
levels(f_data$landcover2015_100m) <- landcover_levels
# Verify the new levels
levels(f_data$landcover2015_100m)

str(f_data$ccfiredec)
f_data$ccfiredec <- factor(f_data$ccfiredec)
levels(f_data$ccfiredec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(f_data), 0.8*nrow(f_data))
f_train = f_data[n,]
f_test = f_data[-n,] 

#######################################################################
#Univariate Model Selection
#Landcover
#Landcover
lc30_1 <- glm(presence~cclandcover2015, family=binomial, data=f_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=f_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=f_train)
#LC 100 is the best based on AIC model selection. large, negative z scores.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Lichen
lichbi30 <- glm(presence~cclichen_bi, family=binomial, data=f_train)
lichenbi100 <- glm(presence~cclichen_bi100, family=binomial, data=f_train)
lichenbi250 <- glm(presence~cclichen_bi250, family=binomial, data=f_train)
#summary() # Lichen percent cover 30 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~ccconifer_bi, family=binomial, data=f_train)
conbi100 <- glm(presence~ccconifer_bi100, family=binomial, data=f_train)
conbi250 <- glm(presence~ccconifer_bi250, family=binomial, data=f_train)
#summary(conbi30) #conbi30 is best based on AIC scores. High, negative z score.
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~ccdecid_bi, family=binomial, data=f_train)
decidbi100 <- glm(presence~ccdecid_bi100, family=binomial, data=f_train)
decidbi250 <- glm(presence~ccdecid_bi250, family=binomial, data=f_train)
summary(decidbi250) #deciduous binomial 100 m is the best based on AIC model selection. High, negative z score. 
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~ccgram_bi, family=binomial, data=f_train)
grambi100 <- glm(presence~ccgram_bi100, family=binomial, data=f_train)
grambi250 <- glm(presence~ccgram_bi250, family=binomial, data=f_train)
summary(grambi30) #graminoid binomial 100 m is the best based on AIC scores. High, positive z score. Gram 30% is best percent cover covariate.
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)

#Inclue: dem (q), slope (q), LC100, Lichenbi 30m, Coniferbi 30m, Decidbi 100m, Grambi 100m

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Fall Rut:
#add quadratic terms
f_data$ccslope2 <- f_data$ccslope^2 #quadratic transformation
f_data$ccdem2 <- f_data$ccdem^2

#missing_cols <- setdiff(c("ccdem", "ccdem2", "ccslope", "ccslope2", "ccrough", "ndvi_f", "waterd", "waterd2", "landcover2015_250m", "cclichen_bi100", "ccconifer_bi", "ccdecid_bi250", "ccgram_bi", "cceverg_bi", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018d", "dist2018d2", "dist2018_250", "dist2018_500", "dist2018_100", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018d", "mine2018d2", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018d", "road2018d2", "road2018_250", "road2018_500", "road2018_1000"), colnames(f_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.

f_dat.cor = f_data[,c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m", "cclichen_bi", "ccconifer_bi", "ccdecid_bi100", "ccgram_bi100", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")]
names(f_dat.cor) = c("ccdem", "ccdem2", "ccslope", 'ccslope2', "landcover2015_100m", "cclichen_bi", "ccconifer_bi", "ccdecid_bi100", "ccgram_bi100", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")
cor_matrix <- cor(f_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(f_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/f_correlation_matrix_core.csv", row.names = TRUE)
#slope and roughness are highly correlated
#the fire variables are highly correlated with each other
#disturbance density, LF density highly correlated
#buffered disturbance and roads highly correlated
#distance to variables highly correlated.

#####################################################################
#Fall Rut GLMs
#Some experimentation
#exclude fire and waterd
#core data: dem (q), slope (q), LC100, Lichenbi 30m, Coniferbi 30m, Decidbi 100m, Grambi 100m
#conifer and elevation are negatively correlated and will need to go in separate models
#elevation mods
#all variables
f_test1 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#no binomial veg
f_test2 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m, family=binomial,data=f_train)
#no landcover
f_test3 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#conifer
#all variables
f_test4 <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)

summary(f_test4)

#Elevation is primary driver of model
#lichen is not significant
#shrublands and deciduous shrub have opposite effect directions. Shrublands has larger effect size (positive)
#conifer cover has significant negative effect. Smaller effect size than elevation. 
#lichen is significant in conifer model and both slope terms are significant

#Unique model for each season
# Phase 1: Find the best environmental base model
#INCLUDE: dem (q), slope (q), LC100, Lichenbi 30m, Coniferbi 30m, Decidbi 100m, Grambi 100m
#conifer cover and elevation are negatively correlated
#Starting with glm (no mixed effect), standardized predictors:

#elevation
#global model
f_glm_1e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train) #All top performing predictors from evaluations, no fire
#All top performing variables from evaluation, no graminoid
f_glm_2e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi100), family=binomial,data=f_train) # all top performing variable,s no graminoid
#All top performing variables from evaluations, no deciduous shrub
f_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi100), family=binomial,data=f_train) #All top performing variables from evaluations, no deciduous shrub
#All top performing variables from evaluations, no lichen
f_glm_4e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train) #All top performing variables from evaluations, no lichen                 
#All top performing variables from evaluations, no landcover
f_glm_5e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train) #All top performing variables from evaluations, no landcover
#All top performing variables from evaluations, no slope
f_glm_6e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train) #All top performing variables from evaluations, no slope
#all top performing variables from evaluations, no elevation or conifer
f_glm_7e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train) #all top performing variables from evaluations, no elevation
#conifer 
#global
f_glm_8e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#no gram
f_glm_9e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccdecid_bi100), family=binomial,data=f_train)
#no decid
f_glm_10e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi100), family=binomial,data=f_train)
#no lichen
f_glm_11e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccconifer_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#no landcover
f_glm_12e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#no slope
f_glm_13e <- glm(presence~landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#null model
f_glm_14e <- glm(presence~1,family=binomial,data=f_train) #null model
#no conifer, no waterd
#f_glm_14e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi100) + scale(ccdecid_bi250) + scale(ccgram_bi) + ccfiredec, family=binomial,data=f_train) #no conifer, no waterd
#no conifer, no waterd, no lichen
#f_glm_15e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(ccdecid_bi250) + scale(ccgram_bi) + ccfiredec, family=binomial,data=f_train) #no conifer, no water, no graminoid

f_e_glms <- list(f_glm_1e, f_glm_2e, f_glm_3e, f_glm_4e, f_glm_5e, f_glm_6e, f_glm_7e, f_glm_8e, f_glm_9e, f_glm_10e, f_glm_11e, f_glm_12e, f_glm_13e, f_glm_14e)
f_e_glms_names <- c('f_glm_1e', 'f_glm_2e', 'f_glm_3e', 'f_glm_4e', 'f_glm_5e', 'f_glm_6e', 'f_glm_7e', 'f_glm_8e', 'f_glm_9e', 'f_glm_10e', 'f_glm_11e', 'f_glm_12e', 'f_glm_13e', 'f_glm_14e')
aictab(cand.set = f_e_glms, modnames = f_e_glms_names)
#removing lichen improves model fit for elevation model
#elevation is the top predictor
#some low levels of correlation between elevation and some of the other predictor variables (slope, lichen)
#elevation models performed better than the conifer models
#of the conifer models, the global one fit the best. 
#graminoid is important for model fit

summary(f_glm_4e)
summary()

report(lw_glm_5e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(f_glm_4e, newdata = f_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = f_test$presence,       # Binary response variable in the test data
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
predicted_probs_test <- predict(f_glm_4e, newdata = f_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = f_test$presence,       # Binary response variable
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

#Fair calibration. model performs reasonably well, some deviation. Deviation increases >0.75 predicted probabilities
#Predicted probabilities generally align with the observed rates across the 10 bins. 
#the model tends to slightly under-predict

#Evaluation with AUC / ROC
#unscaled model:
f_glm_4e_us <- glm(presence~ccdem + I(ccdem^2) + ccslope + I(ccslope^2) + landcover2015_100m + ccdecid_bi100 + ccgram_bi100, family=binomial,data=f_train)
#evaluate
f_glm_4e_eval <- pa_evaluate(predict(f_glm_4e_us, f_test[f_test$presence==1, ]), predict(f_glm_4e_us, f_test[f_test$presence==0, ]))
print(f_glm_4e_eval)
plot(f_glm_4e_eval, "ROC")
#AUC=0.918 #high AUC score
##Model performance with easystats
r2(f_glm_4e_us) #Tjur's R2 = 0.489

#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m2.tif') %>%
 # resample(ccmask, method="near")

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 2100m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16)
#cclichen_bi100 <- rast100m$cclichen_bi100 |>
  #resample(ccmask, method='near')
ccdecid_bi100 <- rast100m$ccdecid_bi100 |>
  resample(ccmask, method='near')
ccgram_bi100 <- rast100m$ccgram_bi100 |>
  resample(ccmask, method='near')

rasters <- c(ccdem, ccslope, landcover2015_100m, ccdecid_bi100, ccgram_bi100) 
#run predict function
p1 <- predict(rasters, f_glm_4e_us, type="response")
plot(p1) #plots!

#Add gps points to map
#filtered used points w/o covariates extracted to them:
f_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/GPS_collar_data/CC_gps_2018_fallrut_core.csv')
st_crs(f_points)
f_points <- st_as_sf(f_points, coords = c("location_long", "location_lat"), crs = 4326)
f_points <- st_transform(f_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(f_points)
plot(st_geometry(f_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation. 
#Try out interaction term for variables that are correlated with elevation
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 
#model 7e + disturbance 30m buff
f_glm_1d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_30), family=binomial,data=f_train)
#model 7e + disturbance 250m buff
f_glm_2d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_250), family=binomial,data=f_train)
#model 7e + disturbance 500m buff
f_glm_3d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_500), family=binomial,data=f_train)
#model 7e + disturbance 1000m buff
f_glm_4d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_1000), family=binomial,data=f_train)
#dist 1000m buff + interaction with elevation
f_glm_4di <- glm(presence~scale(ccdem) * scale(dist2018_1000) + scale(I(ccdem^2)) * scale(dist2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e + disturbance 2000m buff
f_glm_5d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_2000), family=binomial,data=f_train)
#dist 2000 m buff + interaction with elevation
f_glm_5di <- glm(presence~scale(ccdem) * scale(dist2018_2000) + scale(I(ccdem^2)) * scale(dist2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e + disturbance 3000m buff
f_glm_6d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(dist2018_3000), family=binomial,data=f_train)
#dist 3000 m buff + interaction with elevation
f_glm_6di <- glm(presence~scale(ccdem) * scale(dist2018_3000) + scale(I(ccdem^2)) * scale(dist2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e + mining 30 m buff
f_glm_7d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_30), family=binomial,data=f_train)
#model 7e + mining 250 m buff
f_glm_8d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_250), family=binomial,data=f_train)
#model 7e + mining 500 m buff
f_glm_9d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_500), family=binomial,data=f_train)
#model 7e + mining 1000 m buff
f_glm_10d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_1000), family=binomial,data=f_train)
#model 7e + mining 2000 m buff
f_glm_11d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_2000), family=binomial,data=f_train)
#model 7e + mining 3000 m buff
f_glm_12d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_3000), family=binomial,data=f_train)
#model 7e + roads and trails 30 m buff
f_glm_13d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(road2018_30), family=binomial,data=f_train)
#model 7e + roads and trails 250 m buff
f_glm_14d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(road2018_250), family=binomial,data=f_train)
#model 7e + roads and trails 500 m buff
f_glm_15d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(road2018_500), family=binomial,data=f_train)
#model 7e + roads and trails 1000 m buff
f_glm_16d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(road2018_1000), family=binomial,data=f_train)
#road 1000 m buff + interaction with elevation
f_glm_16di <- glm(presence~scale(ccdem) * scale(road2018_1000) + scale(I(ccdem^2)) * scale(road2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e + disturbance density 250 m radii
f_glm_17d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(distden2018_250), family=binomial,data=f_train)
#model 7e + disturbance density 500 m radii
f_glm_18d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(distden2018_500), family=binomial,data=f_train)
#model 7e + disturbance density 1000 m radii
f_glm_19d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(distden2018_1000), family=binomial,data=f_train)
#model 7e + disturbance density 2000 m radii
f_glm_20d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(distden2018_2000), family=binomial,data=f_train)
#model 7e + disturbance density 3000 m radii
f_glm_21d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(distden2018_3000), family=binomial,data=f_train)
#model 7e + LF density 250 m radii
f_glm_22d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(lineden2018_250), family=binomial,data=f_train)
#model 7e + LF density 500 m radii
f_glm_23d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(lineden2018_500), family=binomial,data=f_train)
#model 7e + LF density 1000 m radii
f_glm_24d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(lineden2018_1000), family=binomial,data=f_train)
#model 7e + LF density 2000 m radii
f_glm_25d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(lineden2018_2000), family=binomial,data=f_train)
#line den 2000 m + interaction with elevation
f_glm_25di <- glm(presence~scale(ccdem) * scale(lineden2018_2000) + scale(I(ccdem^2)) * scale(lineden2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e + LF density 3000 m radii
f_glm_26d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(lineden2018_3000), family=binomial,data=f_train)
#line den 3000 m + interaction with elevation
f_glm_26di <- glm(presence~scale(ccdem) * scale(lineden2018_3000) + scale(I(ccdem^2)) * scale(lineden2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#model 7e
f_glm_4e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100), family=binomial,data=f_train)
#null model
f_glm_27d <- glm(presence~1,family=binomial,data=f_train)
#distden 1000 + mining 500m buff
f_glm_28d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccdecid_bi100) + scale(ccgram_bi100) + scale(mine2018_500) + scale(distden2018_1000), family=binomial,data=f_train)
#distden 1000 + road 30m buff
#f_glm_29d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi100) + scale(ccdecid_bi250) + scale(ccgram_bi) + scale(distden2018_1000) + scale(road2018_30), family=binomial,data=f_train)
#AIC
f_d_glms <- list(f_glm_1d, f_glm_2d, f_glm_3d, f_glm_4d, f_glm_4di, f_glm_5d, f_glm_5di, f_glm_6d, f_glm_6di, f_glm_7d, f_glm_8d, f_glm_9d, f_glm_10d, f_glm_11d, f_glm_12d, f_glm_13d, f_glm_14d, f_glm_15d, f_glm_16d, f_glm_16di, f_glm_17d, f_glm_18d, f_glm_19d, f_glm_20d, f_glm_21d, f_glm_22d, f_glm_23d, f_glm_24d, f_glm_25d, f_glm_25di, f_glm_26d, f_glm_26di, f_glm_27d, f_glm_28d, f_glm_4e)
f_d_glm_names <- c('f_glm_1d', 'f_glm_2d', 'f_glm_3d', 'f_glm_4d', 'f_glm_4di', 'f_glm_5d', 'f_glm_5di', 'f_glm_6d', 'f_glm_6di', 'f_glm_7d', 'f_glm_8d', 'f_glm_9d', 'f_glm_10d', 'f_glm_11d', 'f_glm_12d', 'f_glm_13d', 'f_glm_14d', 'f_glm_15d', 'f_glm_16d', 'f_glm_16di', 'f_glm_17d', 'f_glm_18d', 'f_glm_19d', 'f_glm_20d', 'f_glm_21d', 'f_glm_22d', 'f_glm_23d', 'f_glm_24d', 'f_glm_25d', 'f_glm_25di', 'f_glm_26d', 'f_glm_26di', 'f_glm_27d', 'f_glm_28d', 'f_glm_4e')
aictab(cand.set = f_d_glms, modnames = f_d_glm_names)

summary(f_glm_9d)

summary(f_glm_26di)
report()

#the only interaction to improve fit better than the base model was lineden 3000 * elevation, but lineden not significant and ccdem:lineden not significant. 
#mining 500 m buff was top disturbance variable. Negative effect but not significant. 
#disturbance density 1000 m was second top disturbance variable, negative effect and just under significance (0.09)
#mining 500 m buff and distden 1000 not correlated so tried in combined model. Better fit based on delta AIC but neither variable significant.
#top models are equivalent based on delta AIC
#Other disturbance variables generally have the expected negative effect based on covariate distributions. None are significant. Some have low levels of correlation with elevation but doesn't seem to affect coefficients.
#only roads and trails 250m buff and distden 2000 and 3000 had (small) positive effects (not significant)
#proceeding with GLM 19d because the disturbance variable is close to significant. GLM 9d (mining 500m buff) has a larger effect size but larger p value. 

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(f_glm_9d, newdata = f_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = f_test$presence,       # Binary response variable in the test data
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
predicted_probs_test <- predict(f_glm_9d, newdata = f_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = f_test$presence,       # Binary response variable
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
#the model tends to slightly under-predict

#Evaluation with AUC / ROC
#unscaled model:
f_glm_9d_us <- glm(presence~ccdem + I(ccdem^2) + ccslope + I(ccslope^2) + landcover2015_100m + ccdecid_bi100 + ccgram_bi100 + mine2018_500, family=binomial,data=f_train)
#save unscaled model
saveRDS(f_glm_9d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_fr_gps_mod.rds")
f_glm_9d_copy <- readRDS(file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_fr_gps_mod.rds")

#evaluate
f_glm_9d_eval <- pa_evaluate(predict(f_glm_9d_us, f_test[f_test$presence==1, ]), predict(f_glm_9d_us, f_test[f_test$presence==0, ]))
print(f_glm_9d_eval)
plot(f_glm_9d_eval, "ROC")
#AUC=0.93 #high AUC score. Improved from base environmental model. 
####################################
#Model performance with easystats
r2(f_glm_9d_us) #Tjur's R2 = 0.491

windows()
check_model(f_glm_9d_us)
#effect size plot with see pkg from easystats
# Rename levels of the landcover2015_100m factor
#levels(f_train$landcover2015_100m) <- c(
 # "Conifer Forest 100 m", "Deciduous/Mixed 100 m", "Shrublands 100 m", "Grasslands 100 m", "Wetlands 100 m", "Non-vegetated 100 m"
)

plot(effectsize(f_glm_9d_us)) +
  ggplot2::labs(title = "Clear Creek - Fall Rut 2018") +
  ggplot2::scale_y_discrete(labels = c(
    "ccdem" = "Elevation",
    "I(ccdem^2)" = "Elevation^2",
    "ccslope" = "Slope",
    "I(ccslope^2)" = "Slope^2",
    "landcover2015_100m5" = "Deciduous/mixed 100 m",
    "landcover2015_100m8" = "Shrublands 100 m",
    "landcover2015_100m10" = "Grasslands 100 m",
    "landcover2015_100m14" = "Wetlands 100 m",
    "landcover2015_100m16" = "Non-vegetated 100 m",
    "ccdecid_bi100" = "Deciduous Shrub 100 m",
    "ccgram_bi100" = "Graminoid 100 m",
    "mine2018_500" = "Mining 500 m"
  )) 

plot(parameters(f_glm_9d_us)) +
  ggplot2::labs(title = "Fall Rut")

#Save models
saveRDS(f_glm_9d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_f_scaled.rds")
saveRDS(f_glm_9d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_f_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
#ndvirast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccndvi.tif')
#firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccfirerast2018_1.tif')
minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccminebuffrast2018_core.tif')
distdenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistdenrast2018_core.tif')
#minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccminebuffrast2018_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m_core.tif') %>%
  #resample(ccmask, method="near")

#landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
#landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
#landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
#levels(landcover2015_100m)
#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 2100m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16)
#cclichen_bi100 <- rast100m$cclichen_bi100 |>
#resample(ccmask, method='near')
ccdecid_bi100 <- rast100m$ccdecid_bi100 |>
  resample(ccmask, method='near')
ccgram_bi100 <- rast100m$ccgram_bi100 |>
  resample(ccmask, method='near')
mine2018_500 <- minebuffrast$mine2018_500
distden2018_1000 <- distdenrast$distden2018_1000

rasters <- c(ccdem, ccslope, landcover2015_100m, ccdecid_bi100, ccgram_bi100, mine2018_500) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, f_glm_9d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_f_mod9d_prediction_core.tif')
plot(st_geometry(f_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

####################################################
####Calculate peak of quadratic relationship for elevation
windows()
summary(f_glm_9d_us) #get unscaled coefficients for elevation
#dem = 0.03683, dem2 = -1.039*10^-5
f_dem_peak = -(0.03683 / (2 * (-1.039*10^-5)))
f_dem_peak # peak occurs at 1772.377 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -33.64          # Intercept
beta_1 <- 0.03683         # Coefficient for dem
beta_2 <- -1.039*10^-5      # Coefficient for I(dem^2)

# Create a sequence of dem values around the observed range
dem_values <- seq(1000, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * dem_values + beta_2 * dem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(dem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Fall Rut")
abline(v = 1772.377, col = "red", lty = 2)  # Add vertical line at the peak
