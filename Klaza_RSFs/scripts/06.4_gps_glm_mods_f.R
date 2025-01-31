#Script for gps Klaza glm models
## FALL RUT

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
library(ggplot2)
library(car)
library(AICcmodavg)

set.seed(123)

#Import used and available points with covariates previously extracted to them
f_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_f_data3.csv')

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
#Landcover
str(f_data$landcover2015_100m)
f_data$landcover2015_100m <- factor(f_data$landcover2015_100m)
levels(f_data$landcover2015_100m)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(f_data), 0.8*nrow(f_data))
f_train = f_data[n,]
f_test = f_data[-n,] 

######################################\
## FALL RUT
#Univariate Model Selection
#Landcover
lc30_1 <- glm(presence~landcover2015, family=binomial, data=f_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=f_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=f_train)
#LC 100 is the best based on AIC model selection. Other scales have closer delta AIC.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Lichen
lichbi30 <- glm(presence~lichen_bi, family=binomial, data=f_train)
lichenbi100 <- glm(presence~lichen_bi100, family=binomial, data=f_train)
lichenbi250 <- glm(presence~lichen_bi250, family=binomial, data=f_train)
#summary() # Lichen percent cover 100 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi, family=binomial, data=f_train)
conbi100 <- glm(presence~conifer_bi100, family=binomial, data=f_train)
conbi250 <- glm(presence~conifer_bi250, family=binomial, data=f_train)

#summary(conbi30) #conbi30 is best based on AIC scores. High, negative z score.
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi, family=binomial, data=f_train)
decidbi100 <- glm(presence~decid_bi100, family=binomial, data=f_train)
decidbi250 <- glm(presence~decid_bi250, family=binomial, data=f_train)

summary(decidbi250) #deciduous binomial 30 m is the best based on AIC model selection. decid bi 30 is very close second (delta AIC 0.43)
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~gram_bi, family=binomial, data=f_train)
grambi100 <- glm(presence~gram_bi100, family=binomial, data=f_train)
grambi250 <- glm(presence~gram_bi250, family=binomial, data=f_train)

summary(grambi30) #graminoid binomial 30 m is the best based on AIC scores. High, positive z score. Gram 30% is best percent cover covariate.
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)

#Inclue: dem (q), slope (q), LC100, Lichenbi 100m, Coniferbi 30m, Decidbi 30m, Grambi 30m
#excluding fire - probably not biologically relevant in this season

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Fall Rut:
#add quadratic terms
f_data$slope2 <- f_data$slope^2 #quadratic transformation
f_data$kdem2 <- f_data$kdem^2

#missing_cols <- setdiff(c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi", "conifer_bi100", "decid_bi", "gram_bi100", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000"), colnames(f_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.

f_dat.cor = f_data[,c("kdem", "kdem2", "slope", "slope2", "landcover2015_100m", "lichen_bi100", "conifer_bi", "decid_bi", "gram_bi", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")]
names(f_dat.cor) = c("kdem", "kdem2", "slope", "slope2", "landcover2015_100m", "lichen_bi100", "conifer_bi", "decid_bi", "gram_bi", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")
cor_matrix <- cor(f_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(f_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/f_correlation_matrix.csv", row.names = TRUE)

#############################################################################
##Fall Rut GLM Model Selection##

#Initial model exploration
#all variables
f_test1 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train)
#no elevation
f_test2 <- glm(presence~scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train)
#landcover only
f_test3 <- glm(presence~landcover2015_100m, family=binomial,data=f_train)

summary(f_test3)
#Neither slope term significant in full model
#deciduous forest and wetlands landcover class not significant
#graminoid cover not significant in full model
#Elevation and lichen are primary drivers of the model
#slope and graminoid significant when elevation is removed. Graminoid has a positive effect.


#Starting with GLM models first  --- UPDATE MODELS BASED ON COVARIATE EXPLORATION
#Unique model for each season
# Phase 1: Find the best environmental base model
#INCLUDE: dem (q), slope (q), LC100, Lichenbi 100m, Coniferbi 30m, Decidbi 30m, Grambi 30m
#exluding waterd and fire for this season
#Starting with glm (no mixed effect), standardized predictors:

#All top performing predictors from evaluations 
f_glm_1e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) 
#All top performing variables from evaluation, no graminoid
f_glm_2e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi), family=binomial,data=f_train) # all top performing variable,s no graminoid
#All top performing variables from evaluations, no deciduous shrub
f_glm_3e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(gram_bi), family=binomial,data=f_train) #All top performing variables from evaluations, no deciduous shrub
#All top performing variables from evaluations, no conifer
f_glm_4e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) #All top performing variables from evaluations, no conifer
#All top performing variables from evaluations, no lichen
f_glm_5e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) #All top performing variables from evaluations, no lichen                 
#All top performing variables from evaluations, no landcover
f_glm_6e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) #All top performing variables from evaluations, no landcover
#All top performing variables from evaluations, no slope
f_glm_7e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) #All top performing variables from evaluations, no slope
#all top performing variables from evaluations, no elevation
f_glm_8e <- glm(presence~scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(gram_bi), family=binomial,data=f_train) #all top performing variables from evaluations, no elevation
#null model
f_glm_9e <- glm(presence~1,family=binomial,data=f_train) #null model
#no graminoid, no slope
f_glm_10e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi), family=binomial,data=f_train)
#no graminoid, no conifer
f_glm_11e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_100m + scale(lichen_bi100) + scale(decid_bi), family=binomial,data=f_train)
#no graminoid, no conifer, no slope
f_glm_12e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(decid_bi), family=binomial,data=f_train)

f_e_glms <- list(f_glm_1e, f_glm_2e, f_glm_3e, f_glm_4e, f_glm_5e, f_glm_6e, f_glm_7e, f_glm_8e, f_glm_9e, f_glm_10e, f_glm_11e, f_glm_12e)
f_e_glms_names <- c('f_glm_1e', 'f_glm_2e', 'f_glm_3e', 'f_glm_4e', 'f_glm_5e', 'f_glm_6e', 'f_glm_7e', 'f_glm_8e', 'f_glm_9e', 'f_glm_10e', 'f_glm_11e', 'f_glm_12e')
aictab(cand.set = f_e_glms, modnames = f_e_glms_names)
#graminoid reduces model performance, not significant
#slope second least important predictor, not significant
# elevation, landcover, lichen are most important variables
#glm 10e without graminoid and slope had best fit based on delta AIC - simpler model and neither slope term is significant. 

summary(f_glm_10e) # Top performing model. 
summary(f_glm_8e)
summary()

report(lw_glm_5e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(f_glm_10e, newdata = f_test, type = "response")

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
predicted_probs_test <- predict(f_glm_10e, newdata = f_test, type = "response")

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

#Fair calibration. model has some deviation in the lower ranges. 
#the model tends to slightly under-predict

#Evaluation with AUC / ROC
#unscaled model:
f_glm_10e_us <- glm(presence~kdem + I(kdem^2) + landcover2015_100m + lichen_bi100 + conifer_bi + decid_bi, family=binomial,data=f_train)
#evaluate
f_glm_10e_eval <- pa_evaluate(predict(f_glm_10e_us, f_test[f_test$presence==1, ]), predict(f_glm_10e_us, f_test[f_test$presence==0, ]))
print(f_glm_10e_eval)
plot(f_glm_10e_eval, "ROC")
#AUC=0.971 #high AUC score
###########################################
#Model performance with easystats
r2(f_glm_10e_us) #Tjur's R2 = 0.638


#predicted occurrence map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif') %>%
  resample(mask, method="near")

#lw_rasters <- c(envrast, ndvirast, firerast, distdenrast, minebuffrast, rast100m, rast250m)
#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
landcover2015_100m <- as.factor(rast100m$landcover2015_100m) 
lichen_bi100 <- rast100m$lichen_bi100 |>
  resample(mask, method='near')
conifer_bi <- envrast$conifer_bi
decid_bi <- envrast$decid_bi

rasters <- c(kdem, landcover2015_100m, lichen_bi100, conifer_bi, decid_bi) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, f_glm_10e_us, type="response")
plot(p1) #plots!

#Add gps points to map
#filtered used points w/o covariates extracted to them:
f_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_f_daily.shp')

# Plot the map and add points on top
plot(st_geometry(f_points), add = TRUE, col = "red", pch = 19, cex = 0.6)


###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation. 
#distance to variables have a quadratic distribution in fall.
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 
#model 10e + disturbance 30m buff
f_glm_1d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_30), family=binomial,data=f_train)
#model 10e + disturbance 250m buff
f_glm_2d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_250), family=binomial,data=f_train)
#model 10e + disturbance 500m buff
f_glm_3d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_500), family=binomial,data=f_train)
#model 10e + disturbance 1000m buff
f_glm_4d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_1000), family=binomial,data=f_train)
#model 10e + disturbance 2000m buff
f_glm_5d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_2000), family=binomial,data=f_train)
#model 10e + disturbance 3000m buff
f_glm_6d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(dist2015_3000), family=binomial,data=f_train)
#model 10e + mining 30 m buff
f_glm_7d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_30), family=binomial,data=f_train)
#model 10e + mining 250 m buff
f_glm_8d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_250), family=binomial,data=f_train)
#model 10e + mining 500 m buff
f_glm_9d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_500), family=binomial,data=f_train)
#model 10e + mining 1000 m buff
f_glm_10d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_1000), family=binomial,data=f_train)
#model 10e + mining 2000 m buff
f_glm_11d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_2000), family=binomial,data=f_train)
#model 10e + mining 3000 m buff
f_glm_12d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(mine2015_3000), family=binomial,data=f_train)
#model 10e + roads and trails 30 m buff
f_glm_13d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(road2015_30), family=binomial,data=f_train)
#model 10e + roads and trails 250 m buff
f_glm_14d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(road2015_250), family=binomial,data=f_train)
#model 10e + roads and trails 500 m buff
f_glm_15d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(road2015_500), family=binomial,data=f_train)
#model 10e + roads and trails 1000 m buff
f_glm_16d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(road2015_1000), family=binomial,data=f_train)
#model 10e + disturbance density 250 m radii
f_glm_17d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(distden2015_250), family=binomial,data=f_train)
#model 10e + disturbance density 500 m radii
f_glm_18d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(distden2015_500), family=binomial,data=f_train)
#model 10e + disturbance density 1000 m radii
f_glm_19d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(distden2015_1000), family=binomial,data=f_train)
#model 10e + disturbance density 2000 m radii
f_glm_20d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(distden2015_2000), family=binomial,data=f_train)
#model 10e + disturbance density 3000 m radii
f_glm_21d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(distden2015_3000), family=binomial,data=f_train)
#model 10e + LF density 250 m radii
f_glm_22d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(lineden2015_250), family=binomial,data=f_train)
#model 10e + LF density 500 m radii
f_glm_23d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(lineden2015_500), family=binomial,data=f_train)
#model 10e + LF density 1000 m radii
f_glm_24d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(lineden2015_1000), family=binomial,data=f_train)
#model 10e + LF density 2000 m radii
f_glm_25d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(lineden2015_2000), family=binomial,data=f_train)
#model 10e + LF density 3000 m radii
f_glm_26d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi) + scale(lineden2015_3000), family=binomial,data=f_train)
#null model
f_glm_27d <- glm(presence~1,family=binomial,data=f_train)
#glm 10e
f_glm_10e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_100m + scale(lichen_bi100) + scale(conifer_bi) + scale(decid_bi), family=binomial,data=f_train)

#AIC
f_d_glms <- list(f_glm_1d, f_glm_2d, f_glm_3d, f_glm_4d, f_glm_5d, f_glm_6d, f_glm_7d, f_glm_8d, f_glm_9d, f_glm_10d, f_glm_11d, f_glm_12d, f_glm_13d, f_glm_14d, f_glm_15d, f_glm_16d, f_glm_17d, f_glm_18d, f_glm_19d, f_glm_20d, f_glm_21d, f_glm_22d, f_glm_23d, f_glm_24d, f_glm_25d, f_glm_26d, f_glm_27d, f_glm_10e)
f_d_glm_names <- c('f_glm_1d', 'f_glm_2d', 'f_glm_3d', 'f_glm_4d', 'f_glm_5d', 'f_glm_6d', 'f_glm_7d', 'f_glm_8d', 'f_glm_9d', 'f_glm_10d', 'f_glm_11d', 'f_glm_12d', 'f_glm_13d', 'f_glm_14d', 'f_glm_15d', 'f_glm_16d', 'f_glm_17d', 'f_glm_18d', 'f_glm_19d', 'f_glm_20d', 'f_glm_21d', 'f_glm_22d', 'f_glm_23d', 'f_glm_24d', 'f_glm_25d', 'f_glm_26d', 'f_glm_27d', 'f_glm_10e')
aictab(cand.set = f_d_glms, modnames = f_d_glm_names)

summary(f_glm_23d)
summary(f_glm_24d)
report()

#lineden 1000 is top variable
#dist 1000m buff is second top variable
#roads and trails 1000 m buff is third top variable
#buffered mining variables didn't perform well

#lineden 1000 is correlated with dist 1000 m buff and road 1000 m buff  (0.57 and 0.58). Originally ran in combined model but too much correlation in a model with so many variables. 
#dist 1000 m buff and road 1000 m buff are correlated

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(f_glm_24d, newdata = f_test, type = "response")

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
predicted_probs_test <- predict(f_glm_24d, newdata = f_test, type = "response")

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
#the model generallly tends to slightly under-predict

#Evaluation with AUC / ROC
#unscaled model:
f_glm_24d_us <- glm(presence~kdem + I(kdem^2) + landcover2015_100m + lichen_bi100 + conifer_bi + decid_bi + lineden2015_1000, family=binomial,data=f_train)
#evaluate
f_glm_24d_eval <- pa_evaluate(predict(f_glm_24d_us, f_test[f_test$presence==1, ]), predict(f_glm_24d_us, f_test[f_test$presence==0, ]))
print(f_glm_24d_eval)
plot(f_glm_24d_eval, "ROC")
#AUC=0.972 #high AUC score. Excellent performance and improved from base environmental model. 
###########################################
#Model performance with easystats
r2(f_glm_24d_us) #Tjur's R2 = 0.648

windows()
check_model(f_glm_24d_us)
#effect size plot with see pkg from easystats
plot(effectsize(f_glm_24d_us)) +
  ggplot2::labs(title = "Klaza - Fall Rut 2015") +
  ggplot2::scale_y_discrete(labels = c(
    "kdem" = "Elevation",
    "I(kdem^2)" = "Elevation^2",
    "landcover2015_100m5" = "Deciduous/mixed 100 m",
    "landcover2015_100m8" = "Shrublands 100 m",
    "landcover2015_100m10" = "Grasslands 100 m",
    "landcover2015_100m14" = "Wetlands 100 m",
    "landcover2015_100m16" = "Non-vegetated 100 m",
    "lichen_bi100" = "Lichen 100 m",
    "conifer_bi" = "Conifer 30 m",
    "decid_bi" = "Deciduous shrub 30 m",
    "lineden2015_1000" = "Line density 1000 m"
  )) 

plot(parameters(f_glm_24d_us)) +
  ggplot2::labs(title = "Fall Rut")

#Save models
saveRDS(f_glm_24d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_f_scaled.rds")
saveRDS(f_glm_24d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_f_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
linedenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/klinedenrast2015_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
landcover2015_100m <- as.factor(rast100m$landcover2015_100m) 
lichen_bi100 <- rast100m$lichen_bi100 |>
  resample(mask, method='near')
conifer_bi <- envrast$conifer_bi
decid_bi <- envrast$decid_bi |>
  resample(mask, method='near')
lineden2015_1000 <- linedenrast$lineden2015_1000

rasters <- c(kdem, landcover2015_100m, lichen_bi100, conifer_bi, decid_bi, lineden2015_1000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, f_glm_24d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_f_mod24d_prediction_3.tif')

#reload prediction raster
#f_prediction <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_f_mod24d_prediction.tif')
#Add gps points to map
#filtered used points w/o covariates extracted to them:
f_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_f_daily.shp')

# Plot the map and add points on top
#plot(f_prediction)
plot(st_geometry(f_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

####Calculate peak of quadratic relationship for elevation##

summary(f_glm_24d_us) #get unscaled coefficients for kdem & kdem2

f_dem_peak = -(0.04103 / (2 * (-1.304*10^-5)))
f_dem_peak # peak occurs at 1573.236 metres 
summary(f_train$kdem)
summary(f_test$kdem)

#Plot Relationship
# Coefficients from the model
beta_0 <- -35.58          # Intercept
beta_1 <- 0.04103         # Coefficient for kdem
beta_2 <- -1.304e-05      # Coefficient for I(kdem^2)

# Create a sequence of kdem values around the observed range
kdem_values <- seq(500, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * kdem_values + beta_2 * kdem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(kdem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Fall Rut")
abline(v = 1573.236, col = "red", lty = 2)  # Add vertical line at the peak
