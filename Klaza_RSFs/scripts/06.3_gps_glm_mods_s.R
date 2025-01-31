#Script for gps Klaza glm models
## SUMMER

#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

#Reviewed Nov 20, 2024
install.packages("easystats")
easystats::install_suggested()
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
s_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_s_data3.csv')

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(s_data$landcover2015_250m)
s_data$landcover2015_250m <- factor(s_data$landcover2015_250m)
levels(s_data$landcover2015_250m)

str(s_data$firedec)
s_data$firedec <- factor(s_data$firedec)
levels(s_data$firedec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(s_data), 0.8*nrow(s_data))
s_train = s_data[n,]
s_test = s_data[-n,] 

#########################################################################
#Univariate Model Selection

#Run univariate model selection on covariates that were unclear from above evaluations
#Landcover
lc30_1 <- glm(presence~landcover2015, family=binomial, data=s_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=s_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=s_train)

summary(lc100_1) #LC 250 is the best based on AIC model selection. large, positive z scores.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Lichen
lichbi30 <- glm(presence~lichen_bi, family=binomial, data=s_train)
lichenbi100 <- glm(presence~lichen_bi100, family=binomial, data=s_train)
lichenbi250 <- glm(presence~lichen_bi250, family=binomial, data=s_train)
summary(lichbi30) # Lichen binomial 250 m is the best based on AIC model selection. High, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi, family=binomial, data=s_train)
conbi100 <- glm(presence~conifer_bi100, family=binomial, data=s_train)
conbi250 <- glm(presence~conifer_bi250, family=binomial, data=s_train)
summary(conbi100) #conbi250 is best based on AIC scores. High, negative z score.
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi, family=binomial, data=s_train)
decidbi100 <- glm(presence~decid_bi100, family=binomial, data=s_train)
decidbi250 <- glm(presence~decid_bi250, family=binomial, data=s_train)
summary(decid250_1) #deciduous bi 100 m is the best based on AIC model selection. High, negative z score. 
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~gram_bi, family=binomial, data=s_train)
grambi100 <- glm(presence~gram_bi100, family=binomial, data=s_train)
grambi250 <- glm(presence~gram_bi250, family=binomial, data=s_train)
summary(grambi30) #graminoid binomial 30 m is the best based on AIC scores. High, positive z score. Gram 30% is best percent cover covariate.
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)

#Include: dem (q), slope (q), LC250, Lichenbi250, coniferbi250, decidbi100, grambi30
#excluding fire from summer models

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Summer:
#add quadratic terms
s_data$slope2 <- s_data$slope^2 #quadratic transformation
s_data$kdem2 <- s_data$kdem^2

#only doing 30 m versions and 250 m radii versions of disturbance variables
#missing_cols <- setdiff(c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi", "conifer_bi100", "decid_bi", "gram_bi100", dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000"), colnames(s_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.


s_dat.cor = s_data[,c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi250", "conifer_bi250", "decid_bi100", "gram_bi", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")]
names(s_dat.cor) = c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi250", "conifer_bi250", "decid_bi100", "gram_bi", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")
cor_matrix <- cor(s_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(s_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/s_correlation_matrix.csv", row.names = TRUE)

##############################################################################
## SUMMER GLM MODEL SELECTION

#Inital model exploration
#all variables
s_test1 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train)
#elevation and landcover
s_test2 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_250m, family=binomial,data=s_train)
#no elevation
s_test3 <- glm(presence~scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train)

summary(s_test3)
#Graminoid is not significant and has very small positive effect in full model
#deciduous/mixed wood forest and wetlands landcover classes are not significant in full model
#coefficients of environmental variables make sense and are in agreement
#elevation and landcover are primary drivers of model
#graminoid is significant (positive) when elevation is removed. Some correlation (0.4). Graminoid has smaller effect size, will likely be removed based on AIC.


#Starting with GLM models first  --- UPDATE MODELS BASED ON COVARIATE EXPLORATION
#Unique model for each season
# Phase 1: Find the best environmental base model
#INCLUDE: #Include: dem (q), slope (q), LC250, Lichenbi250, coniferbi250, decidbi100, grambi30
#exluding fire and waterd for summer season
#Starting with glm (no mixed effect), standardized predictors:

#All top performing predictors from evaluations
s_glm_1e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #All top performing predictors from evaluations
#All top performing variables from evaluation, no graminoid
s_glm_2e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100), family=binomial,data=s_train) # all top performing variable,s no graminoid
#All top performing variables from evaluations, no deciduous shrub
s_glm_3e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(gram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no deciduous shrub
#All top performing variables from evaluations, no conifer
s_glm_4e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no conifer
#All top performing variables from evaluations, no lichen
s_glm_5e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no lichen                 
#All top performing variables from evaluations, no landcover
s_glm_6e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no landcover
#All top performing variables from evaluations, no slope
s_glm_7e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #All top performing variables from evaluations, no slope
#all top performing variables from evaluations, no elevation
s_glm_8e <- glm(presence~scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train) #all top performing variables from evaluations, no elevation
#null model
s_glm_9e <- glm(presence~1,family=binomial,data=s_train) #null model
#no slope, no graminoid
#s_glm_11e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(gram_bi), family=binomial,data=s_train)

s_e_glms <- list(s_glm_1e, s_glm_2e, s_glm_3e, s_glm_4e, s_glm_5e, s_glm_6e, s_glm_7e, s_glm_8e, s_glm_9e)
s_e_glms_names <- c('s_glm_1e', 's_glm_2e', 's_glm_3e', 's_glm_4e', 's_glm_5e', 's_glm_6e', 's_glm_7e', 's_glm_8e', 's_glm_9e')
aictab(cand.set = s_e_glms, modnames = s_e_glms_names)

#top performing model doesn't have graminoid cover. Graminoid wasn't significant. 
#elevation, lichen, landcover, conifer all important variables 
#slope not very important
summary(s_glm_2e)
summary(s_glm_1e)
summary()

report()

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(s_glm_2e, newdata = s_test, type = "response")

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
predicted_probs_test <- predict(s_glm_2e, newdata = s_test, type = "response")

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

#Evaluation with AUC / ROC
#unscaled model:
s_glm_2e_us <- glm(presence~kdem + I(kdem^2) + slope + I(slope^2) + landcover2015_250m + lichen_bi250 + conifer_bi250 + decid_bi100, family=binomial,data=s_train)
#evaluate
s_glm_2e_eval <- pa_evaluate(predict(s_glm_2e_us, s_test[s_test$presence==1, ]), predict(s_glm_2e_us, s_test[s_test$presence==0, ]))
print(s_glm_2e_eval)
plot(s_glm_2e_eval, "ROC")
#AUC=0.953 #high AUC score, excellent model performance.
###########################################
#Model performance with easystats
r2(s_glm_2e_us) #Tjur's R2 = 0.586

#predicted occurrence map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif') %>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
landcover2015_250m <- as.factor(rast250m$landcover2015_250m) |>
  resample(mask, method='near')
slope <- envrast$slope
lichen_bi250 <- rast250m$lichen_bi250 |>
  resample(mask, method='near')
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
decid_bi100 <- rast100m$decid_bi100 |>
  resample(mask, method='near')

rasters <- c(kdem, slope, landcover2015_250m, lichen_bi250, conifer_bi250, decid_bi100) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, s_glm_2e_us, type="response")
plot(p1) #plots!

#add gps points
#filtered used points w/o covariates extracted to them:
s_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_s_daily.shp')

#Plot the map and add points on top
plot(st_geometry(s_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation. 
#distance to variables have a quadratic distribution in summer.
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 
#model 2e + disturbance 30m buff
s_glm_1d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_30), family=binomial,data=s_train)
#model 2e + disturbance 250m buff
s_glm_2d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_250), family=binomial,data=s_train)
#model 2e + disturbance 500m buff
s_glm_3d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_500), family=binomial,data=s_train)
#model 2e + disturbance 1000m buff
s_glm_4d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_1000), family=binomial,data=s_train)
#model 2e + disturbance 2000m buff
s_glm_5d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_2000), family=binomial,data=s_train)
#model 2e + disturbance 3000m buff
s_glm_5d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(dist2015_3000), family=binomial,data=s_train)
#model 2e + mining 30 m buff
s_glm_6d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_30), family=binomial,data=s_train)
#model 2e + mining 250 m buff
s_glm_7d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_250), family=binomial,data=s_train)
#model 2e + mining 500 m buff
s_glm_8d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_500), family=binomial,data=s_train)
#model 2e + mining 1000 m buff
s_glm_9d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_1000), family=binomial,data=s_train)
#model 2e + mining 2000 m buff
s_glm_10d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_2000), family=binomial,data=s_train)
#model 2e + mining 3000 m buff
s_glm_11d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(mine2015_3000), family=binomial,data=s_train)
#model 2e + roads and trails 30 m buff
s_glm_12d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(road2015_30), family=binomial,data=s_train)
#model 2e + roads and trails 250 m buff
s_glm_13d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(road2015_250), family=binomial,data=s_train)
#model 2e + roads and trails 500 m buff
s_glm_14d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(road2015_500), family=binomial,data=s_train)
#model 2e + roads and trails 1000 m buff
s_glm_15d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(road2015_1000), family=binomial,data=s_train)
#model 2e + disturbance density 250 m radii
s_glm_16d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_250), family=binomial,data=s_train)
#model 2e + disturbance density 500 m radii
s_glm_17d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_500), family=binomial,data=s_train)
#model 2e + disturbance density 1000 m radii
s_glm_18d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_1000), family=binomial,data=s_train)
#model 2e + disturbance density 2000 m radii
s_glm_19d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_2000), family=binomial,data=s_train)
#model 2e + disturbance density 3000 m radii
s_glm_20d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_3000), family=binomial,data=s_train)
#model 2e + LF density 250 m radii
s_glm_21d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(lineden2015_250), family=binomial,data=s_train)
#model 2e + LF density 500 m radii
s_glm_22d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(lineden2015_500), family=binomial,data=s_train)
#model 2e + LF density 1000 m radii
s_glm_23d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(lineden2015_1000), family=binomial,data=s_train)
#model 2e + LF density 2000 m radii
s_glm_24d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(lineden2015_2000), family=binomial,data=s_train)
#model 2e + LF density 3000 m radii
s_glm_25d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(lineden2015_3000), family=binomial,data=s_train)
#model 2e
s_glm_2e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100), family=binomial,data=s_train)
#null model
s_glm_26d <- glm(presence~1,family=binomial,data=s_train)
#model 2e + distden 1000 + mining 500 m buff
s_glm_27d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi250) + scale(conifer_bi250) + scale(decid_bi100) + scale(distden2015_1000) + scale(mine2015_500), family=binomial,data=s_train)

#AIC
s_d_glms <- list(s_glm_1d, s_glm_2d, s_glm_3d, s_glm_4d, s_glm_5d, s_glm_6d, s_glm_7d, s_glm_8d, s_glm_9d, s_glm_10d, s_glm_11d, s_glm_12d, s_glm_13d, s_glm_14d, s_glm_15d, s_glm_16d, s_glm_17d, s_glm_18d, s_glm_19d, s_glm_20d, s_glm_21d, s_glm_22d, s_glm_23d, s_glm_24d, s_glm_25d, s_glm_26d, s_glm_27d, s_glm_2e)
s_d_glm_names <- c('s_glm_1d', 's_glm_2d', 's_glm_3d', 's_glm_4d', 's_glm_5d', 's_glm_6d', 's_glm_7d', 's_glm_8d', 's_glm_9d', 's_glm_10d', 's_glm_11d', 's_glm_12d', 's_glm_13d', 's_glm_14d', 's_glm_15d', 's_glm_16d', 's_glm_17d', 's_glm_18d', 's_glm_19d', 's_glm_20d', 's_glm_21d', 's_glm_22d', 's_glm_23d', 's_glm_24d', 's_glm_25d', 's_glm_26d', 's_glm_27d', 's_glm_2e')
aictab(cand.set = s_d_glms, modnames = s_d_glm_names)

summary(s_glm_18d)
summary(s_glm_27d)
report()

#distden 1000 m is top predictor
#lineden 250 is second top predictor, followed by lineden 3000 
#mine 500 m buffer is top buffered mining variable but not significant
#buffered roads and trails didn't perform very well
#top variables are correlated, will try distden 1000 and road 500 m buff in same model


## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(s_glm_27d, newdata = s_test, type = "response")

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
predicted_probs_test <- predict(s_glm_27d, newdata = s_test, type = "response")

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

#Good calibration. model performs well. Predicted probabilities align very closely with the observed rates across the 10 bins. 

#Evaluation with AUC / ROC
#unscaled model:
s_glm_27d_us <- glm(presence~kdem + I(kdem^2) + slope + I(slope^2) + landcover2015_250m + lichen_bi250 + conifer_bi250 + decid_bi100 + distden2015_1000 + mine2015_500, family=binomial,data=s_train)
#evaluate
s_glm_27d_eval <- pa_evaluate(predict(s_glm_27d_us, s_test[s_test$presence==1, ]), predict(s_glm_27d_us, s_test[s_test$presence==0, ]))
print(s_glm_27d_eval)
plot(s_glm_27d_eval, "ROC")
#AUC=0.955 #Excellent performance. Improved from base environmental model. 

###########################################
#Model performance and plots with easystats
r2(s_glm_27d_us) #Tjur's R2 = 0.591

windows()
check_model(s_glm_27d_us)
#effect size plot with see pkg from easystats
plot(effectsize(s_glm_27d_us)) +
  ggplot2::labs(title = "Klaza - Summer 2015") +
  ggplot2::scale_y_discrete(labels = c(
    "kdem" = "Elevation",
    "I(kdem^2)" = "Elevation^2",
    "slope" = "Slope",
    "I(slope^2)" = "Slope^2",
    "landcover2015_250m5" = "Deciduous/mixed 250 m",
    "landcover2015_250m8" = "Shrublands 250 m",
    "landcover2015_250m10" = "Grasslands 250 m",
    "landcover2015_250m14" = "Wetlands 250 m",
    "landcover2015_250m16" = "Non-vegetated 250 m",
    "lichen_bi250" = "Lichen 250 m",
    "conifer_bi250" = "Conifer 250 m",
    "decid_bi100" = "Deciduous shrub 100 m",
    "distden2015_1000" = "Disturbance density 1000 m",
    "mine2015_500" = "Mining 500 m"
  )) 


plot(parameters(s_glm_27d_us)) +
  ggplot2::labs(title = "Summer")

#Save models
saveRDS(s_glm_27d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_s_scaled.rds")
saveRDS(s_glm_27d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_s_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kminebuffrast2015_1.tif')
distdenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kdistdenrast2015_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_100m2.tif') %>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
landcover2015_250m <- as.factor(rast250m$landcover2015_250m) |>
  resample(mask, method='near')
slope <- envrast$slope
lichen_bi250 <- rast250m$lichen_bi250 |>
  resample(mask, method='near')
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
decid_bi100 <- rast100m$decid_bi100 |>
  resample(mask, method='near')
distden2015_1000 <- distdenrast$distden2015_1000
mine2015_500 <- minebuffrast$mine2015_500

rasters <- c(kdem, slope, landcover2015_250m, lichen_bi250, conifer_bi250, decid_bi100, distden2015_1000, mine2015_500) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, s_glm_27d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_s_mod27d_prediction_3.tif')

#reload prediction map
#s_prediction <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_s_mod27d_prediction.tif')
#Add gps points to map
#filtered used points w/o covariates extracted to them:
s_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_s_daily.shp')

# Plot the map and add points on top
plot(st_geometry(s_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

####Calculate peak of quadratic relationship for elevation##

summary(s_glm_27d_us) #get unscaled coefficients for kdem & kdem2

s_dem_peak = -(0.03083 / (2 * (-9.508*10^-6)))
s_dem_peak # peak occurs at 1621.266 metres 
summary(s_train$kdem) #max elevation in dataset is 1941

#Plot Relationship
# Coefficients from the model
beta_0 <- -27.16          # Intercept
beta_1 <- 0.03083         # Coefficient for kdem
beta_2 <- -9.508e-06      # Coefficient for I(kdem^2)

# Create a sequence of kdem values around the observed range
kdem_values <- seq(500, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * kdem_values + beta_2 * kdem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(kdem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Summer")
abline(v = 1621.266, col = "red", lty = 2)  # Add vertical line at the peak
